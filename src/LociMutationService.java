import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Central mutation coordinator for locus operations (split/merge/resize/undo).
 * All mutations are synchronized and produce a standardized MutationResult.
 */
public class LociMutationService {

    public static class MutationResult {
        public boolean ok;
        public String error;
        public String mutationId;
        public String mutationType;
        public Map<String, List<String>> retiredToNewIds = new LinkedHashMap<>();
        public List<String> affectedIds = new ArrayList<>();
        public String manifestJson;
    }

    static class JournalEntry {
        String mutationId;
        String type;
        long timestamp;
        Map<String, List<String>> retiredToNewIds;
        List<String> backupFiles = new ArrayList<>();
        String backupDir;
    }

    private final Config config;
    private final GffParser gff;
    private final List<Locus> loci;
    private final List<LocusOutput> outputs;
    private final Deque<JournalEntry> undoJournal = new ArrayDeque<>();
    private static final int MAX_UNDO = 20;

    public LociMutationService(Config config, GffParser gff,
                                List<Locus> loci, List<LocusOutput> outputs) {
        this.config = config;
        this.gff = gff;
        this.loci = loci;
        this.outputs = outputs;
    }

    public synchronized MutationResult resize(int locusIndex, long newStart, long newEnd) {
        MutationResult mr = new MutationResult();
        mr.mutationId = UUID.randomUUID().toString();
        mr.mutationType = "resize";
        try {
            Locus original = findByIndex(locusIndex);
            if (original == null) { mr.error = "Locus " + locusIndex + " not found"; return mr; }

            snapshotForUndo(mr.mutationId, "resize", Arrays.asList(locusIndex));

            LocusUpdater.UpdateResult ur = LocusUpdater.update(
                locusIndex, newStart, newEnd, config, gff, loci);
            if (!ur.ok) { mr.error = ur.error; return mr; }

            mr.affectedIds.add(original.id);
            mr.ok = true;
            mr.manifestJson = readManifestJson();
            addJournalEntry(mr);
        } catch (Exception e) {
            mr.error = e.getMessage();
            e.printStackTrace();
        }
        return mr;
    }

    public synchronized MutationResult split(int locusIndex, List<LocusUpdater.SplitRegion> regions) {
        MutationResult mr = new MutationResult();
        mr.mutationId = UUID.randomUUID().toString();
        mr.mutationType = "split";
        try {
            Locus original = findByIndex(locusIndex);
            if (original == null) { mr.error = "Locus " + locusIndex + " not found"; return mr; }

            snapshotForUndo(mr.mutationId, "split", allIndices());

            String retiredId = original.id;
            LocusUpdater.UpdateResult ur = LocusUpdater.split(
                locusIndex, regions, config, gff, loci, outputs);
            if (!ur.ok) { mr.error = ur.error; return mr; }

            List<String> newIds = new ArrayList<>();
            for (LocusOutput lo : outputs) {
                if (lo.id != null && !lo.id.equals(retiredId)) {
                    newIds.add(lo.id);
                }
            }
            mr.retiredToNewIds.put(retiredId, newIds);
            mr.affectedIds.addAll(newIds);

            reExportManifest();
            mr.ok = true;
            mr.manifestJson = readManifestJson();
            addJournalEntry(mr);
        } catch (Exception e) {
            mr.error = e.getMessage();
            e.printStackTrace();
        }
        return mr;
    }

    public synchronized MutationResult delete(int locusIndex) {
        MutationResult mr = new MutationResult();
        mr.mutationId = UUID.randomUUID().toString();
        mr.mutationType = "delete";
        try {
            Locus target = findByIndex(locusIndex);
            if (target == null) { mr.error = "Locus " + locusIndex + " not found"; return mr; }

            snapshotForUndo(mr.mutationId, "delete", allIndices());

            String retiredId = target.id;
            mr.retiredToNewIds.put(retiredId, Collections.emptyList());
            mr.affectedIds.add(retiredId);

            loci.removeIf(l -> l.index == locusIndex);
            outputs.removeIf(o -> o.locusIndex == locusIndex);

            String dataDir = config.outputDir + "/data";
            new File(dataDir + "/locus_" + locusIndex + ".json").delete();
            new File(dataDir + "/locus_" + locusIndex + ".js").delete();

            reExportManifest();
            mr.ok = true;
            mr.manifestJson = readManifestJson();
            addJournalEntry(mr);

            System.out.printf("[LociMutationService] Deleted locus %d (id=%s)%n", locusIndex, retiredId);

        } catch (Exception e) {
            mr.error = e.getMessage();
            e.printStackTrace();
        }
        return mr;
    }

    public synchronized MutationResult merge(List<Integer> locusIndices, String mergedName) {
        MutationResult mr = new MutationResult();
        mr.mutationId = UUID.randomUUID().toString();
        mr.mutationType = "merge";
        try {
            if (locusIndices == null || locusIndices.size() < 2) {
                mr.error = "Need at least 2 loci to merge"; return mr;
            }

            List<Locus> toMerge = new ArrayList<>();
            for (int idx : locusIndices) {
                Locus l = findByIndex(idx);
                if (l == null) { mr.error = "Locus " + idx + " not found"; return mr; }
                toMerge.add(l);
            }

            // Validate: all same chromosome
            Set<String> chrs = new HashSet<>();
            for (Locus l : toMerge) chrs.add(l.chr);
            if (chrs.size() > 1) { mr.error = "Cannot merge loci on different chromosomes"; return mr; }

            snapshotForUndo(mr.mutationId, "merge", allIndices());

            // Compute merged boundaries
            String chr = toMerge.get(0).chr;
            long mergedStart = Long.MAX_VALUE, mergedEnd = Long.MIN_VALUE;
            for (Locus l : toMerge) {
                mergedStart = Math.min(mergedStart, l.start);
                mergedEnd = Math.max(mergedEnd, l.end);
            }

            // Track retired IDs
            List<String> retiredIds = new ArrayList<>();
            for (Locus l : toMerge) retiredIds.add(l.id);

            // Remove originals
            Set<Integer> removeIndices = new HashSet<>(locusIndices);
            loci.removeIf(l -> removeIndices.contains(l.index));
            outputs.removeIf(o -> removeIndices.contains(o.locusIndex));

            // Delete old JSON files
            String dataDir = config.outputDir + "/data";
            for (int idx : locusIndices) {
                new File(dataDir + "/locus_" + idx + ".json").delete();
                new File(dataDir + "/locus_" + idx + ".js").delete();
            }

            // Create merged locus
            if (mergedName == null || mergedName.isEmpty())
                mergedName = "Merged Locus";

            LocusUpdater.UpdateResult ur = LocusUpdater.create(
                chr, mergedStart, mergedEnd, mergedName, config, gff, loci, outputs);
            if (!ur.ok) { mr.error = "Failed creating merged locus: " + ur.error; return mr; }

            // Find the new locus ID
            LocusOutput newOutput = outputs.get(outputs.size() - 1);
            String newId = newOutput.id;

            for (String retiredId : retiredIds) {
                mr.retiredToNewIds.put(retiredId, Collections.singletonList(newId));
            }
            mr.affectedIds.add(newId);

            reExportManifest();
            mr.ok = true;
            mr.manifestJson = readManifestJson();
            addJournalEntry(mr);

            System.out.printf("[LociMutationService] Merged %d loci into %s: chr%s:%d-%d%n",
                locusIndices.size(), mergedName, chr, mergedStart, mergedEnd);

        } catch (Exception e) {
            mr.error = e.getMessage();
            e.printStackTrace();
        }
        return mr;
    }

    public synchronized MutationResult undo() {
        MutationResult mr = new MutationResult();
        mr.mutationId = UUID.randomUUID().toString();
        mr.mutationType = "undo";
        try {
            if (undoJournal.isEmpty()) {
                mr.error = "Nothing to undo"; return mr;
            }
            JournalEntry entry = undoJournal.pop();

            // Restore backed-up files
            if (entry.backupDir != null) {
                File backupDir = new File(entry.backupDir);
                File dataDir = new File(config.outputDir + "/data");
                if (backupDir.isDirectory()) {
                    // Clear current data dir of locus files and manifest
                    for (File f : dataDir.listFiles()) {
                        if (f.getName().startsWith("locus_") || f.getName().startsWith("manifest.")) {
                            f.delete();
                        }
                    }
                    // Copy backup files to data dir
                    for (File f : backupDir.listFiles()) {
                        Files.copy(f.toPath(), new File(dataDir, f.getName()).toPath(),
                            StandardCopyOption.REPLACE_EXISTING);
                    }
                    // Clean up backup
                    deleteDirectory(backupDir);
                }
            }

            // Reload state from disk
            reloadStateFromDisk();

            mr.ok = true;
            mr.manifestJson = readManifestJson();
            System.out.printf("[LociMutationService] Undo '%s' (mutation %s) successful%n",
                entry.type, entry.mutationId);

        } catch (Exception e) {
            mr.error = e.getMessage();
            e.printStackTrace();
        }
        return mr;
    }

    public boolean hasUndo() {
        return !undoJournal.isEmpty();
    }

    // ── Internal helpers ───────────────────────────────────────────

    private Locus findByIndex(int index) {
        for (Locus l : loci) if (l.index == index) return l;
        return null;
    }

    private Locus findById(String id) {
        for (Locus l : loci) if (l.id.equals(id)) return l;
        return null;
    }

    private List<Integer> allIndices() {
        List<Integer> indices = new ArrayList<>();
        for (Locus l : loci) indices.add(l.index);
        return indices;
    }

    private void snapshotForUndo(String mutationId, String type, List<Integer> affectedIndices) {
        try {
            String backupDir = config.outputDir + "/undo/" + mutationId;
            new File(backupDir).mkdirs();
            String dataDir = config.outputDir + "/data";

            // Back up manifest
            copyIfExists(dataDir + "/manifest.json", backupDir + "/manifest.json");
            copyIfExists(dataDir + "/manifest.js", backupDir + "/manifest.js");

            // Back up affected locus files
            for (int idx : affectedIndices) {
                copyIfExists(dataDir + "/locus_" + idx + ".json", backupDir + "/locus_" + idx + ".json");
                copyIfExists(dataDir + "/locus_" + idx + ".js", backupDir + "/locus_" + idx + ".js");
            }

            JournalEntry entry = new JournalEntry();
            entry.mutationId = mutationId;
            entry.type = type;
            entry.timestamp = System.currentTimeMillis();
            entry.backupDir = backupDir;
            undoJournal.push(entry);

            while (undoJournal.size() > MAX_UNDO) {
                JournalEntry old = undoJournal.removeLast();
                if (old.backupDir != null) deleteDirectory(new File(old.backupDir));
            }
        } catch (Exception e) {
            System.err.println("[LociMutationService] Snapshot failed: " + e.getMessage());
        }
    }

    private void addJournalEntry(MutationResult mr) {
        // Journal entry already created in snapshotForUndo; update with ID mappings
        if (!undoJournal.isEmpty()) {
            JournalEntry entry = undoJournal.peek();
            if (entry.mutationId.equals(mr.mutationId)) {
                entry.retiredToNewIds = mr.retiredToNewIds;
            }
        }
    }

    private void reExportManifest() throws IOException {
        JsonExporter.exportManifest(outputs, config);
    }

    private String readManifestJson() throws IOException {
        File mf = new File(config.outputDir + "/data/manifest.json");
        if (mf.exists()) return new String(Files.readAllBytes(mf.toPath()), "UTF-8");
        return "{}";
    }

    private void reloadStateFromDisk() throws IOException {
        // Re-read manifest to rebuild loci/outputs lists
        // This is a simplified reload — re-parse manifest and reload each locus JSON
        String dataDir = config.outputDir + "/data";
        File manifestFile = new File(dataDir + "/manifest.json");
        if (!manifestFile.exists()) return;

        String manifestStr = new String(Files.readAllBytes(manifestFile.toPath()), "UTF-8");

        // Clear and rebuild from files on disk
        loci.clear();
        outputs.clear();

        // Parse manifest to find locus indices
        List<Integer> indices = new ArrayList<>();
        int lociArr = manifestStr.indexOf("\"loci\"");
        if (lociArr >= 0) {
            int start = manifestStr.indexOf('[', lociArr);
            int end = manifestStr.lastIndexOf(']');
            if (start >= 0 && end > start) {
                String arrContent = manifestStr.substring(start + 1, end);
                int pos = 0;
                while (pos < arrContent.length()) {
                    int idxKey = arrContent.indexOf("\"index\":", pos);
                    if (idxKey < 0) break;
                    int valStart = idxKey + 8;
                    while (valStart < arrContent.length() && arrContent.charAt(valStart) == ' ') valStart++;
                    int valEnd = valStart;
                    while (valEnd < arrContent.length() && Character.isDigit(arrContent.charAt(valEnd))) valEnd++;
                    if (valEnd > valStart) {
                        indices.add(Integer.parseInt(arrContent.substring(valStart, valEnd)));
                    }
                    pos = valEnd;
                }
            }
        }

        // For each index, try to load the locus JSON and reconstruct minimal state
        for (int idx : indices) {
            File locusFile = new File(dataDir + "/locus_" + idx + ".json");
            if (!locusFile.exists()) continue;
            // We only need to reconstruct enough state for subsequent mutations
            // The full state is read from the JSON files by the frontend
            String json = new String(Files.readAllBytes(locusFile.toPath()), "UTF-8");
            String chr = extractSimpleField(json, "chr");
            String startStr = extractSimpleField(json, "start");
            String endStr = extractSimpleField(json, "end");
            String id = extractSimpleField(json, "id");
            if (chr != null && startStr != null && endStr != null) {
                Locus l = id != null
                    ? new Locus(id, idx, chr, Long.parseLong(startStr), Long.parseLong(endStr), config.locusPadding)
                    : new Locus(idx, chr, Long.parseLong(startStr), Long.parseLong(endStr), config.locusPadding);
                loci.add(l);

                LocusOutput lo = new LocusOutput();
                lo.id = l.id;
                lo.locusIndex = idx;
                lo.locusName = "Locus " + idx;
                lo.chr = chr;
                lo.start = l.start;
                lo.end = l.end;
                lo.paddedStart = l.paddedStart;
                lo.paddedEnd = l.paddedEnd;

                // Parse top_snp from JSON
                int topSnpIdx = json.indexOf("\"top_snp\":");
                if (topSnpIdx >= 0) {
                    String topId = extractNestedField(json, topSnpIdx, "id");
                    String topPos = extractNestedField(json, topSnpIdx, "pos");
                    String topP = extractNestedField(json, topSnpIdx, "pvalue");
                    if (topId != null && topPos != null && topP != null) {
                        try {
                            Snp top = new Snp(topId, chr, Long.parseLong(topPos),
                                Double.parseDouble(topP), "", "");
                            lo.topSnp = top;
                        } catch (NumberFormatException ignore) {}
                    }
                }

                // Parse nearest_genes
                lo.nearestGenes = extractStringArray(json, "nearest_genes");
                lo.genes = gff.overlapping(chr, l.paddedStart, l.paddedEnd);

                outputs.add(lo);
            }
        }

        loci.sort(Comparator.comparingInt(Locus::chrInt).thenComparingLong(l -> l.start));
        outputs.sort(Comparator.comparingInt((LocusOutput o) -> Locus.chrToInt(o.chr))
            .thenComparingLong(o -> o.start));
    }

    private static void copyIfExists(String src, String dest) {
        try {
            File f = new File(src);
            if (f.exists()) Files.copy(f.toPath(), new File(dest).toPath(),
                StandardCopyOption.REPLACE_EXISTING);
        } catch (IOException e) {
            System.err.println("[LociMutationService] Copy failed: " + src + " → " + e.getMessage());
        }
    }

    private static void deleteDirectory(File dir) {
        File[] files = dir.listFiles();
        if (files != null) {
            for (File f : files) {
                if (f.isDirectory()) deleteDirectory(f);
                else f.delete();
            }
        }
        dir.delete();
    }

    private static String extractSimpleField(String json, String key) {
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return null;
        i += marker.length();
        while (i < json.length() && json.charAt(i) == ' ') i++;
        if (i >= json.length()) return null;
        if (json.charAt(i) == '"') {
            StringBuilder sb = new StringBuilder(); i++;
            while (i < json.length()) {
                char c = json.charAt(i++);
                if (c == '"') break;
                if (c == '\\' && i < json.length()) { sb.append(json.charAt(i++)); continue; }
                sb.append(c);
            }
            return sb.toString();
        }
        if (json.startsWith("null", i)) return null;
        int end = i;
        while (end < json.length() && ",}]".indexOf(json.charAt(end)) < 0) end++;
        return json.substring(i, end).trim();
    }

    private static String extractNestedField(String json, int startFrom, String key) {
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker, startFrom);
        if (i < 0 || i > startFrom + 500) return null;
        return extractSimpleField(json.substring(i - marker.length() + marker.length()), key);
    }

    private static List<String> extractStringArray(String json, String key) {
        List<String> result = new ArrayList<>();
        String marker = "\"" + key + "\":";
        int i = json.indexOf(marker);
        if (i < 0) return result;
        i += marker.length();
        while (i < json.length() && json.charAt(i) != '[') i++;
        if (i >= json.length()) return result;
        int end = json.indexOf(']', i);
        if (end < 0) return result;
        String arrContent = json.substring(i + 1, end);
        for (String part : arrContent.split(",")) {
            part = part.trim();
            if (part.startsWith("\"") && part.endsWith("\"")) {
                result.add(part.substring(1, part.length() - 1));
            }
        }
        return result;
    }
}
