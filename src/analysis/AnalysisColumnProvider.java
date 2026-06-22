
import java.io.*;
import java.util.*;
import export.ColumnSpec;
import export.SnpColumnProvider;
import export.SnpContext;

/**
 * Exposes analysis results (from any tool's unified result.tsv) as columns
 * in the Excel export via the existing SnpColumnProvider registry.
 *
 * Scans the locus's runs/ directory for the latest successful result,
 * reads result.manifest.json to discover columns, and joins result.tsv
 * to the export by stable SNP id.
 */
public class AnalysisColumnProvider implements SnpColumnProvider {

    private final String projectDir;
    private Map<String, Map<String, String>> snpResults;
    private List<ColumnSpec> columnSpecs;
    private boolean loaded = false;
    private String lastLocusId = null;

    public AnalysisColumnProvider(String projectDir) {
        this.projectDir = projectDir;
    }

    @Override
    public List<ColumnSpec> columns() {
        if (columnSpecs == null) return Collections.emptyList();
        return columnSpecs;
    }

    @Override
    public Object value(SnpContext ctx, ColumnSpec col) {
        ensureLoaded(ctx);
        if (snpResults == null) return null;

        // Build stable SNP ID from context
        String stableId = StableSnpId.of(ctx.chr, ctx.pos, ctx.ea, ctx.nea).id;

        Map<String, String> row = snpResults.get(stableId);
        if (row == null) return null;

        String val = row.get(col.key);
        if (val == null || val.equals("NA")) return null;

        try {
            switch (col.type) {
                case DOUBLE:
                case SCIENTIFIC: return Double.parseDouble(val);
                case INT: return (int) Double.parseDouble(val);
                case BOOL: return Boolean.parseBoolean(val);
                default: return val;
            }
        } catch (NumberFormatException e) {
            return val;
        }
    }

    @Override
    public boolean isAvailable(String projectDir) {
        File analysisDir = new File(projectDir, "loci_analysis");
        return analysisDir.isDirectory();
    }

    private void ensureLoaded(SnpContext ctx) {
        // Reload when locus changes
        String locusId = findLocusId(ctx);
        if (locusId == null) return;
        if (loaded && locusId.equals(lastLocusId)) return;

        lastLocusId = locusId;
        loaded = true;
        snpResults = new LinkedHashMap<>();
        columnSpecs = new ArrayList<>();

        File runsDir = new File(projectDir, "loci_analysis/" + locusId + "/runs");
        if (!runsDir.isDirectory()) return;

        // Find latest run with a valid result
        File[] runDirs = runsDir.listFiles(File::isDirectory);
        if (runDirs == null || runDirs.length == 0) return;
        Arrays.sort(runDirs, Comparator.comparingLong(File::lastModified).reversed());

        for (File rd : runDirs) {
            File resultFile = new File(rd, "result.tsv");
            File manifestFile = new File(rd, "result.manifest.json");
            if (!resultFile.exists() || !manifestFile.exists()) continue;

            try {
                // Parse manifest for column declarations
                OutputContractValidator.ValidationResult vr = OutputContractValidator.validate(rd);
                if (!vr.valid) continue;

                int order = 900; // analysis columns come after standard columns
                for (OutputContractValidator.ColumnDecl cd : vr.columns) {
                    ColumnSpec.ColType ct = ColumnSpec.ColType.DOUBLE;
                    if ("int".equals(cd.type)) ct = ColumnSpec.ColType.INT;
                    else if ("string".equals(cd.type)) ct = ColumnSpec.ColType.TEXT;
                    else if ("boolean".equals(cd.type)) ct = ColumnSpec.ColType.BOOL;
                    String group = cd.method != null ? cd.method.toUpperCase() : "Analysis";
                    columnSpecs.add(new ColumnSpec(cd.name, cd.name, group, ct, order++));
                }

                // Load result data
                try (BufferedReader br = new BufferedReader(new FileReader(resultFile))) {
                    String header = br.readLine();
                    if (header == null) continue;
                    String[] headers = header.split("\t", -1);

                    String line;
                    while ((line = br.readLine()) != null) {
                        if (line.trim().isEmpty()) continue;
                        String[] fields = line.split("\t", -1);
                        if (fields.length < 2) continue;
                        String snpId = fields[0].trim();
                        Map<String, String> row = new LinkedHashMap<>();
                        for (int i = 1; i < Math.min(fields.length, headers.length); i++) {
                            row.put(headers[i].trim(), fields[i].trim());
                        }
                        snpResults.put(snpId, row);
                    }
                }

                // Use only the latest successful run
                break;
            } catch (IOException e) {
                System.err.printf("[AnalysisColumnProvider] Error reading run %s: %s%n",
                    rd.getName(), e.getMessage());
            }
        }
    }

    private String findLocusId(SnpContext ctx) {
        // Try to find locus ID from manifest by matching locus index
        File analysisDir = new File(projectDir, "loci_analysis");
        if (!analysisDir.isDirectory()) return null;

        File[] locusDirs = analysisDir.listFiles(File::isDirectory);
        if (locusDirs == null) return null;

        // Match by checking if the locus contains this SNP's position
        for (File ld : locusDirs) {
            if (ld.getName().equals("runs")) continue;
            File baseManifest = new File(ld, "base/base.manifest.json");
            if (baseManifest.exists()) {
                StepManifest sm = StepManifest.read(new File(ld, "base"), "base");
                if (sm != null && sm.inputs.containsKey("locus_id")) {
                    // Check if this locus matches the context's locus index
                    String region = sm.inputs.get("region");
                    if (region != null && region.startsWith(ctx.locusChr + ":")) {
                        return ld.getName();
                    }
                }
            }
        }
        return null;
    }
}
