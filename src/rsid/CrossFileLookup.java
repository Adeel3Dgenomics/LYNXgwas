package rsid;

import java.io.*;
import java.util.*;
import java.util.regex.*;

/**
 * Strategy 3 for rsID recovery: after local dbSNP matching and (optional) API
 * completion, cross-references remaining unmatched chr:pos SNPs against every
 * project's GWAS file that has a genuine rsID column configured — including
 * the current project's own file, in case it's one of the sources. If the
 * same variant appears in another cohort's summary stats with a real rsID
 * already resolved, that's used to fill the gap. Whatever's still unmatched
 * after this goes to the cross-project "missing rsIDs" list.
 */
public class CrossFileLookup {

    private static final Pattern RSID_PATTERN = Pattern.compile("^rs\\d+$", Pattern.CASE_INSENSITIVE);

    /**
     * Builds a chr:pos -> rsid index from every project under projectsRoot whose
     * config.properties has a non-empty col.rsid pointing at a real rsID column.
     */
    public static Map<String, String> buildIndex(File projectsRoot) {
        Map<String, String> index = new HashMap<>();
        File[] dirs = projectsRoot.listFiles(File::isDirectory);
        if (dirs == null) return index;

        for (File dir : dirs) {
            try {
                File cfgFile = new File(dir, "config.properties");
                if (!cfgFile.exists()) continue;
                Map<String, String> props = readProperties(cfgFile);
                String gwasFile = props.get("gwas.file");
                String colChr = props.getOrDefault("col.chr", "chrom");
                String colPos = props.getOrDefault("col.pos", "pos");
                String colRsid = props.get("col.rsid");

                if (colRsid == null || colRsid.isEmpty()) continue;
                if (gwasFile == null || !new File(gwasFile).exists()) continue;

                int added = indexOneFile(gwasFile, colChr, colPos, colRsid, index);
                System.out.printf("[CrossFileLookup] %s: indexed %,d rsIDs from %s%n",
                    dir.getName(), added, gwasFile);
            } catch (Exception e) {
                System.err.println("[CrossFileLookup] Skipping " + dir.getName() + ": " + e.getMessage());
            }
        }
        System.out.printf("[CrossFileLookup] Total cross-file index size: %,d%n", index.size());
        return index;
    }

    /** Minimal config.properties reader — only needs a handful of scalar keys here. */
    private static Map<String, String> readProperties(File file) throws IOException {
        Map<String, String> props = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#") || line.startsWith("!")) continue;
                int eq = line.indexOf('=');
                if (eq < 0) continue;
                String key = line.substring(0, eq).trim();
                String val = line.substring(eq + 1).trim().replace("\\", "/");
                props.put(key, val);
            }
        }
        return props;
    }

    private static int indexOneFile(String gwasFile, String colChr, String colPos, String colRsid,
                                    Map<String, String> index) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile), 1 << 20)) {
            String header = br.readLine();
            if (header == null) return 0;
            String[] cols = header.trim().split("\t");
            int iChr = colIdx(cols, colChr);
            int iPos = colIdx(cols, colPos);
            int iRsid = colIdx(cols, colRsid);
            if (iChr < 0 || iPos < 0 || iRsid < 0) return 0;
            int maxIdx = Math.max(iChr, Math.max(iPos, iRsid));

            int added = 0;
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t", -1);
                if (f.length <= maxIdx) continue;
                String rsid = f[iRsid].trim();
                if (!RSID_PATTERN.matcher(rsid).matches()) continue;
                String key = f[iChr].trim() + ":" + f[iPos].trim();
                if (index.putIfAbsent(key, rsid) == null) added++;
            }
            return added;
        }
    }

    private static int colIdx(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++)
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        return -1;
    }
}
