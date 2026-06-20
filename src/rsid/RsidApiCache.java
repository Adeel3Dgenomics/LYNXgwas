package rsid;

import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Disk-backed cache for API rsID lookups, keyed by build:chr:pos:ref:alt.
 */
public class RsidApiCache {
    private final Map<String, String> cache = new LinkedHashMap<>();
    private final Path cacheFile;

    public RsidApiCache(String outputDir) {
        this.cacheFile = Paths.get(outputDir, "rsid_api_cache.tsv");
        load();
    }

    public String get(String build, String chr, long pos, String ea, String nea) {
        return cache.get(key(build, chr, pos, ea, nea));
    }

    public void put(String build, String chr, long pos, String ea, String nea, String rsid) {
        cache.put(key(build, chr, pos, ea, nea), rsid);
    }

    public void save() {
        try {
            Files.createDirectories(cacheFile.getParent());
            try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(cacheFile.toFile())))) {
                for (Map.Entry<String, String> e : cache.entrySet())
                    pw.println(e.getKey() + "\t" + e.getValue());
            }
        } catch (IOException e) {
            System.err.println("[RsidApiCache] Save failed: " + e.getMessage());
        }
    }

    private void load() {
        if (!Files.exists(cacheFile)) return;
        try (BufferedReader br = new BufferedReader(new FileReader(cacheFile.toFile()))) {
            String line;
            while ((line = br.readLine()) != null) {
                int tab = line.indexOf('\t');
                if (tab > 0) cache.put(line.substring(0, tab), line.substring(tab + 1));
            }
        } catch (IOException e) {
            System.err.println("[RsidApiCache] Load failed: " + e.getMessage());
        }
    }

    private static String key(String build, String chr, long pos, String ea, String nea) {
        return build + ":" + chr + ":" + pos + ":" + ea.toUpperCase() + ":" + nea.toUpperCase();
    }
}
