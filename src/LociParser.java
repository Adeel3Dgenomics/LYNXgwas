import java.io.*;
import java.util.*;

public class LociParser {

    public static List<Locus> parse(Config config) throws IOException {
        List<Locus> loci = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(config.lociFile))) {
            String header = br.readLine();
            if (header == null) throw new IOException("loci.txt is empty");

            // Map column names to indices
            String[] cols = header.trim().split("\t");
            int idxChr   = indexOf(cols, "meta_chr");
            int idxStart = indexOf(cols, "meta_start");
            int idxEnd   = indexOf(cols, "meta_end");
            if (idxChr < 0 || idxStart < 0 || idxEnd < 0) {
                throw new IOException("loci.txt must have columns: meta_chr, meta_start, meta_end");
            }

            // Load optional top-SNP overrides
            Map<Integer, String> topSnpOverrides = loadTopSnpFile(config.topSnpFile);

            String line;
            int idx = 1;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] f = line.trim().split("\t");
                String chr = f[idxChr].trim();
                long start = Long.parseLong(f[idxStart].trim());
                long end   = Long.parseLong(f[idxEnd].trim());
                Locus locus = new Locus(idx, chr, start, end, config.locusPadding);
                loci.add(locus);
                idx++;
            }
        }
        System.out.printf("[LociParser] Loaded %d loci%n", loci.size());
        return loci;
    }

    private static int indexOf(String[] cols, String name) {
        for (int i = 0; i < cols.length; i++) {
            if (cols[i].trim().equalsIgnoreCase(name)) return i;
        }
        return -1;
    }

    private static Map<Integer, String> loadTopSnpFile(String path) throws IOException {
        Map<Integer, String> map = new HashMap<>();
        if (path == null || path.isEmpty()) return map;
        File f = new File(path);
        if (!f.exists()) { System.err.println("[WARN] top.snp.file not found: " + path); return map; }
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty() || line.startsWith("#")) continue;
                String[] parts = line.trim().split("\t");
                if (parts.length >= 2) {
                    map.put(Integer.parseInt(parts[0].trim()), parts[1].trim());
                }
            }
        }
        return map;
    }
}
