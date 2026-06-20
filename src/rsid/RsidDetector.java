package rsid;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class RsidDetector {

    private static final Pattern RSID_PATTERN = Pattern.compile("^rs\\d+$", Pattern.CASE_INSENSITIVE);
    private static final double THRESHOLD = 0.5;

    public static class DetectionResult {
        public final boolean present;
        public final String column;
        public final double confidence;

        public DetectionResult(boolean present, String column, double confidence) {
            this.present    = present;
            this.column     = column;
            this.confidence = confidence;
        }
    }

    /**
     * Pure function: detect whether sample rows contain a usable rsID column.
     * A column qualifies if >50% of its non-empty values match ^rs\d+$.
     */
    public static DetectionResult detectRsidColumn(String[] headerCols, List<String[]> sampleRows) {
        if (headerCols == null || headerCols.length == 0 || sampleRows.isEmpty())
            return new DetectionResult(false, null, 0);

        String bestCol = null;
        double bestConf = 0;

        for (int c = 0; c < headerCols.length; c++) {
            int total = 0, matches = 0;
            for (String[] row : sampleRows) {
                if (c >= row.length) continue;
                String val = row[c].trim();
                if (val.isEmpty() || val.equals(".") || val.equalsIgnoreCase("NA")) continue;
                total++;
                if (RSID_PATTERN.matcher(val).matches()) matches++;
            }
            if (total == 0) continue;
            double conf = (double) matches / total;
            if (conf > bestConf) {
                bestConf = conf;
                bestCol = headerCols[c].trim();
            }
        }

        if (bestConf >= THRESHOLD && bestCol != null)
            return new DetectionResult(true, bestCol, bestConf);
        return new DetectionResult(false, null, bestConf);
    }

    /**
     * Detect from a GWAS file directly. Reads header + first ~1000 rows.
     */
    public static DetectionResult detectFromFile(String gwasFile) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(gwasFile))) {
            String headerLine = br.readLine();
            if (headerLine == null) return new DetectionResult(false, null, 0);
            String[] header = headerLine.trim().split("\t");

            List<String[]> sample = new ArrayList<>();
            String line;
            int count = 0;
            while ((line = br.readLine()) != null && count < 1000) {
                if (line.isEmpty()) continue;
                sample.add(line.split("\t", -1));
                count++;
            }
            return detectRsidColumn(header, sample);
        }
    }
}
