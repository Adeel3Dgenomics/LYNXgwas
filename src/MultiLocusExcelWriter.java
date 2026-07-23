import java.io.*;
import java.util.*;
import export.XlsxWriter;

/** Builds the Locus Matrix export workbook from an already-computed MultiLocusResult. */
public class MultiLocusExcelWriter {

    public static void write(MultiLocusResult result, OutputStream out) throws IOException {
        try (XlsxWriter xlsx = new XlsxWriter(out)) {
            XlsxWriter.Sheet s = xlsx.addSheet("Locus Matrix").freezeHeader().autoFilter();

            List<Object> header = new ArrayList<>();
            header.add("Locus");
            header.add("Chr");
            header.add("Start");
            header.add("End");
            header.add("Size_bp");
            header.add("Nearest/Overlapping Gene");
            header.add("Gene Dist");
            header.add("Overall Best SNP");
            header.add("Overall Best P");
            header.add("Overall Best Dataset");
            for (MultiLocusResult.DatasetInfo d : result.datasets) {
                header.add(d.name + " Best SNP");
                header.add(d.name + " EA");
                header.add(d.name + " NEA");
                header.add(d.name + " P");
                header.add(d.name + " " + d.effectType);
                header.add(d.name + " SE");
                header.add(d.name + " MAF");
                header.add(d.name + " N");
                header.add(d.name + " Info");
                header.add(d.name + " N p<5e-5");
                header.add(d.name + " N p<5e-8");
                for (String extraCol : d.extraColumns) {
                    header.add(d.name + " " + extraCol);
                }
            }
            s.addRow(header.toArray());

            for (MultiLocusResult.LocusRow l : result.loci) {
                List<Object> row = new ArrayList<>();
                row.add("L" + (l.index + 1));
                row.add(l.chr);
                row.add(l.start);
                row.add(l.end);
                row.add(l.sizeBp);
                row.add(dash(l.nearestGene));
                row.add(l.geneDist >= 0 ? l.geneDist : "-");
                row.add(dash(l.overallBestSnpId));
                row.add(num(l.overallBestP));
                row.add(dash(l.overallBestDatasetId));

                for (MultiLocusResult.DatasetInfo d : result.datasets) {
                    MultiLocusResult.DatasetLocusStat stat = l.cells.get(d.id);
                    if (stat == null) {
                        row.add("-"); row.add("-"); row.add("-"); row.add("-");
                        row.add("-"); row.add("-"); row.add("-"); row.add("-");
                        row.add("-"); row.add(0);   row.add(0);
                        for (int i = 0; i < d.extraColumns.size(); i++) row.add("-");
                        continue;
                    }
                    row.add(dash(stat.bestSnpId));
                    row.add(dash(stat.ea));
                    row.add(dash(stat.nea));
                    row.add(num(stat.bestP));
                    row.add(num("OR".equals(d.effectType) ? stat.or : stat.beta));
                    row.add(num(stat.se));
                    row.add(num(stat.maf));
                    row.add(num(stat.n));
                    row.add(num(stat.info));
                    row.add(stat.nP5e5);
                    row.add(stat.nP5e8);
                    for (String extraCol : d.extraColumns) {
                        String v = stat.extra.get(extraCol);
                        row.add((v == null || v.isEmpty()) ? "-" : v);
                    }
                }
                s.addRow(row.toArray());
            }

            xlsx.finish();
        }
    }

    private static Object dash(String s) {
        return (s == null || s.isEmpty()) ? "-" : s;
    }

    private static Object num(double d) {
        return Double.isNaN(d) ? "-" : d;
    }
}
