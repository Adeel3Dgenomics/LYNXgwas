package export;

import java.io.*;
import java.util.*;
import java.util.zip.*;

/**
 * Minimal xlsx writer using Java's built-in ZIP + XML.
 * Supports: multiple sheets, text/number/boolean cells, header freeze, auto-filter.
 */
public class XlsxWriter implements Closeable {

    private final ZipOutputStream zos;
    private final List<SheetData> sheets = new ArrayList<>();
    private final List<String> sharedStrings = new ArrayList<>();
    private final Map<String, Integer> ssIndex = new LinkedHashMap<>();

    public XlsxWriter(OutputStream out) {
        this.zos = new ZipOutputStream(out);
    }

    public Sheet addSheet(String name) {
        Sheet s = new Sheet(name, sheets.size());
        sheets.add(s.data);
        return s;
    }

    public void finish() throws IOException {
        // Write each sheet
        for (int i = 0; i < sheets.size(); i++) {
            writeSheet(i, sheets.get(i));
        }
        writeSharedStrings();
        writeStyles();
        writeWorkbook();
        writeContentTypes();
        writeRels();
        writeWorkbookRels();
        zos.finish();
    }

    @Override
    public void close() throws IOException { zos.close(); }

    private int ssi(String s) {
        return ssIndex.computeIfAbsent(s, k -> { sharedStrings.add(k); return sharedStrings.size() - 1; });
    }

    // ── Sheet building ───────────────────────────────────────────

    static class SheetData {
        String name;
        int index;
        List<Object[]> rows = new ArrayList<>();
        int frozenRows = 0;
        boolean autoFilter = false;
        int colCount = 0;
    }

    public class Sheet {
        final SheetData data;
        Sheet(String name, int index) {
            data = new SheetData();
            data.name = name;
            data.index = index;
        }
        public Sheet freezeHeader() { data.frozenRows = 1; return this; }
        public Sheet autoFilter() { data.autoFilter = true; return this; }
        public void addRow(Object... cells) {
            data.rows.add(cells);
            data.colCount = Math.max(data.colCount, cells.length);
        }
    }

    // ── XML writing ──────────────────────────────────────────────

    private void writeSheet(int idx, SheetData sd) throws IOException {
        zos.putNextEntry(new ZipEntry("xl/worksheets/sheet" + (idx + 1) + ".xml"));
        // Write directly to stream to avoid OOM on large sheets
        Writer w = new OutputStreamWriter(zos, "UTF-8");

        w.write("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n");
        w.write("<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">\n");

        // Dimension
        String lastCell = sd.rows.isEmpty() ? "A1" : colRef(Math.max(0, sd.colCount - 1)) + sd.rows.size();
        w.write("<dimension ref=\"A1:"); w.write(lastCell); w.write("\"/>\n");

        // OOXML element order: dimension → sheetViews → cols → sheetData → autoFilter

        // 1. Freeze panes (sheetViews must come first after dimension)
        if (sd.frozenRows > 0) {
            w.write("<sheetViews><sheetView tabSelected=\"");
            w.write(idx == 0 ? "1" : "0");
            w.write("\" workbookViewId=\"0\"><pane ySplit=\"");
            w.write(String.valueOf(sd.frozenRows));
            w.write("\" topLeftCell=\"A");
            w.write(String.valueOf(sd.frozenRows + 1));
            w.write("\" activePane=\"bottomLeft\" state=\"frozen\"/></sheetView></sheetViews>\n");
        }

        // sheetFormatPr (required between sheetViews and cols)
        w.write("<sheetFormatPr defaultRowHeight=\"15\"/>\n");

        // 2. Column widths
        if (sd.colCount > 0) {
            w.write("<cols>");
            for (int c = 0; c < sd.colCount; c++) {
                w.write("<col min=\""); w.write(String.valueOf(c + 1));
                w.write("\" max=\""); w.write(String.valueOf(c + 1));
                w.write("\" width=\"15\" bestFit=\"1\" customWidth=\"1\"/>");
            }
            w.write("</cols>\n");
        }

        // 3. Sheet data (rows)
        w.write("<sheetData>\n");
        for (int r = 0; r < sd.rows.size(); r++) {
            Object[] cells = sd.rows.get(r);
            w.write("<row r=\""); w.write(String.valueOf(r + 1)); w.write("\">");
            for (int c = 0; c < cells.length; c++) {
                String ref = colRef(c) + (r + 1);
                Object val = cells[c];
                if (val == null || "".equals(val)) {
                    continue;
                } else if (val instanceof Number) {
                    double d = ((Number) val).doubleValue();
                    if (val instanceof Double && (Double.isNaN(d) || Double.isInfinite(d))) {
                        w.write("<c r=\""); w.write(ref); w.write("\" t=\"s\"><v>");
                        w.write(String.valueOf(ssi("-"))); w.write("</v></c>");
                    } else {
                        String style = (r == 0) ? " s=\"1\"" :
                            (val instanceof Double && Math.abs(d) < 0.01 && d != 0) ? " s=\"2\"" : "";
                        w.write("<c r=\""); w.write(ref); w.write("\""); w.write(style);
                        w.write("><v>"); w.write(formatNum(d)); w.write("</v></c>");
                    }
                } else if (val instanceof Boolean) {
                    w.write("<c r=\""); w.write(ref); w.write("\" t=\"b\"><v>");
                    w.write((Boolean) val ? "1" : "0"); w.write("</v></c>");
                } else {
                    String s = val.toString();
                    w.write("<c r=\""); w.write(ref); w.write("\" t=\"s\"");
                    if (r == 0) w.write(" s=\"1\"");
                    w.write("><v>"); w.write(String.valueOf(ssi(s))); w.write("</v></c>");
                }
            }
            w.write("</row>\n");
            if (r % 10000 == 0) w.flush();
        }
        w.write("</sheetData>\n");

        // 4. Auto-filter (must come after sheetData)
        if (sd.autoFilter && !sd.rows.isEmpty()) {
            String range = "A1:" + colRef(sd.colCount - 1) + sd.rows.size();
            w.write("<autoFilter ref=\""); w.write(range); w.write("\"/>\n");
        }

        w.write("</worksheet>");
        w.flush();
        zos.closeEntry();
    }

    private String formatNum(double d) {
        if (d == (long) d && Math.abs(d) < 1e15) return String.valueOf((long) d);
        return String.valueOf(d);
    }

    private void writeSharedStrings() throws IOException {
        zos.putNextEntry(new ZipEntry("xl/sharedStrings.xml"));
        StringBuilder xml = new StringBuilder();
        xml.append("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n");
        xml.append("<sst xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" count=\"")
           .append(sharedStrings.size()).append("\" uniqueCount=\"").append(sharedStrings.size()).append("\">\n");
        for (String s : sharedStrings) {
            xml.append("<si><t>").append(escXml(s)).append("</t></si>\n");
        }
        xml.append("</sst>");
        zos.write(xml.toString().getBytes("UTF-8"));
        zos.closeEntry();
    }

    private void writeStyles() throws IOException {
        zos.putNextEntry(new ZipEntry("xl/styles.xml"));
        String xml = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n" +
            "<styleSheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">\n" +
            "<numFmts count=\"1\"><numFmt numFmtId=\"164\" formatCode=\"0.00E+00\"/></numFmts>\n" +
            "<fonts count=\"2\"><font><sz val=\"11\"/><name val=\"Calibri\"/></font>" +
            "<font><b/><sz val=\"11\"/><name val=\"Calibri\"/></font></fonts>\n" +
            "<fills count=\"2\"><fill><patternFill patternType=\"none\"/></fill>" +
            "<fill><patternFill patternType=\"gray125\"/></fill></fills>\n" +
            "<borders count=\"1\"><border><left/><right/><top/><bottom/><diagonal/></border></borders>\n" +
            "<cellStyleXfs count=\"1\"><xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\"/></cellStyleXfs>\n" +
            "<cellXfs count=\"3\">" +
            "<xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\"/>" +
            "<xf numFmtId=\"0\" fontId=\"1\" fillId=\"0\" borderId=\"0\" xfId=\"0\" applyFont=\"1\"/>" +
            "<xf numFmtId=\"164\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\" applyNumberFormat=\"1\"/>" +
            "</cellXfs>\n</styleSheet>";
        zos.write(xml.getBytes("UTF-8"));
        zos.closeEntry();
    }

    private void writeWorkbook() throws IOException {
        zos.putNextEntry(new ZipEntry("xl/workbook.xml"));
        StringBuilder xml = new StringBuilder();
        xml.append("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n");
        xml.append("<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" " +
            "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\">\n<sheets>\n");
        for (int i = 0; i < sheets.size(); i++) {
            xml.append("<sheet name=\"").append(escXml(sheets.get(i).name))
               .append("\" sheetId=\"").append(i + 1)
               .append("\" r:id=\"rId").append(i + 1).append("\"/>\n");
        }
        xml.append("</sheets>\n</workbook>");
        zos.write(xml.toString().getBytes("UTF-8"));
        zos.closeEntry();
    }

    private void writeContentTypes() throws IOException {
        zos.putNextEntry(new ZipEntry("[Content_Types].xml"));
        StringBuilder xml = new StringBuilder();
        xml.append("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n");
        xml.append("<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">\n");
        xml.append("<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>\n");
        xml.append("<Default Extension=\"xml\" ContentType=\"application/xml\"/>\n");
        xml.append("<Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>\n");
        xml.append("<Override PartName=\"/xl/sharedStrings.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml\"/>\n");
        xml.append("<Override PartName=\"/xl/styles.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml\"/>\n");
        for (int i = 0; i < sheets.size(); i++) {
            xml.append("<Override PartName=\"/xl/worksheets/sheet").append(i + 1)
               .append(".xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>\n");
        }
        xml.append("</Types>");
        zos.write(xml.toString().getBytes("UTF-8"));
        zos.closeEntry();
    }

    private void writeRels() throws IOException {
        zos.putNextEntry(new ZipEntry("_rels/.rels"));
        String xml = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n" +
            "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\n" +
            "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/>\n" +
            "</Relationships>";
        zos.write(xml.getBytes("UTF-8"));
        zos.closeEntry();
    }

    private void writeWorkbookRels() throws IOException {
        zos.putNextEntry(new ZipEntry("xl/_rels/workbook.xml.rels"));
        StringBuilder xml = new StringBuilder();
        xml.append("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n");
        xml.append("<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\n");
        for (int i = 0; i < sheets.size(); i++) {
            xml.append("<Relationship Id=\"rId").append(i + 1)
               .append("\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet")
               .append(i + 1).append(".xml\"/>\n");
        }
        xml.append("<Relationship Id=\"rId").append(sheets.size() + 1)
           .append("\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/sharedStrings\" Target=\"sharedStrings.xml\"/>\n");
        xml.append("<Relationship Id=\"rId").append(sheets.size() + 2)
           .append("\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles\" Target=\"styles.xml\"/>\n");
        xml.append("</Relationships>");
        zos.write(xml.toString().getBytes("UTF-8"));
        zos.closeEntry();
    }

    private static String colRef(int col) {
        StringBuilder sb = new StringBuilder();
        col++;
        while (col > 0) { col--; sb.insert(0, (char) ('A' + col % 26)); col /= 26; }
        return sb.toString();
    }

    private static String escXml(String s) {
        if (s == null) return "";
        return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
                .replace("\"", "&quot;").replace("'", "&apos;");
    }
}
