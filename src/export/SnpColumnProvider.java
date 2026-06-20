package export;

import java.util.List;

/**
 * Provides one or more columns of per-SNP data for Excel export.
 * To add a new annotation column: implement this interface and register it
 * in ExportRegistry — no change to ExcelExporter needed.
 */
public interface SnpColumnProvider {
    List<ColumnSpec> columns();
    Object value(SnpContext ctx, ColumnSpec c);
    boolean isAvailable(String projectDir);
}
