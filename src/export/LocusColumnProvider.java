package export;

import java.util.List;

/**
 * Provides one or more columns of per-locus data for Excel export.
 * To add a new locus annotation: implement this interface and register it.
 */
public interface LocusColumnProvider {
    List<ColumnSpec> columns();
    Object value(LocusContext ctx, ColumnSpec c);
    boolean isAvailable(String projectDir);
}
