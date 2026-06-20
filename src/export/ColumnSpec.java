package export;

public class ColumnSpec {
    public final String key;
    public final String header;
    public final String group;
    public final ColType type;
    public final int order;

    public enum ColType { TEXT, INT, DOUBLE, SCIENTIFIC, BOOL }

    public ColumnSpec(String key, String header, String group, ColType type, int order) {
        this.key = key; this.header = header; this.group = group;
        this.type = type; this.order = order;
    }
}
