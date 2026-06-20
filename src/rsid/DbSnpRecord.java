package rsid;

import java.util.List;

public class DbSnpRecord {
    public final String rsid;
    public final String ref;
    public final List<String> altList;
    public final Integer dbSnpBuild;

    public DbSnpRecord(String rsid, String ref, List<String> altList, Integer dbSnpBuild) {
        this.rsid       = rsid;
        this.ref        = ref.toUpperCase();
        this.altList    = altList;
        this.dbSnpBuild = dbSnpBuild;
    }

    @Override
    public String toString() {
        return rsid + " ref=" + ref + " alt=" + String.join(",", altList)
            + (dbSnpBuild != null ? " build=" + dbSnpBuild : "");
    }
}
