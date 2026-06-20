package export;

import java.util.*;

public class SnpContext {
    public String id;
    public String chr;
    public long pos;
    public double pvalue;
    public double negLog10P;
    public double r2WithLead = -1;
    public String ea, nea;
    public double beta = Double.NaN, oddsRatio = Double.NaN, se = Double.NaN;
    public double sampleN = Double.NaN, maf = Double.NaN, info = Double.NaN;
    public boolean isLead;
    public long distanceToLead;
    public String nearestGene = "";
    public long nearestGeneDist = Long.MAX_VALUE;
    // Locus context
    public int locusIndex;
    public String locusName;
    public String locusChr;
    public long locusStart, locusEnd;
    // Dynamic annotations (key → value)
    public Map<String, String> annotations = new LinkedHashMap<>();
}
