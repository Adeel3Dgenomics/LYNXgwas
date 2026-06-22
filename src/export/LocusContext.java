package export;

import java.util.*;

public class LocusContext {
    public int index;
    public String name;
    public String chr;
    public long start, end;
    public long sizeBp;
    public String leadSnpId = "";
    public long leadSnpPos;
    public double leadP = Double.NaN;
    public double leadBeta = Double.NaN;
    public double leadOr = Double.NaN;
    public String nearestGene = "";
    public long nearestGeneDist = Long.MAX_VALUE;
    public int nSnps;
    public int nP5e3, nP5e5, nP5e8;
    public String refPanel = "";
    public boolean ldComputed;
    // Inter-marker distance: gap from previous locus end to this locus start
    // "start" if first locus on this chromosome
    public String imd = "start";
    // Dynamic annotations
    public Map<String, String> annotations = new LinkedHashMap<>();
}
