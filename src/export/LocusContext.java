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
    public String leadEa = "";
    public String leadNea = "";
    public double leadSe = Double.NaN;
    public double leadN = Double.NaN;
    public double leadMaf = Double.NaN;
    public double leadInfo = Double.NaN;
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
    // Lead SNP's GWAS columns not part of the standardized pipeline (e.g. METAL's
    // Direction/HetISq/HetChiSq/HetDf/HetPVal, MR-MEGA's per-cohort stats)
    public Map<String, String> leadExtra = new LinkedHashMap<>();
}
