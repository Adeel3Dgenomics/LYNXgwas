package rsid;

import java.util.concurrent.atomic.*;

public class RsidProgress {
    public volatile String currentStep = "Initializing";
    public volatile int stepIndex   = 0;
    public volatile int totalSteps  = 5;
    public volatile boolean done    = false;
    public volatile boolean cancelled = false;
    public volatile String error    = null;

    // Live counters
    public final AtomicInteger matched  = new AtomicInteger(0);
    public final AtomicInteger forward  = new AtomicInteger(0);
    public final AtomicInteger reverse  = new AtomicInteger(0);
    public final AtomicInteger unmatched = new AtomicInteger(0);
    public volatile int totalSnps  = 0;
    public volatile int currentLocus = 0;
    public volatile int totalLoci   = 0;

    // Result
    public volatile double recoveryRate = 0;
    public volatile String outputFile   = "";

    public int pct() {
        if (done) return 100;
        double phaseWeight = 100.0 / totalSteps;
        double base = Math.max(0, stepIndex - 1) * phaseWeight;
        double within = 0;
        if (totalLoci > 0) within = currentLocus * phaseWeight / totalLoci;
        return Math.min(99, (int)(base + within));
    }

    public String toJson() {
        return String.format(
            "{\"pct\":%d,\"current_step\":\"%s\",\"step_index\":%d,\"total_steps\":%d," +
            "\"done\":%s,\"error\":%s," +
            "\"matched\":%d,\"forward\":%d,\"reverse\":%d,\"unmatched\":%d," +
            "\"total_snps\":%d,\"current_locus\":%d,\"total_loci\":%d," +
            "\"recovery_rate\":%.1f,\"output_file\":\"%s\"}",
            pct(), esc(currentStep), stepIndex, totalSteps,
            done, error != null ? "\"" + esc(error) + "\"" : "null",
            matched.get(), forward.get(), reverse.get(), unmatched.get(),
            totalSnps, currentLocus, totalLoci,
            recoveryRate, esc(outputFile));
    }

    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
