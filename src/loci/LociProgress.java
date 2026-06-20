package loci;

public class LociProgress {
    public volatile String currentStep = "Initializing";
    public volatile int stepIndex      = 0;
    public volatile int totalSteps     = 6;
    public volatile boolean done       = false;
    public volatile boolean cancelled  = false;
    public volatile String error       = null;

    public volatile int candidateSnps     = 0;
    public volatile int seedSnps          = 0;
    public volatile int lociFound         = 0;
    public volatile int totalChromosomes  = 0;
    public volatile int currentChromosome = 0;

    public int pct() {
        if (done) return 100;
        double phaseWeight = 100.0 / totalSteps;
        double base = Math.max(0, stepIndex - 1) * phaseWeight;
        double within = totalChromosomes > 0 ? currentChromosome * phaseWeight / totalChromosomes : 0;
        return Math.min(99, (int)(base + within));
    }

    public String toJson() {
        return String.format(
            "{\"pct\":%d,\"current_step\":\"%s\",\"step_index\":%d,\"total_steps\":%d," +
            "\"done\":%s,\"error\":%s," +
            "\"candidate_snps\":%d,\"seed_snps\":%d,\"loci_found\":%d," +
            "\"total_chromosomes\":%d,\"current_chromosome\":%d}",
            pct(), esc(currentStep), stepIndex, totalSteps,
            done, error != null ? "\"" + esc(error) + "\"" : "null",
            candidateSnps, seedSnps, lociFound,
            totalChromosomes, currentChromosome);
    }

    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
