import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/** Thread-safe pipeline progress shared between Main and LocalServer's /progress endpoint. */
public class ProgressTracker {
    public volatile String phase      = "Initializing";
    public volatile int    locusIndex = 0;
    public volatile int    totalLoci  = 0;
    public volatile boolean done      = false;
    public final Set<Integer> completedLoci = ConcurrentHashMap.newKeySet();

    // Overall progress: each phase gets a weight, progress within phase is per-locus
    private volatile int    phaseIndex  = 0;
    private volatile int    totalPhases = 1;

    public void setPhases(int total) {
        this.totalPhases = Math.max(1, total);
        this.phaseIndex  = 0;
    }

    public void update(String phase, int locusIndex, int total) {
        this.phase      = phase;
        this.locusIndex = locusIndex;
        this.totalLoci  = total;
    }

    public void nextPhase(String phase, int total) {
        this.phase      = phase;
        this.locusIndex = 0;
        this.totalLoci  = total;
        this.phaseIndex++;
    }

    public void advance(int locusIndex) {
        this.locusIndex = locusIndex;
    }

    /**
     * Overall pipeline percentage (0-99).
     * Each phase gets an equal share. Within a phase, progress is per-locus.
     */
    public int pct() {
        if (done) return 100;
        double phaseWeight = 100.0 / totalPhases;
        double completedPhases = Math.max(0, phaseIndex - 1) * phaseWeight;
        double withinPhase = totalLoci > 0
            ? (locusIndex * phaseWeight / totalLoci)
            : 0;
        return Math.min(99, (int)(completedPhases + withinPhase));
    }
}
