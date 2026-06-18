/** Thread-safe pipeline progress shared between Main and LocalServer's /progress endpoint. */
public class ProgressTracker {
    public volatile String phase      = "Initializing";
    public volatile int    locusIndex = 0;
    public volatile int    totalLoci  = 0;
    public volatile boolean done      = false;

    public void update(String phase, int locusIndex, int total) {
        this.phase      = phase;
        this.locusIndex = locusIndex;
        this.totalLoci  = total;
    }

    public int pct() {
        return totalLoci > 0 ? Math.min(99, (int)(locusIndex * 100.0 / totalLoci)) : 0;
    }
}
