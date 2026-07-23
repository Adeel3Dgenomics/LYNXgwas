import loci.LociProgress;

/** Progress tracker for a Locus Matrix job: merge -> identify -> scan -> done. */
public class MultiLocusProgress {
    public volatile String phase          = "merge"; // merge | identify | scan | done | error
    public volatile int datasetIndex      = 0;
    public volatile int datasetTotal      = 0;
    public volatile String currentDataset = "";
    public volatile boolean done          = false;
    public volatile String error          = null;
    public volatile int lociFound         = 0;

    public final LociProgress identifyProgress = new LociProgress();

    public int pct() {
        if (done) return 100;
        double band;
        double within;
        switch (phase) {
            case "merge":
                band = 0;
                within = datasetTotal > 0 ? (double) datasetIndex / datasetTotal : 0;
                break;
            case "identify":
                band = 1;
                within = identifyProgress.pct() / 100.0;
                break;
            case "scan":
                band = 2;
                within = datasetTotal > 0 ? (double) datasetIndex / datasetTotal : 0;
                break;
            default:
                band = 3;
                within = 0;
        }
        return Math.min(99, (int) ((band + within) * 100.0 / 3));
    }

    public String toJson() {
        return String.format(
            "{\"pct\":%d,\"phase\":\"%s\",\"dataset_index\":%d,\"dataset_total\":%d," +
            "\"current_dataset\":\"%s\",\"loci_found\":%d,\"done\":%s,\"error\":%s}",
            pct(), esc(phase), datasetIndex, datasetTotal,
            esc(currentDataset), lociFound, done,
            error != null ? "\"" + esc(error) + "\"" : "null");
    }

    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
