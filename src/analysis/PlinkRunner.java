import java.io.*;
import java.util.*;
import java.util.concurrent.Semaphore;

/**
 * Runs PLINK for the per-locus base steps with a bounded memory footprint.
 *
 * Without --memory, PLINK 1.9 reserves half of physical RAM per process. When the base build
 * runs several loci in parallel next to the server's own heap, those reservations exhaust RAM,
 * PLINK exits with "Out of memory" before writing any output, and the step used to report a
 * misleading "no variants" error. Each call here gets an explicit --memory / --threads 1 and
 * draws its megabytes from a shared budget, so parallel builds queue instead of failing.
 *
 * Scanning the full reference panel (tens of millions of .bim lines) needs ~3 GB; per-locus
 * subsets need far less.
 */
public class PlinkRunner {

    /** Memory for a PLINK call that reads the whole reference panel .bim. */
    public static final int FULL_PANEL_MB = 4096;
    /** Memory for a PLINK call on a per-locus subset. */
    public static final int SUBSET_MB = 2048;

    private static final int OUTPUT_TAIL_LINES = 12;
    private static final Semaphore BUDGET = new Semaphore(budgetMb(), true);

    public static class Result {
        public int exitCode = -1;
        /** Last lines of PLINK's combined stdout/stderr, for error messages. */
        public String outputTail = "";

        public boolean ok() { return exitCode == 0; }

        /** Short reason for a failure, e.g. "exit 7: Error: Out of memory." */
        public String reason() {
            String tail = outputTail.trim();
            String last = "";
            for (String line : tail.split("\n"))
                if (line.startsWith("Error") || line.contains("rror:")) { last = line.trim(); break; }
            if (last.isEmpty() && !tail.isEmpty()) last = tail.substring(tail.lastIndexOf('\n') + 1).trim();
            return "exit " + exitCode + (last.isEmpty() ? "" : ": " + last);
        }
    }

    /**
     * Run PLINK with {@code args} (everything after the binary). Adds --memory and --threads 1
     * unless the caller already set them, and blocks until {@code memoryMb} fits in the budget.
     */
    public static Result run(String plinkBin, List<String> args, int memoryMb) throws IOException {
        List<String> cmd = new ArrayList<>();
        cmd.add(plinkBin);
        cmd.addAll(args);
        if (!args.contains("--memory")) { cmd.add("--memory"); cmd.add(String.valueOf(memoryMb)); }
        if (!args.contains("--threads")) { cmd.add("--threads"); cmd.add("1"); }

        int permits = Math.min(memoryMb, budgetMb());
        Result result = new Result();
        try {
            BUDGET.acquire(permits);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            result.outputTail = "interrupted while waiting for memory budget";
            return result;
        }
        try {
            Process proc = new ProcessBuilder(cmd).redirectErrorStream(true).start();
            Deque<String> tail = new ArrayDeque<>();
            try (BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()))) {
                String line;
                while ((line = br.readLine()) != null) {
                    tail.addLast(line);
                    if (tail.size() > OUTPUT_TAIL_LINES) tail.removeFirst();
                }
            }
            result.exitCode = proc.waitFor();
            result.outputTail = String.join("\n", tail);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            result.outputTail = "interrupted";
        } finally {
            BUDGET.release(permits);
        }
        return result;
    }

    /**
     * MB available to concurrent PLINK processes: physical RAM minus this JVM's max heap and a
     * 2 GB reserve for the OS, never below one full-panel call. LYNXGWAS_PLINK_BUDGET_MB overrides.
     */
    static int budgetMb() {
        String env = System.getenv("LYNXGWAS_PLINK_BUDGET_MB");
        if (env != null) {
            try { return Math.max(SUBSET_MB, Integer.parseInt(env.trim())); } catch (NumberFormatException ignored) {}
        }
        long totalMb = 0;
        try {
            java.lang.management.OperatingSystemMXBean os =
                java.lang.management.ManagementFactory.getOperatingSystemMXBean();
            if (os instanceof com.sun.management.OperatingSystemMXBean)
                totalMb = ((com.sun.management.OperatingSystemMXBean) os).getTotalMemorySize() >> 20;
        } catch (Throwable ignored) {}
        if (totalMb <= 0) return FULL_PANEL_MB;
        long heapMb = Runtime.getRuntime().maxMemory() >> 20;
        long budget = totalMb - heapMb - 2048;
        return (int) Math.max(FULL_PANEL_MB, Math.min(budget, Integer.MAX_VALUE));
    }
}
