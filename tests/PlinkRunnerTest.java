import java.util.*;

/** PlinkRunner: failure reasons surface PLINK's own error line; the memory budget is released on failure. */
public class PlinkRunnerTest {
    static int failures = 0;

    static void check(boolean ok, String msg) {
        System.out.println((ok ? "PASS: " : "FAIL: ") + msg);
        if (!ok) failures++;
    }

    public static void main(String[] args) {
        PlinkRunner.Result oom = new PlinkRunner.Result();
        oom.exitCode = 7;
        oom.outputTail = "32454 MB RAM detected; reserving 1024 MB for main workspace.\nError: Out of memory.  The --memory flag may be helpful.";
        check(oom.reason().equals("exit 7: Error: Out of memory.  The --memory flag may be helpful."),
              "reason() picks PLINK's Error line: " + oom.reason());
        check(!oom.ok(), "non-zero exit is not ok");

        PlinkRunner.Result silent = new PlinkRunner.Result();
        silent.exitCode = 2;
        check(silent.reason().equals("exit 2"), "reason() with no output is just the exit code: " + silent.reason());

        int budget = PlinkRunner.budgetMb();
        check(budget >= PlinkRunner.FULL_PANEL_MB, "budget always fits one full-panel call (" + budget + " MB)");

        // A launch failure must give its memory back; otherwise later calls would block forever.
        // Asking for more than the whole budget each time means any leak blocks the next call.
        for (int i = 0; i < 3; i++) {
            try {
                PlinkRunner.run("no-such-plink-binary-xyz", Arrays.asList("--version"), Integer.MAX_VALUE);
                check(false, "missing binary should throw (attempt " + (i + 1) + ")");
            } catch (java.io.IOException e) {
                check(true, "missing binary throws IOException and releases the budget (attempt " + (i + 1) + ")");
            }
        }

        if (failures == 0) System.out.println("PASS: all PlinkRunner tests passed");
        System.exit(failures == 0 ? 0 : 1);
    }
}
