import analysis.SnakemakeSubmitter;
import rsid.GlobalConfig;
import java.util.*;

/**
 * Standalone regression test (no external test framework — this project has none) for
 * SnakemakeSubmitter. Only exercises pure logic (Snakefile generation, poll-interval selection,
 * remote-command argument-array construction, shell-metacharacter validation) — this class has
 * deliberately never been run against a real SSH host, per the hard safety constraint in
 * DECISIONS.md, so submit()/pollStatus() against a live target remain untested.
 */
public class SnakemakeSubmitterTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testSnakefileContent();
        failures += testPollIntervalScalesWithSize();
        failures += testPollIntervalRespectsConfiguredBounds();
        failures += testScpCommandIsArgumentArray();
        failures += testShellMetacharValidationRejectsInjection();
        failures += testShellMetacharValidationAllowsNormalValues();
        failures += testWrapperScriptContent();

        if (failures == 0) {
            System.out.println("PASS: all SnakemakeSubmitter tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    private static int testSnakefileContent() {
        String snake = SnakemakeSubmitter.generateSnakefile(
            Arrays.asList(1, 2, 3), "/refpanels/eur", "output/plink_subsets", "output/ld_results",
            Collections.singletonMap("plink", "plink/1.9"));
        int failures = 0;
        failures += assertContains(snake, "LOCI = [1, 2, 3]", "locus list");
        failures += assertContains(snake, "module load plink/1.9", "module load line");
        failures += assertContains(snake, "rule subset_locus", "subset rule");
        failures += assertContains(snake, "rule ld_locus", "LD rule");
        failures += assertContains(snake, "/refpanels/eur", "ref panel path");
        if (failures == 0) System.out.println("PASS: generated Snakefile contains expected rules/values");
        return failures;
    }

    private static int testPollIntervalScalesWithSize() {
        GlobalConfig.HpcConfig hpc = new GlobalConfig.HpcConfig();
        hpc.pollIntervalMinMinutes = 5;
        hpc.pollIntervalMaxMinutes = 20;
        int small = SnakemakeSubmitter.pickPollIntervalMinutes(hpc, 5);
        int large = SnakemakeSubmitter.pickPollIntervalMinutes(hpc, 500);
        if (!(small < large)) {
            System.out.println("FAIL: expected a larger job to pick a longer poll interval (small=" +
                small + ", large=" + large + ")");
            return 1;
        }
        System.out.println("PASS: poll interval scales up with job size (small=" + small + "min, large=" + large + "min)");
        return 0;
    }

    private static int testPollIntervalRespectsConfiguredBounds() {
        GlobalConfig.HpcConfig hpc = new GlobalConfig.HpcConfig();
        hpc.pollIntervalMinMinutes = 5;
        hpc.pollIntervalMaxMinutes = 20;
        int tiny = SnakemakeSubmitter.pickPollIntervalMinutes(hpc, 0);
        int huge = SnakemakeSubmitter.pickPollIntervalMinutes(hpc, 100_000);
        int failures = 0;
        if (tiny < 5 || tiny > 20) { System.out.println("FAIL: interval for 0 loci out of bounds: " + tiny); failures++; }
        if (huge < 5 || huge > 20) { System.out.println("FAIL: interval for huge job out of bounds: " + huge); failures++; }
        if (failures == 0) System.out.println("PASS: poll interval always stays within the configured [min,max] bounds");
        return failures;
    }

    private static int testScpCommandIsArgumentArray() throws Exception {
        GlobalConfig.HpcConfig hpc = new GlobalConfig.HpcConfig();
        hpc.sshHost = "cluster.example.edu"; hpc.sshUser = "jdoe"; hpc.sshKeyPath = "/home/jdoe/.ssh/id_ed25519";
        List<String> cmd = SnakemakeSubmitter.scpCommand(hpc, "/local/Snakefile", "/remote/work/Snakefile");
        int failures = 0;
        failures += check(cmd.get(0), "scp", "argv[0]");
        failures += check(cmd.contains("-i") ? "yes" : "no", "yes", "includes -i flag for key path");
        failures += check(cmd.get(cmd.size() - 1), "jdoe@cluster.example.edu:/remote/work/Snakefile", "destination arg");
        // No single element should itself contain a space-joined shell command - confirms this is
        // a real argument array, not one big string that a shell would re-split/interpret.
        for (String arg : cmd) {
            if (arg.contains(";") || arg.contains("&&") || arg.contains("|")) {
                System.out.println("FAIL: scp argument contains shell metacharacters: " + arg);
                failures++;
            }
        }
        if (failures == 0) System.out.println("PASS: scpCommand() builds a clean argument array, no shell string");
        return failures;
    }

    private static int testShellMetacharValidationRejectsInjection() {
        int failures = 0;
        String[] malicious = {"/scratch/x; rm -rf /", "/scratch/x && curl evil.com", "host`whoami`", "$(reboot)"};
        for (String val : malicious) {
            try {
                SnakemakeSubmitter.validateNoShellMetachars("test", val);
                System.out.println("FAIL: expected rejection of \"" + val + "\" but it was accepted");
                failures++;
            } catch (SnakemakeSubmitter.ValidationException expected) {
                // good
            }
        }
        if (failures == 0) System.out.println("PASS: shell-metacharacter validation rejects injection attempts");
        return failures;
    }

    private static int testShellMetacharValidationAllowsNormalValues() {
        String[] normal = {"/scratch/jdoe/lynxgwas", "login.cluster.example.edu", "jdoe", "plink/1.9"};
        for (String val : normal) {
            try {
                SnakemakeSubmitter.validateNoShellMetachars("test", val);
            } catch (SnakemakeSubmitter.ValidationException e) {
                System.out.println("FAIL: normal value wrongly rejected: \"" + val + "\" — " + e.getMessage());
                return 1;
            }
        }
        System.out.println("PASS: shell-metacharacter validation allows ordinary paths/hostnames/usernames");
        return 0;
    }

    private static int testWrapperScriptContent() throws Exception {
        GlobalConfig.HpcConfig hpc = new GlobalConfig.HpcConfig();
        hpc.remoteWorkDir = "/scratch/jdoe/lynxgwas";
        String script = SnakemakeSubmitter.generateRemoteWrapperScript(hpc, 8, 4);
        int failures = 0;
        failures += assertContains(script, "cd /scratch/jdoe/lynxgwas", "cd to work dir");
        failures += assertContains(script, "snakemake --cores 8 --jobs 4", "snakemake invocation");
        failures += assertContains(script, "set -e", "fails fast on error");
        if (failures == 0) System.out.println("PASS: remote wrapper script contains the expected fixed commands");
        return failures;
    }

    private static int assertContains(String haystack, String needle, String label) {
        if (!haystack.contains(needle)) {
            System.out.println("FAIL: expected " + label + " to contain \"" + needle + "\"");
            return 1;
        }
        return 0;
    }

    private static int check(String actual, String expected, String label) {
        if (!Objects.equals(actual, expected)) {
            System.out.println("FAIL: " + label + " — expected \"" + expected + "\", got \"" + actual + "\"");
            return 1;
        }
        return 0;
    }
}
