import rsid.GlobalConfig;

/**
 * Standalone regression test (no external test framework — this project has none) for the
 * SharedStorage/HpcConfig sections added to GlobalConfig's hand-rolled JSON round-trip.
 */
public class GlobalConfigTest {

    public static void main(String[] args) {
        int failures = 0;
        failures += testRoundTrip();
        failures += testDefaultsAreSafe();

        if (failures == 0) {
            System.out.println("PASS: all GlobalConfig tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    private static int testRoundTrip() {
        GlobalConfig gc = new GlobalConfig();
        gc.sharedStorage.mode = "git";
        gc.sharedStorage.folderPath = "/lab/shared/refdata";
        gc.sharedStorage.gitUrl = "git@github.com:example/lynxgwas-refdata.git";
        gc.sharedStorage.gitLocalClone = "/home/user/.lynxgwas/refdata";

        gc.hpc.enabled = true;
        gc.hpc.sshHost = "login.cluster.example.edu";
        gc.hpc.sshUser = "jdoe";
        gc.hpc.sshKeyPath = "/home/user/.ssh/id_ed25519";
        gc.hpc.remoteWorkDir = "/scratch/jdoe/lynxgwas";
        gc.hpc.pollIntervalMinMinutes = 7;
        gc.hpc.pollIntervalMaxMinutes = 15;
        gc.hpc.moduleNames.put("plink", "plink/1.9");
        gc.hpc.moduleNames.put("r", "R/4.3.1");

        String json = gc.toJson();
        GlobalConfig back = GlobalConfig.parseForTest(json);

        int failures = 0;
        failures += check("sharedStorage.mode", back.sharedStorage.mode, gc.sharedStorage.mode);
        failures += check("sharedStorage.folderPath", back.sharedStorage.folderPath, gc.sharedStorage.folderPath);
        failures += check("sharedStorage.gitUrl", back.sharedStorage.gitUrl, gc.sharedStorage.gitUrl);
        failures += check("sharedStorage.gitLocalClone", back.sharedStorage.gitLocalClone, gc.sharedStorage.gitLocalClone);
        failures += check("hpc.enabled", String.valueOf(back.hpc.enabled), String.valueOf(gc.hpc.enabled));
        failures += check("hpc.sshHost", back.hpc.sshHost, gc.hpc.sshHost);
        failures += check("hpc.sshUser", back.hpc.sshUser, gc.hpc.sshUser);
        failures += check("hpc.sshKeyPath", back.hpc.sshKeyPath, gc.hpc.sshKeyPath);
        failures += check("hpc.remoteWorkDir", back.hpc.remoteWorkDir, gc.hpc.remoteWorkDir);
        failures += check("hpc.pollIntervalMinMinutes", String.valueOf(back.hpc.pollIntervalMinMinutes), "7");
        failures += check("hpc.pollIntervalMaxMinutes", String.valueOf(back.hpc.pollIntervalMaxMinutes), "15");
        failures += check("hpc.moduleNames.size", String.valueOf(back.hpc.moduleNames.size()), "2");
        failures += check("hpc.moduleNames[plink]", back.hpc.moduleNames.get("plink"), "plink/1.9");
        failures += check("hpc.moduleNames[r]", back.hpc.moduleNames.get("r"), "R/4.3.1");

        if (failures == 0) {
            System.out.println("PASS: SharedStorage/HpcConfig round-trip through toJson()/parse() correctly");
        }
        return failures;
    }

    private static int testDefaultsAreSafe() {
        GlobalConfig gc = new GlobalConfig();
        int failures = 0;
        if (gc.hpc.enabled) {
            System.out.println("FAIL: HpcConfig.enabled must default to false (opt-in only)");
            failures++;
        }
        if (!gc.sharedStorage.mode.equals("none")) {
            System.out.println("FAIL: SharedStorage.mode must default to \"none\"");
            failures++;
        }
        if (failures == 0) {
            System.out.println("PASS: new config sections default to fully inert (no behavior change for existing users)");
        }
        return failures;
    }

    private static int check(String label, String actual, String expected) {
        if (!java.util.Objects.equals(actual, expected)) {
            System.out.println("FAIL: " + label + " — expected \"" + expected + "\", got \"" + actual + "\"");
            return 1;
        }
        return 0;
    }
}
