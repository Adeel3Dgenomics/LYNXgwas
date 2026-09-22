import rsid.GithubTokenStore;
import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Regression test for GithubTokenStore. Deliberately never touches the real, already-configured OS
 * credential manager (Git Credential Manager / Windows Credential Manager) — doing so during an
 * automated test run would both pollute the developer's real credential store with a fake test token
 * and risk clobbering the real GitHub credential already relied on elsewhere. Instead this redirects
 * git's credential helper, for these calls only, to git's own real, standard "store" helper backed by
 * a throwaway temp file (via the package-private {@code extraConfigArgs} overloads) — a real git
 * credential helper, just a disposable one, so this genuinely exercises the same
 * approve/reject/fill protocol GithubTokenStore uses in production without touching anything real.
 */
public class GithubTokenStoreTest {

    public static void main(String[] progArgs) throws Exception {
        int failures = 0;

        Path storeFile = Files.createTempFile("lynxgwas_test_credstore_", ".txt");
        Files.deleteIfExists(storeFile); // must not exist yet for a clean test
        // Git config values treat "\" as an escape character, so a raw Windows path (with backslashes)
        // embedded in a "-c credential.helper=store --file=..." value gets silently mangled — found by
        // this test actually failing the first time, not assumed. Forward slashes work fine on Windows
        // git and are not special in a config value.
        String storeFilePath = storeFile.toAbsolutePath().toString().replace('\\', '/');
        List<String> testHelper = Arrays.asList(
            "-c", "credential.helper=", // clear any inherited helper first
            "-c", "credential.helper=store --file=" + storeFilePath
        );

        try {
            failures += check("isConfigured() is false before anything is saved",
                !GithubTokenStore.isConfigured(testHelper));

            GithubTokenStore.Result saveResult = GithubTokenStore.save("ghp_fake_test_token_never_real", testHelper);
            failures += check("save() reports ok", saveResult.ok);
            failures += check("isConfigured() is true after saving", GithubTokenStore.isConfigured(testHelper));

            // The disposable store file is a real credential helper's real backing file — read it
            // directly only to confirm approve() actually persisted something plausible via the real
            // protocol (this is verifying the test's own fixture is real, not LYNXgwas reading back a
            // credential in production code, which GithubTokenStore itself is designed to never do).
            String stored = Files.readString(storeFile);
            failures += check("store file contains the github.com host", stored.contains("github.com"));
            failures += check("store file contains the expected username", stored.contains("x-access-token"));

            GithubTokenStore.Result forgetResult = GithubTokenStore.forget(testHelper);
            failures += check("forget() reports ok", forgetResult.ok);
            failures += check("isConfigured() is false again after forgetting",
                !GithubTokenStore.isConfigured(testHelper));

            GithubTokenStore.Result emptyTokenResult = GithubTokenStore.save("", testHelper);
            failures += check("save('') is rejected rather than silently storing an empty secret",
                !emptyTokenResult.ok);
        } finally {
            Files.deleteIfExists(storeFile);
        }

        if (failures == 0) {
            System.out.println("PASS: all GithubTokenStore tests passed");
        } else {
            System.out.println("FAIL: " + failures + " GithubTokenStore test(s) failed");
            System.exit(1);
        }
    }

    private static int check(String label, boolean cond) {
        if (!cond) { System.out.println("FAIL: " + label); return 1; }
        System.out.println("PASS: " + label);
        return 0;
    }
}
