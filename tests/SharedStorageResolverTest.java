import rsid.GlobalConfig;
import rsid.SharedStorageResolver;
import java.io.*;
import java.nio.file.*;

/**
 * Standalone regression test (no external test framework — this project has none) for
 * SharedStorageResolver. Only exercises paths that need no network access (folder mode, and the
 * config-validation branches of git mode) — an actual git clone/pull is out of scope for an
 * offline test.
 */
public class SharedStorageResolverTest {

    public static void main(String[] args) throws Exception {
        int failures = 0;
        failures += testDirectPathWins();
        failures += testFolderModeFallback();
        failures += testNoneMode();
        failures += testGitModeMissingConfig();

        if (failures == 0) {
            System.out.println("PASS: all SharedStorageResolver tests passed");
        } else {
            System.out.println("FAIL: " + failures + " test(s) failed");
            System.exit(1);
        }
    }

    private static int testDirectPathWins() throws IOException {
        Path tmp = Files.createTempFile("sst-direct", ".txt");
        Files.writeString(tmp, "x");
        GlobalConfig.SharedStorage cfg = new GlobalConfig.SharedStorage();
        cfg.mode = "folder";
        cfg.folderPath = "/some/other/place/that/does/not/matter";

        SharedStorageResolver.Result r = SharedStorageResolver.resolve(cfg, tmp.toString());
        Files.deleteIfExists(tmp);

        if (!r.ok) {
            System.out.println("FAIL: direct existing path should resolve without consulting shared storage: " + r.message);
            return 1;
        }
        System.out.println("PASS: a path that already exists resolves directly, without needing shared storage");
        return 0;
    }

    private static int testFolderModeFallback() throws IOException {
        Path sharedDir = Files.createTempDirectory("sst-shared");
        Path sharedFile = sharedDir.resolve("gencode.gff3");
        Files.writeString(sharedFile, "##gff-version 3\n");

        GlobalConfig.SharedStorage cfg = new GlobalConfig.SharedStorage();
        cfg.mode = "folder";
        cfg.folderPath = sharedDir.toString();

        // User's own path doesn't exist locally, but the same filename exists in the shared folder.
        SharedStorageResolver.Result r = SharedStorageResolver.resolve(cfg,
            "/definitely/does/not/exist/gencode.gff3");

        Files.deleteIfExists(sharedFile);
        Files.deleteIfExists(sharedDir);

        if (!r.ok || !r.resolvedPath.contains("gencode.gff3")) {
            System.out.println("FAIL: expected fallback to shared folder to find gencode.gff3, got ok=" +
                r.ok + " message=" + r.message);
            return 1;
        }
        System.out.println("PASS: missing local path falls back to the shared folder by filename");
        return 0;
    }

    private static int testNoneMode() {
        GlobalConfig.SharedStorage cfg = new GlobalConfig.SharedStorage(); // mode = "none" by default
        SharedStorageResolver.Result r = SharedStorageResolver.resolve(cfg, "/does/not/exist.gff3");
        if (r.ok) {
            System.out.println("FAIL: mode=none must never resolve a missing path");
            return 1;
        }
        System.out.println("PASS: shared storage mode \"none\" never falls back, only the direct path is tried");
        return 0;
    }

    private static int testGitModeMissingConfig() {
        GlobalConfig.SharedStorage cfg = new GlobalConfig.SharedStorage();
        cfg.mode = "git"; // gitUrl/gitLocalClone left empty
        SharedStorageResolver.Result r = SharedStorageResolver.resolve(cfg, "/does/not/exist.gff3");
        if (r.ok) {
            System.out.println("FAIL: git mode with no local clone configured must fail cleanly, not silently succeed");
            return 1;
        }
        System.out.println("PASS: git mode with missing config fails with a clear message, no exception");
        return 0;
    }
}
