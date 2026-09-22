import rsid.GithubProjectSync;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

/**
 * Regression test for GithubProjectSync, against a real local `git init --bare` repository standing
 * in for "GitHub" — fast, no real network dependency, and a genuine exercise of the actual git
 * mechanics (init, remote, add, commit, push, clone, pull) rather than a mock. Uses the
 * *ToAnyRemote/*FromAnyRemote entry points so a plain local path is accepted (the github.com
 * restriction on push()/pullInto() themselves is tested separately, directly, below).
 */
public class GithubProjectSyncTest {

    public static void main(String[] progArgs) throws Exception {
        int failures = 0;

        // ── github.com URL validation (the actual production entry points) ──
        failures += check("accepts a real https://github.com/... URL",
            GithubProjectSync.isValidGithubHttpsUrl("https://github.com/user/repo.git"));
        failures += check("rejects a non-github host",
            !GithubProjectSync.isValidGithubHttpsUrl("https://gitlab.com/user/repo.git"));
        failures += check("rejects a non-https scheme",
            !GithubProjectSync.isValidGithubHttpsUrl("git@github.com:user/repo.git"));
        failures += check("rejects null", !GithubProjectSync.isValidGithubHttpsUrl(null));
        failures += check("push() rejects a non-github URL outright",
            !GithubProjectSync.push(Files.createTempDirectory("x").toFile(), "https://gitlab.com/a/b.git", "main", "msg", null).ok);
        failures += check("pullInto() rejects a non-github URL outright",
            !GithubProjectSync.pullInto(Files.createTempDirectory("x").toFile(), "https://gitlab.com/a/b.git", "main", null).ok);

        // ── Real local bare repo standing in for GitHub ──
        Path bareRepo = Files.createTempDirectory("lynxgwas_test_bare_");
        Files.delete(bareRepo);
        runGit(null, "init", "--bare", "-b", "main", bareRepo.toString());
        String repoUrl = bareRepo.toString().replace('\\', '/'); // git accepts a plain local path as a "URL"

        Path projectDir = Files.createTempDirectory("lynxgwas_test_project_");
        Files.writeString(projectDir.resolve("config.properties"), "gwas.file=input/test.tsv\n");

        GithubProjectSync.Result pushResult = GithubProjectSync.pushToAnyRemote(
            projectDir.toFile(), repoUrl, "main", "Initial push", null);
        failures += check("first push succeeds (init + remote + add + commit + push)", pushResult.ok);

        // ── Extract ("pull") into a brand-new, empty directory == a clone ──
        Path extractDir = Files.createTempDirectory("lynxgwas_test_extract_");
        Files.delete(extractDir);
        GithubProjectSync.Result cloneResult = GithubProjectSync.pullFromAnyRemote(
            extractDir.toFile(), repoUrl, "main", null);
        failures += check("extract into a fresh directory succeeds (real clone)", cloneResult.ok);
        failures += check("the extracted project really has the real file content",
            Files.exists(extractDir.resolve("config.properties")) &&
            Files.readString(extractDir.resolve("config.properties")).contains("input/test.tsv"));

        // ── A second push after a real local change actually updates the remote ──
        Files.writeString(projectDir.resolve("config.properties"), "gwas.file=input/updated.tsv\n");
        GithubProjectSync.Result secondPush = GithubProjectSync.pushToAnyRemote(
            projectDir.toFile(), repoUrl, "main", "Update", null);
        failures += check("second push (real change) succeeds", secondPush.ok);

        GithubProjectSync.Result pullAgain = GithubProjectSync.pullFromAnyRemote(
            extractDir.toFile(), repoUrl, "main", null);
        failures += check("pulling again into the same extracted clone succeeds", pullAgain.ok);
        failures += check("the pulled update is really reflected on disk",
            Files.readString(extractDir.resolve("config.properties")).contains("input/updated.tsv"));

        // ── Pushing with no changes at all is not treated as an error ──
        GithubProjectSync.Result noopPush = GithubProjectSync.pushToAnyRemote(
            projectDir.toFile(), repoUrl, "main", "no-op", null);
        failures += check("pushing with nothing changed still reports ok", noopPush.ok);

        // ── Refuses to clone into a non-empty, unrelated directory ──
        Path nonEmpty = Files.createTempDirectory("lynxgwas_test_nonempty_");
        Files.writeString(nonEmpty.resolve("unrelated.txt"), "pre-existing content");
        GithubProjectSync.Result refused = GithubProjectSync.pullFromAnyRemote(
            nonEmpty.toFile(), repoUrl, "main", null);
        failures += check("refuses to clone into a non-empty unrelated directory", !refused.ok);
        failures += check("the pre-existing file was left alone",
            Files.readString(nonEmpty.resolve("unrelated.txt")).equals("pre-existing content"));

        // ── Refuses to silently repoint an existing project's remote to a different repo ──
        Path secondBareRepo = Files.createTempDirectory("lynxgwas_test_bare2_");
        Files.delete(secondBareRepo);
        runGit(null, "init", "--bare", "-b", "main", secondBareRepo.toString());
        String secondRepoUrl = secondBareRepo.toString().replace('\\', '/');
        GithubProjectSync.Result repointAttempt = GithubProjectSync.pushToAnyRemote(
            projectDir.toFile(), secondRepoUrl, "main", "should be refused", null);
        failures += check("refuses to silently repoint an existing project's origin to a different repo",
            !repointAttempt.ok);

        if (failures == 0) {
            System.out.println("PASS: all GithubProjectSync tests passed");
        } else {
            System.out.println("FAIL: " + failures + " GithubProjectSync test(s) failed");
            System.exit(1);
        }
    }

    private static void runGit(File cwd, String... args) throws IOException, InterruptedException {
        List<String> cmd = new ArrayList<>();
        cmd.add("git");
        cmd.addAll(Arrays.asList(args));
        ProcessBuilder pb = new ProcessBuilder(cmd);
        if (cwd != null) pb.directory(cwd);
        pb.redirectErrorStream(true);
        Process p = pb.start();
        String out = new String(p.getInputStream().readAllBytes(), StandardCharsets.UTF_8);
        p.waitFor();
        if (p.exitValue() != 0) throw new IOException("git " + String.join(" ", args) + " failed: " + out);
    }

    private static int check(String label, boolean cond) {
        if (!cond) { System.out.println("FAIL: " + label); return 1; }
        System.out.println("PASS: " + label);
        return 0;
    }
}
