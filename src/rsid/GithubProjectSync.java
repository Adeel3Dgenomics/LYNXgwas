package rsid;

import java.io.*;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * Pushes a whole project directory to, and pulls ("extracts") one from, a GitHub repository
 * (DECISIONS_PHASE6.md) — the project-level counterpart to {@link SharedStorageResolver}'s existing
 * reference-data git support, following the exact same conventions: shell out to the user's own
 * {@code git} binary via an argument-array {@link ProcessBuilder} (never a concatenated shell string),
 * never touch a credential directly, and never hang waiting on an interactive terminal prompt
 * ({@code GIT_TERMINAL_PROMPT=0} is set on every invocation this class makes).
 *
 * A one-time token (not saved anywhere) is injected only through one child process's own environment
 * block via {@code GIT_ASKPASS} — see {@link #exec} and the static, secret-free askpass helper script
 * this writes once per call. When no one-time token is supplied, git falls back to whatever credential
 * helper is already configured (see {@link GithubTokenStore} for the "save it the safe way" path that
 * feeds that same helper).
 *
 * {@link #pushToAnyRemote}/{@link #pullFromAnyRemote} are public (not just {@link #push}/
 * {@link #pullInto}) so tests can exercise the real git mechanics against a fast local repository
 * instead of a real network dependency, while {@code LocalServer} itself only ever calls the
 * github.com-restricted {@link #push}/{@link #pullInto}.
 */
public class GithubProjectSync {

    private static final long TIMEOUT_SECONDS = 120;

    public static class Result {
        public final boolean ok;
        public final String message;
        Result(boolean ok, String message) { this.ok = ok; this.message = message; }
        static Result ok(String msg) { return new Result(true, msg); }
        static Result fail(String msg) { return new Result(false, msg); }
    }

    /** Rejects anything but a plain https://github.com/... URL — a real safety boundary (predictable
     *  credential host, no accidental operation against an arbitrary remote), not just tidiness. */
    public static boolean isValidGithubHttpsUrl(String url) {
        if (url == null || url.isEmpty()) return false;
        try {
            URI u = URI.create(url);
            return "https".equalsIgnoreCase(u.getScheme()) && "github.com".equalsIgnoreCase(u.getHost());
        } catch (IllegalArgumentException e) {
            return false;
        }
    }

    /** Push a project directory to a GitHub repo: init-if-needed, add, commit (skipped cleanly if
     *  nothing changed), push. */
    public static Result push(File projectDir, String repoUrl, String branch, String commitMessage,
                               String oneTimeToken) {
        if (!isValidGithubHttpsUrl(repoUrl)) return Result.fail("Only https://github.com/... URLs are accepted");
        return pushToAnyRemote(projectDir, repoUrl, branch, commitMessage, oneTimeToken);
    }

    /** Same operation without the github.com restriction — see the class javadoc for why this is public. */
    public static Result pushToAnyRemote(File projectDir, String repoUrl, String branch, String commitMessage,
                                          String oneTimeToken) {
        if (!projectDir.isDirectory()) return Result.fail("Project directory does not exist: " + projectDir);
        String br = (branch == null || branch.isEmpty()) ? "main" : branch;
        String msg = (commitMessage == null || commitMessage.isEmpty()) ? "Update from LYNXgwas" : commitMessage;

        File gitDir = new File(projectDir, ".git");
        try {
            if (!gitDir.isDirectory()) {
                Result r = exec(projectDir, oneTimeToken, "init", "-b", br);
                if (!r.ok) return r;
                r = exec(projectDir, oneTimeToken, "remote", "add", "origin", repoUrl);
                if (!r.ok) return r;
            } else {
                // Already a repo — make sure "origin" points at the requested URL rather than silently
                // pushing an unrelated project to whatever remote happened to be configured before.
                Result remoteCheck = exec(projectDir, oneTimeToken, "remote", "get-url", "origin");
                if (!remoteCheck.ok) {
                    Result r = exec(projectDir, oneTimeToken, "remote", "add", "origin", repoUrl);
                    if (!r.ok) return r;
                } else if (!remoteCheck.message.trim().equals(repoUrl)) {
                    return Result.fail("This project's local .git already has a different 'origin' (" +
                        remoteCheck.message.trim() + ") — refusing to silently repoint it to " + repoUrl);
                }
            }

            Result add = exec(projectDir, oneTimeToken, "add", "-A");
            if (!add.ok) return add;

            Result status = exec(projectDir, oneTimeToken, "status", "--porcelain");
            if (status.ok && status.message.trim().isEmpty()) {
                // Nothing changed since the last push — not an error, just nothing to do beyond that.
                Result pushOnly = exec(projectDir, oneTimeToken, "push", "-u", "origin", br);
                return pushOnly.ok
                    ? Result.ok("Nothing new to commit; pushed existing history to " + repoUrl)
                    : pushOnly;
            }

            Result commit = exec(projectDir, oneTimeToken, "commit", "-m", msg);
            if (!commit.ok) return commit;

            Result push = exec(projectDir, oneTimeToken, "push", "-u", "origin", br);
            if (!push.ok) return push;
            return Result.ok("Pushed to " + repoUrl + " (branch " + br + ")");
        } catch (Exception e) {
            return Result.fail("Unexpected error: " + e.getMessage());
        }
    }

    /** Pull ("extract") a project from a GitHub repo into targetDir: clone if targetDir is missing or
     *  empty, pull if it's already a clone of that same repo, refuse otherwise rather than overwriting
     *  unrelated local content. */
    public static Result pullInto(File targetDir, String repoUrl, String branch, String oneTimeToken) {
        if (!isValidGithubHttpsUrl(repoUrl)) return Result.fail("Only https://github.com/... URLs are accepted");
        return pullFromAnyRemote(targetDir, repoUrl, branch, oneTimeToken);
    }

    /** Same operation without the github.com restriction — see {@link #pushToAnyRemote} for why. */
    public static Result pullFromAnyRemote(File targetDir, String repoUrl, String branch, String oneTimeToken) {
        String br = (branch == null || branch.isEmpty()) ? "main" : branch;

        boolean exists = targetDir.isDirectory();
        boolean empty = !exists || isEmptyDir(targetDir);
        File gitDir = new File(targetDir, ".git");

        if (empty) {
            File parent = targetDir.getParentFile();
            if (parent != null) parent.mkdirs();
            if (exists) targetDir.delete(); // clone requires the target to not already exist (or be genuinely empty)
            Result r = exec(parent != null ? parent : new File("."), oneTimeToken,
                "clone", "--branch", br, repoUrl, targetDir.getAbsolutePath());
            return r.ok ? Result.ok("Cloned " + repoUrl + " into " + targetDir) : r;
        }

        if (gitDir.isDirectory()) {
            Result remoteCheck = exec(targetDir, oneTimeToken, "remote", "get-url", "origin");
            if (remoteCheck.ok && remoteCheck.message.trim().equals(repoUrl)) {
                Result r = exec(targetDir, oneTimeToken, "pull", "origin", br);
                return r.ok ? Result.ok("Pulled latest from " + repoUrl + " into " + targetDir) : r;
            }
            return Result.fail("Target directory already has a different git remote — refusing to overwrite it");
        }

        return Result.fail("Target directory already exists, is not empty, and is not a clone of this repo: " + targetDir);
    }

    private static boolean isEmptyDir(File dir) {
        String[] entries = dir.list();
        return entries == null || entries.length == 0;
    }

    // ── Process execution ──

    private static Result exec(File cwd, String oneTimeToken, String... gitArgs) {
        List<String> cmd = new ArrayList<>();
        cmd.add("git");
        cmd.addAll(Arrays.asList(gitArgs));
        File askpassScript = null;
        try {
            ProcessBuilder pb = new ProcessBuilder(cmd);
            pb.directory(cwd);
            pb.redirectErrorStream(true);
            Map<String, String> env = pb.environment();
            env.put("GIT_TERMINAL_PROMPT", "0"); // never hang waiting on a terminal that doesn't exist
            if (oneTimeToken != null && !oneTimeToken.isEmpty()) {
                askpassScript = writeAskpassHelper();
                env.put("GIT_ASKPASS", askpassScript.getAbsolutePath());
                env.put("LYNX_GIT_TOKEN", oneTimeToken); // lives only in this child process's own environment block
            }
            Process p = pb.start();
            String output = new String(p.getInputStream().readAllBytes(), StandardCharsets.UTF_8);
            boolean finished = p.waitFor(TIMEOUT_SECONDS, TimeUnit.SECONDS);
            if (!finished) { p.destroyForcibly(); return Result.fail("git " + gitArgs[0] + " timed out after " + TIMEOUT_SECONDS + "s"); }
            if (p.exitValue() != 0) {
                return Result.fail("git " + String.join(" ", gitArgs) + " exited with code " + p.exitValue() + ": " +
                    (output.length() > 800 ? output.substring(0, 800) + "..." : output));
            }
            return Result.ok(output);
        } catch (IOException | InterruptedException e) {
            return Result.fail("Failed to run git " + gitArgs[0] + ": " + e.getMessage());
        } finally {
            // The helper script itself never contained the secret (see writeAskpassHelper), but it's
            // still cleaned up promptly rather than left lying around indefinitely.
            if (askpassScript != null) askpassScript.delete();
        }
    }

    /**
     * Writes a small, fixed, secret-free helper script and returns its path. The script contains no
     * token — it only echoes back either the literal username GitHub expects for token auth
     * ("x-access-token") or the {@code LYNX_GIT_TOKEN} environment variable, depending on which prompt
     * git is asking. The real secret is never written to this file; it only ever exists in the calling
     * child process's own environment block (see {@link #exec}), which disappears when that process
     * exits.
     */
    private static File writeAskpassHelper() throws IOException {
        File f = File.createTempFile("lynxgwas_askpass_", ".cmd");
        String script =
            "@echo off\r\n" +
            "echo %~1 | findstr /I \"Username\" >nul\r\n" +
            "if %errorlevel%==0 (\r\n" +
            "  echo x-access-token\r\n" +
            ") else (\r\n" +
            "  echo %LYNX_GIT_TOKEN%\r\n" +
            ")\r\n";
        Files.write(f.toPath(), script.getBytes(StandardCharsets.UTF_8));
        return f;
    }
}
