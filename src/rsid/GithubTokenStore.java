package rsid;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * Saves, forgets, and checks for a GitHub credential — without LYNXgwas itself ever storing, logging,
 * or returning the token as plaintext (DECISIONS_PHASE6.md section 2). All three operations shell out
 * to the user's own {@code git} binary's {@code credential} subcommand, feeding it exactly the same
 * stdin protocol a real git client would use internally
 * (https://git-scm.com/docs/git-credential#IOFMT); the actual secure storage is delegated entirely to
 * whatever credential helper the user's own git is already configured with (on this project's own
 * development machine, {@code credential.helper=manager} — Git Credential Manager, which stores
 * credentials in Windows Credential Manager, DPAPI-encrypted and tied to the Windows login). LYNXgwas
 * never implements its own secret storage.
 */
public class GithubTokenStore {

    private static final String HOST = "github.com";
    private static final long TIMEOUT_SECONDS = 15;

    public static class Result {
        public final boolean ok;
        public final String message;
        Result(boolean ok, String message) { this.ok = ok; this.message = message; }
    }

    /**
     * Hands the token to {@code git credential approve} over stdin so the user's own configured
     * credential helper stores it. The token is never written to a LYNXgwas file, never appears in a
     * command-line argument, and this method never logs it.
     */
    public static Result save(String token) {
        return save(token, Collections.emptyList());
    }

    /** Removes whatever credential is stored for github.com via {@code git credential reject}. */
    public static Result forget() {
        return forget(Collections.emptyList());
    }

    /**
     * Reports only whether a credential currently resolves for github.com — the resolved password is
     * read (git's own {@code fill} protocol requires reading the response) and discarded immediately
     * in this same method; it is never returned, logged, or stored by this class.
     */
    public static boolean isConfigured() {
        return isConfigured(Collections.emptyList());
    }

    // ── Overloads accepting extra leading `-c` git-config arguments, so tests can redirect the
    // credential helper to a disposable local store (git's own real, standard "store" helper against
    // a throwaway file) instead of ever touching the real, already-configured OS credential manager.
    // The zero-arg methods above always pass an empty override list, i.e. "use whatever the user's own
    // git is really configured with" — exactly the production path LocalServer actually calls.

    public static Result save(String token, List<String> extraConfigArgs) {
        if (token == null || token.isEmpty()) return new Result(false, "No token given");
        String input = "protocol=https\nhost=" + HOST + "\nusername=x-access-token\npassword=" + token + "\n\n";
        return run("approve", input, "Saved — your configured git credential helper now has it.", extraConfigArgs);
    }

    public static Result forget(List<String> extraConfigArgs) {
        String input = "protocol=https\nhost=" + HOST + "\n\n";
        return run("reject", input, "Removed.", extraConfigArgs);
    }

    public static boolean isConfigured(List<String> extraConfigArgs) {
        try {
            List<String> cmd = new ArrayList<>();
            cmd.add("git");
            cmd.addAll(extraConfigArgs);
            cmd.add("credential");
            cmd.add("fill");
            ProcessBuilder pb = new ProcessBuilder(cmd);
            pb.redirectErrorStream(false);
            Process p = pb.start();
            try (OutputStream os = p.getOutputStream()) {
                os.write(("protocol=https\nhost=" + HOST + "\n\n").getBytes(StandardCharsets.UTF_8));
            }
            String stdout = new String(p.getInputStream().readAllBytes(), StandardCharsets.UTF_8);
            boolean finished = p.waitFor(TIMEOUT_SECONDS, TimeUnit.SECONDS);
            if (!finished) { p.destroyForcibly(); return false; }
            // Deliberately: check for a non-empty password= line, then let `stdout` (which held the
            // real secret) go out of scope without ever being copied anywhere else.
            boolean has = false;
            for (String line : stdout.split("\n")) {
                if (line.startsWith("password=") && line.length() > "password=".length()) { has = true; break; }
            }
            return has;
        } catch (IOException | InterruptedException e) {
            return false;
        }
    }

    private static Result run(String subcommand, String stdinInput, String successMessage, List<String> extraConfigArgs) {
        try {
            List<String> cmd = new ArrayList<>();
            cmd.add("git");
            cmd.addAll(extraConfigArgs);
            cmd.add("credential");
            cmd.add(subcommand);
            ProcessBuilder pb = new ProcessBuilder(cmd);
            pb.redirectErrorStream(true);
            Process p = pb.start();
            try (OutputStream os = p.getOutputStream()) {
                os.write(stdinInput.getBytes(StandardCharsets.UTF_8));
            }
            // Read and discard stdout without ever logging it — "approve"/"reject" don't print the
            // secret back, but this defensively never assumes that of every possible credential helper.
            p.getInputStream().readAllBytes();
            boolean finished = p.waitFor(TIMEOUT_SECONDS, TimeUnit.SECONDS);
            if (!finished) { p.destroyForcibly(); return new Result(false, "git credential " + subcommand + " timed out"); }
            if (p.exitValue() != 0) return new Result(false, "git credential " + subcommand + " exited with code " + p.exitValue());
            return new Result(true, successMessage);
        } catch (IOException | InterruptedException e) {
            return new Result(false, "Failed to run git credential " + subcommand + ": " + e.getMessage());
        }
    }
}
