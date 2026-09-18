package rsid;

import java.io.*;
import java.util.concurrent.TimeUnit;

/**
 * Resolves a reference-data path against optional shared storage (see GlobalConfig.SharedStorage)
 * so a user doesn't have to re-download the same PLINK reference panel / GFF3 annotation on every
 * machine or project. Never handles git credentials itself — "git" mode only ever shells out to
 * the user's own already-configured git binary via an argument-array ProcessBuilder (never a
 * concatenated shell string), so an SSH deploy key or HTTPS credential helper the user has set up
 * themselves keeps working exactly as it would from their own terminal.
 */
public class SharedStorageResolver {

    public static class Result {
        public boolean ok;
        public String message;
        public String resolvedPath; // set only when ok

        static Result ok(String path, String message) {
            Result r = new Result(); r.ok = true; r.resolvedPath = path; r.message = message; return r;
        }
        static Result fail(String message) {
            Result r = new Result(); r.ok = false; r.message = message; return r;
        }
    }

    /** Resolves userPath: if it already exists as given, use it as-is; otherwise, if shared
     *  storage is configured, look for the same filename under the shared location. Never
     *  performs a network operation (no implicit git clone) — call ensureGitClone() explicitly
     *  first if git mode should fetch/update the local clone. */
    public static Result resolve(GlobalConfig.SharedStorage cfg, String userPath) {
        if (userPath == null || userPath.isEmpty()) return Result.fail("No path given");

        File direct = new File(userPath);
        if (direct.exists()) return Result.ok(direct.getAbsolutePath(), "Found at given path");

        if (cfg == null || cfg.mode == null || cfg.mode.equals("none")) {
            return Result.fail("Not found: " + userPath);
        }

        String filename = new File(userPath).getName();
        if ("folder".equals(cfg.mode)) {
            if (cfg.folderPath.isEmpty()) return Result.fail("Shared folder mode is set but no folder path is configured");
            File candidate = new File(cfg.folderPath, filename);
            if (candidate.exists()) return Result.ok(candidate.getAbsolutePath(), "Found in shared folder");
            return Result.fail("Not found at " + userPath + " or in shared folder " + cfg.folderPath);
        }

        if ("git".equals(cfg.mode)) {
            if (cfg.gitLocalClone.isEmpty()) return Result.fail("Git mode is set but no local clone path is configured");
            File candidate = new File(cfg.gitLocalClone, filename);
            if (candidate.exists()) return Result.ok(candidate.getAbsolutePath(), "Found in git-backed shared storage");
            File cloneDir = new File(cfg.gitLocalClone);
            if (!cloneDir.isDirectory()) {
                return Result.fail("Not found at " + userPath + "; shared git clone does not exist yet at " +
                    cfg.gitLocalClone + " — use the Test/Sync button in Resources & Settings to clone it first");
            }
            return Result.fail("Not found at " + userPath + " or in the cloned shared repo at " + cfg.gitLocalClone);
        }

        return Result.fail("Not found: " + userPath);
    }

    /** Clones cfg.gitUrl to cfg.gitLocalClone if it doesn't exist yet, or runs "git pull" there if
     *  it does. Shells out to the user's own git binary with an argument-array ProcessBuilder —
     *  never a concatenated shell string, never touches credentials directly. */
    public static Result ensureGitClone(GlobalConfig.SharedStorage cfg) {
        if (cfg.gitUrl.isEmpty()) return Result.fail("No git URL configured");
        if (cfg.gitLocalClone.isEmpty()) return Result.fail("No local clone path configured");

        File cloneDir = new File(cfg.gitLocalClone);
        try {
            ProcessBuilder pb;
            if (cloneDir.isDirectory() && new File(cloneDir, ".git").exists()) {
                pb = new ProcessBuilder("git", "-C", cloneDir.getAbsolutePath(), "pull");
            } else {
                File parent = cloneDir.getParentFile();
                if (parent != null) parent.mkdirs();
                pb = new ProcessBuilder("git", "clone", cfg.gitUrl, cloneDir.getAbsolutePath());
            }
            pb.redirectErrorStream(true);
            Process p = pb.start();
            String output = new String(p.getInputStream().readAllBytes(), "UTF-8");
            boolean finished = p.waitFor(120, TimeUnit.SECONDS);
            if (!finished) { p.destroyForcibly(); return Result.fail("git operation timed out after 120s"); }
            if (p.exitValue() != 0) {
                return Result.fail("git exited with code " + p.exitValue() + ": " +
                    (output.length() > 500 ? output.substring(0, 500) + "..." : output));
            }
            return Result.ok(cloneDir.getAbsolutePath(), "git repo ready at " + cloneDir.getAbsolutePath());
        } catch (IOException | InterruptedException e) {
            return Result.fail("Failed to run git: " + e.getMessage());
        }
    }

    /** Validates the current shared-storage config without necessarily fetching anything (used by
     *  the "Test" button) — checks a folder path exists, or that git is reachable/clonable. */
    public static Result test(GlobalConfig.SharedStorage cfg) {
        if (cfg == null || cfg.mode == null || cfg.mode.equals("none")) {
            return Result.fail("Shared storage is not enabled (mode = none)");
        }
        if ("folder".equals(cfg.mode)) {
            if (cfg.folderPath.isEmpty()) return Result.fail("No folder path configured");
            File f = new File(cfg.folderPath);
            if (!f.isDirectory()) return Result.fail("Folder does not exist or is not a directory: " + cfg.folderPath);
            return Result.ok(f.getAbsolutePath(), "Folder is reachable: " + f.getAbsolutePath());
        }
        if ("git".equals(cfg.mode)) {
            return ensureGitClone(cfg);
        }
        return Result.fail("Unknown shared storage mode: " + cfg.mode);
    }
}
