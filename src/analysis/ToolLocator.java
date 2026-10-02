import java.io.File;
import java.util.*;

/**
 * Finds the external programs the analysis tools call, so they run without the user editing PATH.
 *
 * Rscript, first match wins:
 *   1. LYNXGWAS_RSCRIPT environment variable (explicit override);
 *   2. Rscript on PATH;
 *   3. R_HOME/bin;
 *   4. standard Windows installs: %ProgramFiles%\R\R-x.y.z\bin and the per-user
 *      %LOCALAPPDATA%\Programs\R\R-x.y.z\bin, newest version first (compared numerically,
 *      so R-4.10.0 beats R-4.9.3).
 * Falls back to plain "Rscript" so the error message stays the familiar one if R is truly absent.
 *
 * Binaries (gcta64, GWAMA, magma): appRoot/bin/<name>(.exe) first, then PATH; null if not found.
 */
public class ToolLocator {

    private static volatile String cachedRscript;

    private ToolLocator() {}

    public static String rscript() {
        String c = cachedRscript;
        if (c == null) { c = findRscript(System.getenv(), windows()); cachedRscript = c; }
        return c;
    }

    static boolean windows() { return System.getProperty("os.name", "").toLowerCase(Locale.ROOT).contains("win"); }

    static String findRscript(Map<String, String> env, boolean windows) {
        String exe = windows ? "Rscript.exe" : "Rscript";
        String override = env.get("LYNXGWAS_RSCRIPT");
        if (override != null && !override.isEmpty() && new File(override).isFile()) return override;
        String onPath = onPath(env.get("PATH"), exe);
        if (onPath != null) return onPath;
        String rHome = env.get("R_HOME");
        if (rHome != null && !rHome.isEmpty()) {
            File f = new File(new File(rHome, "bin"), exe);
            if (f.isFile()) return f.getAbsolutePath();
        }
        if (windows) {
            List<File> roots = new ArrayList<>();
            for (String k : new String[]{"ProgramFiles", "ProgramW6432", "ProgramFiles(x86)"}) {
                String v = env.get(k); if (v != null && !v.isEmpty()) roots.add(new File(v, "R"));
            }
            String local = env.get("LOCALAPPDATA");
            if (local != null && !local.isEmpty()) roots.add(new File(new File(local, "Programs"), "R"));
            String best = newestRscript(roots, exe);
            if (best != null) return best;
        }
        return "Rscript";
    }

    /** Newest R-x.y.z/bin/Rscript under any of the given R install roots, or null. */
    static String newestRscript(List<File> roots, String exe) {
        File best = null; int[] bestV = null;
        for (File root : roots) {
            File[] dirs = root.listFiles(f -> f.isDirectory() && f.getName().startsWith("R-"));
            if (dirs == null) continue;
            for (File d : dirs) {
                File f = new File(new File(d, "bin"), exe);
                if (!f.isFile()) continue;
                int[] v = version(d.getName().substring(2));
                if (bestV == null || compare(v, bestV) > 0) { best = f; bestV = v; }
            }
        }
        return best == null ? null : best.getAbsolutePath();
    }

    static int[] version(String s) {
        String[] parts = s.split("[^0-9]+");
        int[] v = new int[3];
        for (int i = 0, k = 0; i < parts.length && k < 3; i++) {
            if (parts[i].isEmpty()) continue;
            try { v[k++] = Integer.parseInt(parts[i]); } catch (NumberFormatException ignored) {}
        }
        return v;
    }

    static int compare(int[] a, int[] b) {
        for (int i = 0; i < 3; i++) if (a[i] != b[i]) return Integer.compare(a[i], b[i]);
        return 0;
    }

    /** appRoot/bin/<name>(.exe), then PATH; absolute path or null. */
    public static String binary(File appRoot, String name) {
        String exe = windows() && !name.toLowerCase(Locale.ROOT).endsWith(".exe") ? name + ".exe" : name;
        File inBin = new File(new File(appRoot, "bin"), exe);
        if (inBin.isFile()) return inBin.getAbsolutePath();
        return onPath(System.getenv("PATH"), exe);
    }

    static String onPath(String path, String exe) {
        if (path == null) return null;
        for (String dir : path.split(File.pathSeparator)) {
            if (dir.isEmpty()) continue;
            File f = new File(dir, exe);
            if (f.isFile()) return f.getAbsolutePath();
        }
        return null;
    }
}
