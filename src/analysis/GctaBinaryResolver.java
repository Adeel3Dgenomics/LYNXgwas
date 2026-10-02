
import java.io.File;
import java.io.IOException;

/**
 * GCTA binary resolution for CojoAdapter (--cojo-slct conditional analysis on GWAS summary
 * statistics + a reference panel). (GCTA-GREML, which needs individual-level genotypes and
 * phenotypes, was removed: LYNXgwas runs only methods that work from summary statistics.)
 *
 * Two entry points, matching the two ways callers currently hand in a GCTA path:
 *
 *  - verify(String): the original CojoAdapter behavior, preserved byte-for-byte (same error
 *    message text) — treat the given path as literal (resolved against the working directory
 *    if relative) and fail loudly if it isn't there. This is what LocalServer's endpoints use
 *    today (they already resolve a "gcta_path" param, defaulting to "bin/gcta64.exe", themselves).
 *
 *  - find(File appRoot): a fuller resolution chain for a caller that doesn't already have an
 *    explicit configured path in hand: appRoot/bin/gcta64(.exe) first, then PATH — mirroring
 *    the config-value -> remembered-path -> PATH order described for external tool binaries
 *    elsewhere in this app (PlinkSubsetter.findPlink() only does the first of these; this adds
 *    the PATH fallback PlinkSubsetter does not currently have).
 */
public class GctaBinaryResolver {

    /**
     * Verifies an explicit GCTA path, exactly as CojoAdapter has always done. Kept as its own
     * method (rather than inlined in CojoAdapter) so any future GCTA-based summary-statistics tool
     * reuses the identical check/error text.
     */
    public static File verify(String gctaBin) throws IOException {
        File gctaFile = new File(gctaBin);
        if (!gctaFile.isAbsolute()) gctaFile = gctaFile.getAbsoluteFile();
        if (!gctaFile.exists())
            throw new IOException("GCTA binary not found at: " + gctaFile.getAbsolutePath()
                + ". Place gcta64.exe in the bin/ folder.");
        return gctaFile;
    }

    /**
     * Resolves a GCTA binary without an explicit configured path: appRoot/bin/gcta64(.exe),
     * then PATH. Throws a clear IOException (never returns null) if nothing is found.
     */
    public static String find(File appRoot) throws IOException {
        boolean windows = System.getProperty("os.name", "").toLowerCase().contains("win");
        String exeName = windows ? "gcta64.exe" : "gcta64";

        File inBin = new File(appRoot, "bin" + File.separator + exeName);
        if (inBin.exists()) return inBin.getAbsolutePath();

        String path = System.getenv("PATH");
        if (path != null) {
            for (String dir : path.split(File.pathSeparator)) {
                if (dir.isEmpty()) continue;
                File candidate = new File(dir, exeName);
                if (candidate.exists()) return candidate.getAbsolutePath();
            }
        }

        throw new IOException("GCTA binary not found at: " + inBin.getAbsolutePath()
            + " (or on PATH). Place gcta64.exe in the bin/ folder.");
    }
}
