import java.io.File;
import java.nio.file.Files;
import java.util.*;

/** ToolLocator: Rscript discovery order and numeric version choice, on throwaway fake installs. */
public class ToolLocatorTest {
    static int failures = 0;

    public static void main(String[] args) throws Exception {
        File tmp = Files.createTempDirectory("toollocator").toFile();
        File pf = new File(tmp, "ProgramFiles"), local = new File(tmp, "Local");
        fakeR(new File(pf, "R"), "R-4.9.3");
        fakeR(new File(new File(local, "Programs"), "R"), "R-4.10.0");
        fakeR(new File(new File(local, "Programs"), "R"), "R-4.2.1");

        Map<String, String> env = new HashMap<>();
        env.put("PATH", ""); env.put("ProgramFiles", pf.getPath()); env.put("LOCALAPPDATA", local.getPath());
        String got = ToolLocator.findRscript(env, true);
        check("newest version across both install roots (4.10.0 beats 4.9.3): " + got, got.contains("R-4.10.0"));

        File override = new File(tmp, "custom/Rscript.exe"); override.getParentFile().mkdirs(); override.createNewFile();
        env.put("LYNXGWAS_RSCRIPT", override.getPath());
        check("LYNXGWAS_RSCRIPT override wins", ToolLocator.findRscript(env, true).equals(override.getPath()));
        env.remove("LYNXGWAS_RSCRIPT");

        File pathDir = new File(tmp, "onpath"); pathDir.mkdirs(); new File(pathDir, "Rscript.exe").createNewFile();
        env.put("PATH", pathDir.getPath());
        check("Rscript on PATH beats install folders", ToolLocator.findRscript(env, true).startsWith(pathDir.getAbsolutePath()));

        Map<String, String> none = new HashMap<>(); none.put("PATH", "");
        check("falls back to plain Rscript when nothing is found", ToolLocator.findRscript(none, true).equals("Rscript"));
        check("version parse", Arrays.equals(ToolLocator.version("4.10.2"), new int[]{4, 10, 2}));

        if (failures == 0) System.out.println("PASS: all ToolLocator tests passed");
        else { System.out.println("FAIL: " + failures + " ToolLocator check(s) failed"); System.exit(1); }
    }

    static void fakeR(File root, String name) throws Exception {
        File bin = new File(new File(root, name), "bin"); bin.mkdirs(); new File(bin, "Rscript.exe").createNewFile();
    }

    static void check(String name, boolean ok) {
        System.out.println((ok ? "PASS: " : "FAIL: ") + name);
        if (!ok) failures++;
    }
}
