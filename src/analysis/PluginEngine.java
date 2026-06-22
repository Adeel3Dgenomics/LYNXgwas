
import java.io.*;
import java.util.*;
import java.util.concurrent.*;

/**
 * Discovers tool descriptors, prepares inputs, runs tools via ProcessBuilder,
 * maps raw output to the unified contract, and validates results.
 *
 * Adding a tool = drop a YAML descriptor in the tools/ directory.
 * No engine code changes needed.
 */
public class PluginEngine {

    private static final String TOOLS_DIR = "tools";

    public static class RunRequest {
        public String tool;
        public Map<String, String> params = new LinkedHashMap<>();
        public String locusId;
        public String projectId;
        public String projectDir;
        public String preCreatedRunDir;
    }

    public static class RunResult {
        public boolean ok;
        public String error;
        public String jobId;
        public String tool;
        public String runDir;
        public int resultRows;
        public String diagnosticVerdict;
        public OutputContractValidator.ValidationResult validation;
    }

    /**
     * Discover all available tool descriptors.
     */
    public static List<ToolDescriptor> discoverTools() {
        List<ToolDescriptor> tools = new ArrayList<>();
        File dir = new File(TOOLS_DIR);
        if (!dir.isDirectory()) return tools;

        File[] files = dir.listFiles((d, name) ->
            name.endsWith(".yaml") || name.endsWith(".yml"));
        if (files == null) return tools;

        for (File f : files) {
            try {
                ToolDescriptor td = ToolDescriptor.parse(f);
                if (td.tool != null && !td.tool.isEmpty()) {
                    tools.add(td);
                    System.out.printf("[PluginEngine] Discovered tool: %s v%s (%s)%n",
                        td.tool, td.version, td.language);
                }
            } catch (Exception e) {
                System.err.printf("[PluginEngine] Failed to parse %s: %s%n",
                    f.getName(), e.getMessage());
            }
        }
        return tools;
    }

    /**
     * Find a specific tool descriptor by name.
     */
    public static ToolDescriptor findTool(String toolName) {
        for (ToolDescriptor td : discoverTools()) {
            if (td.tool.equals(toolName)) return td;
        }
        return null;
    }

    /**
     * Execute a tool run: prepare inputs → run command → map output → validate.
     */
    public static RunResult execute(RunRequest req, Config config, Locus locus,
                                     ProgressTracker progress) throws IOException {
        RunResult result = new RunResult();
        result.tool = req.tool;
        result.jobId = UUID.randomUUID().toString();

        ToolDescriptor td = findTool(req.tool);
        if (td == null) {
            result.error = "Tool not found: " + req.tool;
            return result;
        }

        // Verify base artifacts exist
        File analysisRoot = BaseStepPipeline.analysisDir(config, locus);
        File harmonizedDir = new File(analysisRoot, "harmonized");
        File ldDir = new File(analysisRoot, "ld");
        File matchedDir = new File(analysisRoot, "matched");

        for (String req_artifact : td.requires) {
            File check = null;
            switch (req_artifact) {
                case "harmonized_gwas": check = new File(harmonizedDir, "harmonized_gwas.tsv"); break;
                case "matched_ref": check = new File(matchedDir, "matched_ref.bed"); break;
                case "ld_r": check = new File(ldDir, "ld_r.matrix"); break;
                case "ld_r2": check = new File(ldDir, "ld_r2.matrix"); break;
                case "ld_dprime": check = new File(ldDir, "ld_dprime.matrix"); break;
            }
            if (check != null && !check.exists()) {
                result.error = "Required artifact missing: " + req_artifact +
                    ". Run base pipeline first.";
                return result;
            }
        }

        // Use pre-created run directory if adapter already prepared it
        File runDir;
        if (req.preCreatedRunDir != null && new File(req.preCreatedRunDir).isDirectory()) {
            runDir = new File(req.preCreatedRunDir);
            result.jobId = runDir.getName();
        } else {
            runDir = new File(analysisRoot, "runs/" + result.jobId);
            runDir.mkdirs();
        }
        result.runDir = runDir.getAbsolutePath();

        if (progress != null) progress.update("Preparing inputs", 0, 3);

        // Write input manifest
        InputContractWriter.DatasetMeta meta = new InputContractWriter.DatasetMeta();
        meta.projectId = req.projectId;
        meta.locusId = locus.id;
        meta.locusIndex = locus.index;
        meta.chr = locus.chr;
        meta.start = locus.start;
        meta.end = locus.end;
        meta.ancestry = config.refPanelPopulation;
        InputContractWriter.write(runDir, harmonizedDir, ldDir, meta);

        // Resolve command template
        String cmd = resolveCommand(td, req.params, config, locus,
            analysisRoot, harmonizedDir, ldDir, matchedDir, runDir);

        if (progress != null) progress.update("Running " + td.label, 1, 3);

        // Execute
        System.out.printf("[PluginEngine] Running %s: %s%n", td.tool, cmd);

        File logFile = new File(runDir, "run.log");
        int exitCode;
        try {
            exitCode = executeCommand(td.language, cmd, runDir, logFile);
        } catch (Exception e) {
            result.error = "Execution failed: " + e.getMessage();
            return result;
        }

        if (exitCode != 0) {
            result.error = "Tool exited with code " + exitCode + ". Check " + logFile.getAbsolutePath();
            return result;
        }

        if (progress != null) progress.update("Mapping output", 2, 3);

        // Map raw output to unified contract
        if (!td.outputMapping.isEmpty()) {
            mapOutput(td, runDir);
        }

        // Validate output
        OutputContractValidator.ValidationResult vr = OutputContractValidator.validate(runDir);
        result.validation = vr;
        if (!vr.valid) {
            result.error = "Output validation failed: " + String.join("; ", vr.errors);
            return result;
        }
        result.resultRows = vr.rowCount;

        // Write provenance
        writeProvenance(runDir, td, req.params, config, result.jobId);

        result.ok = true;
        if (progress != null) progress.update("Complete", 3, 3);

        System.out.printf("[PluginEngine] %s completed: %d result rows%n",
            td.tool, vr.rowCount);
        return result;
    }

    private static String resolveCommand(ToolDescriptor td, Map<String, String> params,
                                          Config config, Locus locus,
                                          File analysisRoot, File harmonizedDir,
                                          File ldDir, File matchedDir, File runDir) {
        String cmd = td.command;

        // Built-in placeholders
        cmd = cmd.replace("{run_dir}", runDir.getAbsolutePath());
        cmd = cmd.replace("{analysis_dir}", analysisRoot.getAbsolutePath());
        cmd = cmd.replace("{harmonized_gwas}", new File(harmonizedDir, "harmonized_gwas.tsv").getAbsolutePath());
        cmd = cmd.replace("{matched_ref}", new File(matchedDir, "matched_ref").getAbsolutePath());
        cmd = cmd.replace("{ld_r_matrix}", new File(ldDir, "ld_r.matrix").getAbsolutePath());
        cmd = cmd.replace("{ld_r2_matrix}", new File(ldDir, "ld_r2.matrix").getAbsolutePath());
        cmd = cmd.replace("{ld_snp_order}", new File(ldDir, "ld_snp_order.txt").getAbsolutePath());
        cmd = cmd.replace("{chr}", locus.chr);
        cmd = cmd.replace("{start}", String.valueOf(locus.start));
        cmd = cmd.replace("{end}", String.valueOf(locus.end));
        cmd = cmd.replace("{locus_id}", locus.id);
        cmd = cmd.replace("{ref_panel}", config.refPanelPath);
        cmd = cmd.replace("{plink_path}", PlinkSubsetter.findPlink(config) != null ?
            PlinkSubsetter.findPlink(config) : "plink");
        cmd = cmd.replace("{n_samples}", String.valueOf(0)); // from dataset meta

        // User parameters
        for (ToolDescriptor.Param p : td.params) {
            String val = params.getOrDefault(p.name, p.defaultValue);
            if (val != null) cmd = cmd.replace("{" + p.name + "}", val);
        }

        // Global config paths (e.g. {gcta_path})
        // These come from global config or user params
        for (Map.Entry<String, String> e : params.entrySet()) {
            cmd = cmd.replace("{" + e.getKey() + "}", e.getValue());
        }

        return cmd;
    }

    private static int executeCommand(String language, String command,
                                       File workDir, File logFile) throws Exception {
        List<String> cmdList = new ArrayList<>();

        switch (language != null ? language : "binary") {
            case "R":
                cmdList.add("Rscript");
                // If command is source('path/to/script.R'), extract the path and run directly
                if (command.startsWith("source('") && command.endsWith("')")) {
                    String scriptPath = command.substring(8, command.length() - 2);
                    cmdList.clear();
                    cmdList.add("Rscript");
                    cmdList.add(scriptPath);
                } else {
                    cmdList.add("-e");
                    cmdList.add(command);
                }
                break;
            case "python":
                cmdList.add("python");
                cmdList.add("-c");
                cmdList.add(command);
                break;
            case "binary":
            default:
                String os = System.getProperty("os.name", "").toLowerCase();
                if (os.contains("win")) {
                    cmdList.add("cmd.exe");
                    cmdList.add("/c");
                    cmdList.add(command);
                } else {
                    cmdList.add("/bin/sh");
                    cmdList.add("-c");
                    cmdList.add(command);
                }
                break;
        }

        ProcessBuilder pb = new ProcessBuilder(cmdList)
            .directory(workDir)
            .redirectErrorStream(true);

        // Add bin/ to PATH so DLLs (MKL, zlib) are found by child processes
        File binDir = new File("bin");
        if (binDir.isDirectory()) {
            Map<String, String> env = pb.environment();
            String path = env.getOrDefault("PATH", env.getOrDefault("Path", ""));
            env.put("PATH", binDir.getAbsolutePath() + File.pathSeparator + path);
        }

        Process proc = pb.start();

        // Capture stdout+stderr to log file and drain
        try (BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()));
             PrintWriter log = new PrintWriter(new BufferedWriter(new FileWriter(logFile)))) {
            String line;
            while ((line = br.readLine()) != null) {
                log.println(line);
            }
        }

        boolean finished = proc.waitFor(30, TimeUnit.MINUTES);
        if (!finished) {
            proc.destroyForcibly();
            throw new RuntimeException("Tool timed out after 30 minutes");
        }

        return proc.exitValue();
    }

    /**
     * Map raw tool output to unified result.tsv using output_mapping.
     */
    private static void mapOutput(ToolDescriptor td, File runDir) throws IOException {
        // Look for the raw output file (tool writes to run_dir)
        File rawOutput = null;
        for (File f : runDir.listFiles()) {
            if (f.getName().endsWith(".tsv") || f.getName().endsWith(".txt") ||
                f.getName().endsWith(".csv") || f.getName().endsWith(".jma.cojo") ||
                f.getName().endsWith(".cma.cojo")) {
                if (!f.getName().equals("result.tsv") && !f.getName().startsWith("run_input")) {
                    rawOutput = f;
                    break;
                }
            }
        }
        if (rawOutput == null) return;

        // Read raw output and map columns
        try (BufferedReader br = new BufferedReader(new FileReader(rawOutput))) {
            String header = br.readLine();
            if (header == null) return;
            String[] rawHeaders = header.split("\t");

            // Find column indices for mapping
            Map<Integer, ToolDescriptor.OutputColumn> colMap = new LinkedHashMap<>();
            int keyIdx = -1;
            for (int i = 0; i < rawHeaders.length; i++) {
                String h = rawHeaders[i].trim();
                if (h.equals(td.outputKey) || h.equals("snp_id") || h.equals("SNP")) {
                    keyIdx = i;
                }
                for (ToolDescriptor.OutputColumn oc : td.outputMapping) {
                    if (h.equals(oc.raw)) colMap.put(i, oc);
                }
            }

            if (keyIdx < 0 || colMap.isEmpty()) return;

            // Write mapped result.tsv
            File resultFile = new File(runDir, "result.tsv");
            try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(resultFile)))) {
                // Header
                StringBuilder hdr = new StringBuilder("snp_id");
                for (ToolDescriptor.OutputColumn oc : td.outputMapping) {
                    hdr.append('\t').append(oc.as);
                }
                pw.println(hdr);

                // Rows
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.trim().isEmpty()) continue;
                    String[] fields = line.split("\t", -1);
                    if (fields.length <= keyIdx) continue;

                    StringBuilder row = new StringBuilder(fields[keyIdx].trim());
                    for (ToolDescriptor.OutputColumn oc : td.outputMapping) {
                        row.append('\t');
                        Integer idx = null;
                        for (Map.Entry<Integer, ToolDescriptor.OutputColumn> e : colMap.entrySet()) {
                            if (e.getValue() == oc) { idx = e.getKey(); break; }
                        }
                        if (idx != null && idx < fields.length)
                            row.append(fields[idx].trim());
                        else
                            row.append("NA");
                    }
                    pw.println(row);
                }
            }

            // Write result manifest
            List<OutputContractValidator.ColumnDecl> cols = new ArrayList<>();
            for (ToolDescriptor.OutputColumn oc : td.outputMapping) {
                OutputContractValidator.ColumnDecl cd = new OutputContractValidator.ColumnDecl();
                cd.name = oc.as;
                cd.type = oc.type != null ? oc.type : "double";
                cd.scope = oc.scope != null ? oc.scope : "per_snp";
                cd.method = td.tool;
                cd.methodVersion = td.version;
                cols.add(cd);
            }
            OutputContractValidator.writeManifest(runDir, td.tool, td.version,
                new LinkedHashMap<>(), cols);
        }
    }

    private static void writeProvenance(File runDir, ToolDescriptor td,
                                         Map<String, String> params,
                                         Config config, String jobId) throws IOException {
        File provFile = new File(runDir, "provenance.json");
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(provFile)))) {
            pw.println("{");
            pw.printf("  \"job_id\": \"%s\",%n", esc(jobId));
            pw.printf("  \"tool\": \"%s\",%n", esc(td.tool));
            pw.printf("  \"tool_version\": \"%s\",%n", esc(td.version));
            pw.printf("  \"language\": \"%s\",%n", esc(td.language));
            pw.printf("  \"ref_panel\": \"%s\",%n", esc(config.refPanelPath));
            pw.printf("  \"ref_panel_population\": \"%s\",%n", esc(config.refPanelPopulation));
            pw.println("  \"parameters\": {");
            int i = 0;
            for (Map.Entry<String, String> e : params.entrySet()) {
                pw.printf("    \"%s\": \"%s\"%s%n",
                    esc(e.getKey()), esc(e.getValue()), (++i < params.size()) ? "," : "");
            }
            pw.println("  },");
            pw.printf("  \"timestamp\": %d%n", System.currentTimeMillis());
            pw.println("}");
        }
    }

    private static String esc(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
