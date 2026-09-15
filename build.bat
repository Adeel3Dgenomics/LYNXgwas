@echo off
echo === Building LYNXgwas ===
if not exist bin mkdir bin

REM Build classpath for HTSJDK dependencies
set "LIBCP="
for %%f in (lib\*.jar) do call set "LIBCP=%%LIBCP%%;%%f"

javac -encoding UTF-8 -d bin -cp "bin%LIBCP%" ^
  src\rsid\DbSnpRecord.java ^
  src\rsid\MatchResult.java ^
  src\rsid\RsidMatcher.java ^
  src\rsid\RsidRecovery.java ^
  src\rsid\RsidDetector.java ^
  src\rsid\GlobalConfig.java ^
  src\rsid\RsidProgress.java ^
  src\rsid\RateLimiter.java ^
  src\rsid\RsidApiProvider.java ^
  src\rsid\RsidApiCache.java ^
  src\rsid\NcbiDbSnpProvider.java ^
  src\rsid\GnomadProvider.java ^
  src\rsid\RsidApiCompleter.java ^
  src\rsid\CrossFileLookup.java ^
  src\rsid\RsidPipeline.java ^
  src\loci\LociIdentifier.java ^
  src\loci\LociProgress.java ^
  src\analysis\ContentHasher.java ^
  src\analysis\StatsUtil.java ^
  src\analysis\StepManifest.java ^
  src\analysis\LocusGwasExtractor.java ^
  src\analysis\SnpMatcher.java ^
  src\analysis\AlleleHarmonizer.java ^
  src\analysis\LdMatrixComputer.java ^
  src\analysis\LdGwasDiagnostic.java ^
  src\analysis\StableSnpId.java ^
  src\analysis\InputContractWriter.java ^
  src\analysis\OutputContractValidator.java ^
  src\analysis\ToolDescriptor.java ^
  src\analysis\PluginEngine.java ^
  src\analysis\CojoAdapter.java ^
  src\analysis\SusieAdapter.java ^
  src\analysis\FinemapAdapter.java ^
  src\analysis\ColocAdapter.java ^
  src\analysis\GwamaAdapter.java ^
  src\analysis\SusiexAdapter.java ^
  src\analysis\AnalysisColumnProvider.java ^
  src\analysis\BaseStepPipeline.java ^
  src\export\ColumnSpec.java ^
  src\export\GwasSchema.java ^
  src\export\SnpContext.java ^
  src\export\LocusContext.java ^
  src\export\SnpColumnProvider.java ^
  src\export\LocusColumnProvider.java ^
  src\export\ExportRegistry.java ^
  src\export\XlsxWriter.java ^
  src\export\ExcelExporter.java ^
  src\Config.java ^
  src\Locus.java ^
  src\Snp.java ^
  src\Exon.java ^
  src\Transcript.java ^
  src\Gene.java ^
  src\LocusOutput.java ^
  src\ProgressTracker.java ^
  src\LociParser.java ^
  src\GwasParser.java ^
  src\GwasQc.java ^
  src\GffParser.java ^
  src\MultiLocusProgress.java ^
  src\MultiLocusMerger.java ^
  src\MultiLocusScanner.java ^
  src\MultiLocusResult.java ^
  src\MultiLocusExcelWriter.java ^
  src\PlinkSubsetter.java ^
  src\LdCalculator.java ^
  src\SnpAnnotator.java ^
  src\GenomeSkyline.java ^
  src\ProjectMetadata.java ^
  src\JsonExporter.java ^
  src\LocusUpdater.java ^
  src\LociMutationService.java ^
  src\LocalServer.java ^
  src\Main.java
if %ERRORLEVEL% neq 0 (
    echo [FAIL] Compilation failed.
    exit /b 1
)
echo [OK] Compiled to bin\

echo Copying frontend assets to output\...
if not exist output mkdir output
if not exist output\assets mkdir output\assets

copy /Y index.html output\index.html >nul 2>&1
copy /Y viewer.html output\viewer.html >nul 2>&1
copy /Y locus_matrix.html output\locus_matrix.html >nul 2>&1
copy /Y annotations.js output\annotations.js >nul 2>&1
if exist output\annotations.js (
    rem already copied from root
) else (
    rem annotations.js may only exist in output/ already
)
xcopy /Y /E /I assets output\assets >nul 2>&1
if exist projects\config.properties.template (
    if not exist output\projects mkdir output\projects
    copy /Y projects\config.properties.template output\projects\config.properties.template >nul 2>&1
)

echo [OK] Build complete.
