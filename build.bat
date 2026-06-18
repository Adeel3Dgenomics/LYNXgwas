@echo off
echo === Building GWAS Locus Visualization Tool v2 ===
if not exist bin mkdir bin
javac -d bin ^
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
  src\GffParser.java ^
  src\PlinkSubsetter.java ^
  src\LdCalculator.java ^
  src\SnpAnnotator.java ^
  src\GenomeSkyline.java ^
  src\JsonExporter.java ^
  src\LocusUpdater.java ^
  src\LocalServer.java ^
  src\Main.java
if %ERRORLEVEL% neq 0 (
    echo [FAIL] Compilation failed.
    exit /b 1
)
echo [OK] Compiled to bin\
