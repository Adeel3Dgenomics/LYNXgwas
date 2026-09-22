@echo off
echo === Running LYNXgwas tests ===

set "LIBCP="
for %%f in (lib\*.jar) do call set "LIBCP=%%LIBCP%%;%%f"

javac -encoding UTF-8 -d bin -cp "bin%LIBCP%" ^
  src\analysis\StatsUtil.java ^
  src\analysis\AnovaUtil.java ^
  src\analysis\EnrichmentAnalyzer.java ^
  src\analysis\RegulatoryPeakIndex.java ^
  src\analysis\RegulatoryEnrichmentAnalyzer.java ^
  src\analysis\MiniJson.java ^
  src\analysis\LlmClient.java ^
  src\analysis\AgentToolRegistry.java ^
  src\analysis\AgentConfig.java ^
  src\analysis\AgentOrchestrator.java ^
  src\analysis\OllamaPuller.java ^
  tests\MultiLocusScannerTest.java ^
  tests\LdCalculatorTest.java ^
  tests\GwasQcTest.java ^
  tests\GlobalConfigTest.java ^
  tests\SharedStorageResolverTest.java ^
  tests\GithubTokenStoreTest.java ^
  tests\GithubProjectSyncTest.java ^
  tests\SnakemakeSubmitterTest.java ^
  tests\AnovaUtilTest.java ^
  tests\MagmaAdapterTest.java ^
  tests\GctaGremlAdapterTest.java ^
  tests\SusieAdapterTest.java ^
  tests\FinemapAdapterTest.java ^
  tests\CojoAdapterTest.java ^
  tests\ColocAdapterTest.java ^
  tests\GwamaAdapterTest.java ^
  tests\EnrichmentAnalyzerTest.java ^
  tests\GeneConstellationBuilderTest.java ^
  tests\RegionConstellationBuilderTest.java ^
  tests\GlobalSearchIndexTest.java ^
  tests\RegulatoryPeakIndexTest.java ^
  tests\RegulatoryEnrichmentAnalyzerTest.java ^
  tests\MiniJsonTest.java ^
  tests\AgentToolRegistryTest.java ^
  tests\OllamaPullerTest.java
if %ERRORLEVEL% neq 0 (
    echo [FAIL] Test compilation failed.
    exit /b 1
)

set FAILED=0

echo.
echo --- MultiLocusScannerTest ---
java -cp "bin%LIBCP%" MultiLocusScannerTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- LdCalculatorTest ---
java -cp "bin%LIBCP%" LdCalculatorTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- GwasQcTest ---
java -cp "bin%LIBCP%" GwasQcTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- GlobalConfigTest ---
java -cp "bin%LIBCP%" GlobalConfigTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- SharedStorageResolverTest ---
java -cp "bin%LIBCP%" SharedStorageResolverTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- GithubTokenStoreTest ---
java -cp "bin%LIBCP%" GithubTokenStoreTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- GithubProjectSyncTest ---
java -cp "bin%LIBCP%" GithubProjectSyncTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- SnakemakeSubmitterTest ---
java -cp "bin%LIBCP%" SnakemakeSubmitterTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- AnovaUtilTest ---
java -cp "bin%LIBCP%" AnovaUtilTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- MagmaAdapterTest ---
java -cp "bin%LIBCP%" MagmaAdapterTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- GctaGremlAdapterTest ---
java -cp "bin%LIBCP%" GctaGremlAdapterTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- SusieAdapterTest ---
java -cp "bin%LIBCP%" SusieAdapterTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- FinemapAdapterTest ---
java -cp "bin%LIBCP%" FinemapAdapterTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- CojoAdapterTest ---
java -cp "bin%LIBCP%" CojoAdapterTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- ColocAdapterTest ---
java -cp "bin%LIBCP%" ColocAdapterTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- GwamaAdapterTest ---
java -cp "bin%LIBCP%" GwamaAdapterTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- EnrichmentAnalyzerTest ---
java -cp "bin%LIBCP%" EnrichmentAnalyzerTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- GeneConstellationBuilderTest ---
java -cp "bin%LIBCP%" GeneConstellationBuilderTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- RegionConstellationBuilderTest ---
java -cp "bin%LIBCP%" RegionConstellationBuilderTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- GlobalSearchIndexTest ---
java -cp "bin%LIBCP%" GlobalSearchIndexTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- RegulatoryPeakIndexTest ---
java -cp "bin%LIBCP%" RegulatoryPeakIndexTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- RegulatoryEnrichmentAnalyzerTest ---
java -cp "bin%LIBCP%" RegulatoryEnrichmentAnalyzerTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- MiniJsonTest ---
java -cp "bin%LIBCP%" MiniJsonTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- AgentToolRegistryTest ---
java -cp "bin%LIBCP%" AgentToolRegistryTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
echo --- OllamaPullerTest ---
java -cp "bin%LIBCP%" OllamaPullerTest
if %ERRORLEVEL% neq 0 set FAILED=1

echo.
if %FAILED%==0 (
    echo [OK] All test suites passed.
) else (
    echo [FAIL] One or more test suites failed.
    exit /b 1
)
