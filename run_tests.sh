#!/bin/sh
# Cross-platform (Linux/macOS) equivalent of run_tests.bat.
set -e
echo "=== Running LYNXgwas tests ==="

LIBCP=""
for f in lib/*.jar; do
  [ -e "$f" ] && LIBCP="$LIBCP:$f"
done

javac -encoding UTF-8 -d bin -cp "bin$LIBCP" \
  src/analysis/StatsUtil.java \
  src/analysis/AnovaUtil.java \
  src/analysis/EnrichmentAnalyzer.java \
  src/analysis/RegulatoryPeakIndex.java \
  src/analysis/RegulatoryEnrichmentAnalyzer.java \
  src/analysis/MiniJson.java \
  src/analysis/LlmClient.java \
  src/analysis/AgentToolRegistry.java \
  src/analysis/AgentConfig.java \
  src/analysis/AgentOrchestrator.java \
  src/analysis/OllamaPuller.java \
  tests/MultiLocusScannerTest.java \
  tests/LdCalculatorTest.java \
  tests/GwasQcTest.java \
  tests/GlobalConfigTest.java \
  tests/SharedStorageResolverTest.java \
  tests/GithubTokenStoreTest.java \
  tests/GithubProjectSyncTest.java \
  tests/SnakemakeSubmitterTest.java \
  tests/AnovaUtilTest.java \
  tests/MagmaAdapterTest.java \
  tests/GctaGremlAdapterTest.java \
  tests/SusieAdapterTest.java \
  tests/FinemapAdapterTest.java \
  tests/CojoAdapterTest.java \
  tests/ColocAdapterTest.java \
  tests/GwamaAdapterTest.java \
  tests/EnrichmentAnalyzerTest.java \
  tests/GeneConstellationBuilderTest.java \
  tests/RegionConstellationBuilderTest.java \
  tests/GlobalSearchIndexTest.java \
  tests/RegulatoryPeakIndexTest.java \
  tests/RegulatoryEnrichmentAnalyzerTest.java \
  tests/MiniJsonTest.java \
  tests/AgentToolRegistryTest.java \
  tests/OllamaPullerTest.java

FAILED=0

echo
echo "--- MultiLocusScannerTest ---"
java -cp "bin$LIBCP" MultiLocusScannerTest || FAILED=1

echo
echo "--- LdCalculatorTest ---"
java -cp "bin$LIBCP" LdCalculatorTest || FAILED=1

echo
echo "--- GwasQcTest ---"
java -cp "bin$LIBCP" GwasQcTest || FAILED=1

echo
echo "--- GlobalConfigTest ---"
java -cp "bin$LIBCP" GlobalConfigTest || FAILED=1

echo
echo "--- SharedStorageResolverTest ---"
java -cp "bin$LIBCP" SharedStorageResolverTest || FAILED=1

echo
echo "--- GithubTokenStoreTest ---"
java -cp "bin$LIBCP" GithubTokenStoreTest || FAILED=1

echo
echo "--- GithubProjectSyncTest ---"
java -cp "bin$LIBCP" GithubProjectSyncTest || FAILED=1

echo
echo "--- SnakemakeSubmitterTest ---"
java -cp "bin$LIBCP" SnakemakeSubmitterTest || FAILED=1

echo
echo "--- AnovaUtilTest ---"
java -cp "bin$LIBCP" AnovaUtilTest || FAILED=1

echo
echo "--- MagmaAdapterTest ---"
java -cp "bin$LIBCP" MagmaAdapterTest || FAILED=1

echo
echo "--- GctaGremlAdapterTest ---"
java -cp "bin$LIBCP" GctaGremlAdapterTest || FAILED=1

echo
echo "--- SusieAdapterTest ---"
java -cp "bin$LIBCP" SusieAdapterTest || FAILED=1

echo
echo "--- FinemapAdapterTest ---"
java -cp "bin$LIBCP" FinemapAdapterTest || FAILED=1

echo
echo "--- CojoAdapterTest ---"
java -cp "bin$LIBCP" CojoAdapterTest || FAILED=1

echo
echo "--- ColocAdapterTest ---"
java -cp "bin$LIBCP" ColocAdapterTest || FAILED=1

echo
echo "--- GwamaAdapterTest ---"
java -cp "bin$LIBCP" GwamaAdapterTest || FAILED=1

echo
echo "--- EnrichmentAnalyzerTest ---"
java -cp "bin$LIBCP" EnrichmentAnalyzerTest || FAILED=1

echo
echo "--- GeneConstellationBuilderTest ---"
java -cp "bin$LIBCP" GeneConstellationBuilderTest || FAILED=1

echo
echo "--- RegionConstellationBuilderTest ---"
java -cp "bin$LIBCP" RegionConstellationBuilderTest || FAILED=1

echo
echo "--- GlobalSearchIndexTest ---"
java -cp "bin$LIBCP" GlobalSearchIndexTest || FAILED=1

echo
echo "--- RegulatoryPeakIndexTest ---"
java -cp "bin$LIBCP" RegulatoryPeakIndexTest || FAILED=1

echo
echo "--- RegulatoryEnrichmentAnalyzerTest ---"
java -cp "bin$LIBCP" RegulatoryEnrichmentAnalyzerTest || FAILED=1

echo
echo "--- MiniJsonTest ---"
java -cp "bin$LIBCP" MiniJsonTest || FAILED=1

echo
echo "--- AgentToolRegistryTest ---"
java -cp "bin$LIBCP" AgentToolRegistryTest || FAILED=1

echo
echo "--- OllamaPullerTest ---"
java -cp "bin$LIBCP" OllamaPullerTest || FAILED=1

echo
if [ "$FAILED" -eq 0 ]; then
  echo "[OK] All test suites passed."
else
  echo "[FAIL] One or more test suites failed."
  exit 1
fi
