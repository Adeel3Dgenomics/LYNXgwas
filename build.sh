#!/bin/sh
# Cross-platform (Linux/macOS) equivalent of build.bat. Mirrors it file-for-file so the two
# stay in lockstep — if you add a source file to one, add it to the other.
set -e
echo "=== Building LYNXgwas ==="
mkdir -p bin

LIBCP=""
for f in lib/*.jar; do
  [ -e "$f" ] && LIBCP="$LIBCP:$f"
done

javac -encoding UTF-8 -d bin -cp "bin$LIBCP" \
  src/rsid/DbSnpRecord.java \
  src/rsid/MatchResult.java \
  src/rsid/RsidMatcher.java \
  src/rsid/RsidRecovery.java \
  src/rsid/RsidDetector.java \
  src/rsid/GlobalConfig.java \
  src/rsid/RsidProgress.java \
  src/rsid/RateLimiter.java \
  src/rsid/RsidApiProvider.java \
  src/rsid/RsidApiCache.java \
  src/rsid/NcbiDbSnpProvider.java \
  src/rsid/GnomadProvider.java \
  src/rsid/RsidApiCompleter.java \
  src/rsid/CrossFileLookup.java \
  src/rsid/RsidPipeline.java \
  src/rsid/SharedStorageResolver.java \
  src/catalog/GwasCatalogClient.java \
  src/catalog/GwasCatalogLocalIndex.java \
  src/opentargets/OpenTargetsL2GClient.java \
  src/loci/LociIdentifier.java \
  src/loci/LociProgress.java \
  src/analysis/ContentHasher.java \
  src/analysis/StatsUtil.java \
  src/analysis/AnovaUtil.java \
  src/analysis/GeneConstellationResult.java \
  src/analysis/GeneConstellationBuilder.java \
  src/analysis/GlobalSearchIndex.java \
  src/analysis/StepManifest.java \
  src/analysis/LocusGwasExtractor.java \
  src/analysis/SnpMatcher.java \
  src/analysis/AlleleHarmonizer.java \
  src/analysis/LdMatrixComputer.java \
  src/analysis/LdGwasDiagnostic.java \
  src/analysis/StableSnpId.java \
  src/analysis/InputContractWriter.java \
  src/analysis/OutputContractValidator.java \
  src/analysis/ToolDescriptor.java \
  src/analysis/PluginEngine.java \
  src/analysis/SnakemakeSubmitter.java \
  src/analysis/GctaBinaryResolver.java \
  src/analysis/CojoAdapter.java \
  src/analysis/GctaGremlAdapter.java \
  src/analysis/MagmaAdapter.java \
  src/analysis/SusieAdapter.java \
  src/analysis/FinemapAdapter.java \
  src/analysis/ColocAdapter.java \
  src/analysis/GwamaAdapter.java \
  src/analysis/SusiexAdapter.java \
  src/analysis/AnalysisColumnProvider.java \
  src/analysis/BaseStepPipeline.java \
  src/analysis/EnrichmentAnalyzer.java \
  src/analysis/RegulatoryPeakIndex.java \
  src/analysis/RegulatoryEnrichmentAnalyzer.java \
  src/export/ColumnSpec.java \
  src/export/GwasSchema.java \
  src/export/SnpContext.java \
  src/export/LocusContext.java \
  src/export/SnpColumnProvider.java \
  src/export/LocusColumnProvider.java \
  src/export/ExportRegistry.java \
  src/export/XlsxWriter.java \
  src/export/ExcelExporter.java \
  src/Config.java \
  src/Locus.java \
  src/Snp.java \
  src/Exon.java \
  src/Transcript.java \
  src/Gene.java \
  src/LocusOutput.java \
  src/ProgressTracker.java \
  src/LociParser.java \
  src/GwasParser.java \
  src/GwasQc.java \
  src/GffParser.java \
  src/MultiLocusProgress.java \
  src/MultiLocusMerger.java \
  src/MultiLocusScanner.java \
  src/MultiLocusResult.java \
  src/MultiLocusExcelWriter.java \
  src/GeneConstellationExcelWriter.java \
  src/PlinkSubsetter.java \
  src/LdCalculator.java \
  src/SnpAnnotator.java \
  src/GenomeSkyline.java \
  src/GenomeLiftover.java \
  src/ProjectMetadata.java \
  src/JsonExporter.java \
  src/LocusUpdater.java \
  src/LociMutationService.java \
  src/LocalServer.java \
  src/Main.java

echo "[OK] Compiled to bin/"

echo "Copying frontend assets to output/..."
mkdir -p output/assets

cp -f index.html output/index.html 2>/dev/null || true
cp -f viewer.html output/viewer.html 2>/dev/null || true
cp -f gene_constellation.html output/gene_constellation.html 2>/dev/null || true
cp -f annotations.js output/annotations.js 2>/dev/null || true
cp -Rf assets/. output/assets/ 2>/dev/null || true
if [ -f projects/config.properties.template ]; then
  mkdir -p output/projects
  cp -f projects/config.properties.template output/projects/config.properties.template
fi

echo "[OK] Build complete."
