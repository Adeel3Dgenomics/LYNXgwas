# Graph Report - .  (2026-06-17)

## Corpus Check
- Corpus is ~21,657 words - fits in a single context window. You may not need a graph.

## Summary
- 237 nodes · 419 edges · 18 communities (12 shown, 6 thin omitted)
- Extraction: 90% EXTRACTED · 10% INFERRED · 0% AMBIGUOUS · INFERRED: 43 edges (avg confidence: 0.82)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Pipeline Orchestration|Pipeline Orchestration]]
- [[_COMMUNITY_JSON Serialization|JSON Serialization]]
- [[_COMMUNITY_Gene Annotation Parsing|Gene Annotation Parsing]]
- [[_COMMUNITY_LD Computation Engine|LD Computation Engine]]
- [[_COMMUNITY_HTTP Server Endpoints|HTTP Server Endpoints]]
- [[_COMMUNITY_GWAS Locus Parsing|GWAS Locus Parsing]]
- [[_COMMUNITY_Locus Visualization Panel|Locus Visualization Panel]]
- [[_COMMUNITY_LD Block Detection|LD Block Detection]]
- [[_COMMUNITY_Loci File Parser|Loci File Parser]]
- [[_COMMUNITY_SNP Annotation Service|SNP Annotation Service]]
- [[_COMMUNITY_Genome Skyline Generator|Genome Skyline Generator]]
- [[_COMMUNITY_Output Data Models|Output Data Models]]
- [[_COMMUNITY_Transcript Model|Transcript Model]]
- [[_COMMUNITY_Exon Model|Exon Model]]
- [[_COMMUNITY_SNP Model|SNP Model]]
- [[_COMMUNITY_HTML Document Roots|HTML Document Roots]]
- [[_COMMUNITY_App Initialization|App Initialization]]
- [[_COMMUNITY_Empty Document|Empty Document]]

## God Nodes (most connected - your core abstractions)
1. `LocalServer` - 19 edges
2. `GffParser` - 11 edges
3. `Jb` - 10 edges
4. `HttpExchange` - 10 edges
5. `String` - 9 edges
6. `String` - 8 edges
7. `LdCalculator` - 8 edges
8. `fillPanel (Locus Panel Composition)` - 8 edges
9. `JsonExporter` - 7 edges
10. `String` - 7 edges

## Surprising Connections (you probably didn't know these)
- `sanitizeAndParse (JSON Sparse Array Repair)` --semantically_similar_to--> `sanitizeJsonText (Main App JSON Repair)`  [INFERRED] [semantically similar]
  ld_triangle_test.html → output/index.html
- `r2Color (LD Triangle Test Color Mapping)` --semantically_similar_to--> `r2Color (Main App LD Color Mapping)`  [INFERRED] [semantically similar]
  ld_triangle_test.html → output/index.html
- `classifyPairs (Test Page LD Pair Classification)` --semantically_similar_to--> `classifyPairs (Main App LD Pair Classification)`  [INFERRED] [semantically similar]
  ld_triangle_test.html → output/index.html
- `detectLdBlocks (Test Page Gabriel-style Block Detection)` --semantically_similar_to--> `detectLdBlocks (Main App Gabriel-style Block Detection)`  [INFERRED] [semantically similar]
  ld_triangle_test.html → output/index.html
- `drawBlockOutlines (Test Page Block Outline Renderer)` --semantically_similar_to--> `drawBlockOutlines (Main App Block Outline Renderer)`  [INFERRED] [semantically similar]
  ld_triangle_test.html → output/index.html

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Double-Comma Sparse Array Repair Pipeline** — concept_double_comma_serialization_bug, lociCatalogue_ld_triangle_test_sanitizeandparse, output_index_sanitizejsontext, output_index_sanitizeldmatrix [EXTRACTED 1.00]
- **Gabriel-style LD Block Detection System** — concept_gabriel_ld_block_detection, lociCatalogue_ld_triangle_test_classifypairs, lociCatalogue_ld_triangle_test_detectldblocks, output_index_classifypairs, output_index_detectldblocks [EXTRACTED 1.00]
- **GWAS Locus Panel Multi-Track Visualization** — output_index_fillpanel, output_index_rendermanhattan, output_index_rendergenetrack, output_index_renderldtrack, output_index_renderldtriangle, output_index_rendercontextsketch, output_index_renderoverviewdensity [EXTRACTED 1.00]

## Communities (18 total, 6 thin omitted)

### Community 0 - "Pipeline Orchestration"
Cohesion: 0.09
Nodes (18): Config, String, List, LocusOutput, String, Main, Config, Integer (+10 more)

### Community 1 - "JSON Serialization"
Cohesion: 0.19
Nodes (13): Jb, LdTriangle, LocusRef, Runnable, Config, Exon, List, Locus (+5 more)

### Community 2 - "Gene Annotation Parsing"
Cohesion: 0.15
Nodes (10): Gene, Gene, String, GffParser, Config, Exon, List, Map (+2 more)

### Community 3 - "LD Computation Engine"
Cohesion: 0.18
Nodes (12): File, Config, Integer, List, Locus, Map, ProgressTracker, Snp (+4 more)

### Community 4 - "HTTP Server Endpoints"
Cohesion: 0.28
Nodes (4): HttpExchange, InputStream, String, LocalServer

### Community 5 - "GWAS Locus Parsing"
Cohesion: 0.17
Nodes (8): GwasParser, Config, List, Locus, String, Override, String, Locus

### Community 6 - "Locus Visualization Panel"
Cohesion: 0.14
Nodes (16): Double-Comma Sparse Array Serialization Bug, Multi-Track GWAS Locus Panel Architecture, extractLocusObject (Locus JSON/JS Loader), sanitizeAndParse (JSON Sparse Array Repair), buildExonIntervals (Exon Interval Merging), fillPanel (Locus Panel Composition), isExonic (Binary Search Exonic SNP Check), loadLocus (Locus Data Loader with Sanitization) (+8 more)

### Community 7 - "LD Block Detection"
Cohesion: 0.16
Nodes (16): Gabriel-style Proportion-based LD Block Detection, Haploview-style Color Palette for LD Visualization, Reference Panel vs GWAS SNP Filtering, classifyPairs (Test Page LD Pair Classification), detectLdBlocks (Test Page Gabriel-style Block Detection), drawBlockOutlines (Test Page Block Outline Renderer), drawTriangle (Test Page Triangle Canvas Renderer), r2Color (LD Triangle Test Color Mapping) (+8 more)

### Community 8 - "Loci File Parser"
Cohesion: 0.24
Nodes (7): Config, Integer, List, Locus, Map, String, LociParser

### Community 9 - "SNP Annotation Service"
Cohesion: 0.27
Nodes (6): Config, Integer, Map, Snp, String, SnpAnnotator

### Community 10 - "Genome Skyline Generator"
Cohesion: 0.43
Nodes (3): GenomeSkyline, Config, String

### Community 11 - "Output Data Models"
Cohesion: 0.25
Nodes (6): String, LdSnp, LdTriangle, LocusContext, LocusOutput, LocusRef

## Knowledge Gaps
- **59 isolated node(s):** `String`, `String`, `Transcript`, `Config`, `Config` (+54 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **6 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `LocalServer` connect `HTTP Server Endpoints` to `LD Computation Engine`?**
  _High betweenness centrality (0.112) - this node is a cross-community bridge._
- **Why does `Locus` connect `GWAS Locus Parsing` to `Pipeline Orchestration`?**
  _High betweenness centrality (0.098) - this node is a cross-community bridge._
- **What connects `String`, `String`, `Transcript` to the rest of the system?**
  _61 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Pipeline Orchestration` be split into smaller, more focused modules?**
  _Cohesion score 0.0907563025210084 - nodes in this community are weakly interconnected._
- **Should `Locus Visualization Panel` be split into smaller, more focused modules?**
  _Cohesion score 0.14166666666666666 - nodes in this community are weakly interconnected._