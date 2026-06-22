# Graph Report - .  (2026-06-22)

## Corpus Check
- 60 files · ~123,554 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 805 nodes · 2158 edges · 28 communities (26 shown, 2 thin omitted)
- Extraction: 84% EXTRACTED · 16% INFERRED · 0% AMBIGUOUS · INFERRED: 337 edges (avg confidence: 0.8)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_LD & RSID Cache|LD & RSID Cache]]
- [[_COMMUNITY_Allele Matching & DbSNP|Allele Matching & DbSNP]]
- [[_COMMUNITY_Excel Export Engine|Excel Export Engine]]
- [[_COMMUNITY_RSID API Providers|RSID API Providers]]
- [[_COMMUNITY_Local Server HTTP|Local Server HTTP]]
- [[_COMMUNITY_Project Metadata|Project Metadata]]
- [[_COMMUNITY_GWAS Pipeline Config|GWAS Pipeline Config]]
- [[_COMMUNITY_Configuration Parsing|Configuration Parsing]]
- [[_COMMUNITY_Pipeline Orchestration|Pipeline Orchestration]]
- [[_COMMUNITY_Gene Model & GFF|Gene Model & GFF]]
- [[_COMMUNITY_Web Frontend & Docs|Web Frontend & Docs]]
- [[_COMMUNITY_Locus Mutation Service|Locus Mutation Service]]
- [[_COMMUNITY_Loci Identification|Loci Identification]]
- [[_COMMUNITY_SNP Annotation|SNP Annotation]]
- [[_COMMUNITY_JSON Export & LD Triangle|JSON Export & LD Triangle]]
- [[_COMMUNITY_SNP Export Providers|SNP Export Providers]]
- [[_COMMUNITY_GWAS File Parser|GWAS File Parser]]
- [[_COMMUNITY_Server Endpoint Tests|Server Endpoint Tests]]
- [[_COMMUNITY_Export Column Spec|Export Column Spec]]
- [[_COMMUNITY_Provider Value Methods|Provider Value Methods]]
- [[_COMMUNITY_Multi-Project Tests|Multi-Project Tests]]
- [[_COMMUNITY_Locus Export Providers|Locus Export Providers]]
- [[_COMMUNITY_Locus Column Provider|Locus Column Provider]]
- [[_COMMUNITY_SNP Column Provider|SNP Column Provider]]
- [[_COMMUNITY_RSID Detection|RSID Detection]]
- [[_COMMUNITY_Logo & Branding|Logo & Branding]]
- [[_COMMUNITY_Locus Counts Provider|Locus Counts Provider]]
- [[_COMMUNITY_Build Checks Config|Build Checks Config]]

## God Nodes (most connected - your core abstractions)
1. `LocalServer` - 57 edges
2. `HttpExchange` - 41 edges
3. `String` - 23 edges
4. `LociMutationService` - 20 edges
5. `ColumnSpec` - 20 edges
6. `GffParser` - 18 edges
7. `ProjectMetadata` - 18 edges
8. `ExcelExporter` - 17 edges
9. `XlsxWriter` - 17 edges
10. `String` - 16 edges

## Surprising Connections (you probably didn't know these)
- `LYNXgwas Icon SVG` --conceptually_related_to--> `LYNXgwas Home Page`  [INFERRED]
  assets/lynxgwas-icon.svg → index.html
- `LYNXgwas Home Page` --references--> `Loci Identification`  [EXTRACTED]
  index.html → README.md
- `LYNXgwas Home Page` --references--> `LYNXgwas Interactive Viewer`  [EXTRACTED]
  index.html → viewer.html
- `LYNXgwas Home Page` --references--> `rsID Recovery Pipeline`  [EXTRACTED]
  index.html → README.md
- `LYNXgwas Interactive Viewer` --references--> `Multi-Project GWAS Management`  [EXTRACTED]
  viewer.html → README.md

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Core GWAS Pipeline Flow** — src_main, src_lociparser, src_gffparser, src_gwasparser, src_plinksubsetter, src_ldcalculator, src_jsonexporter [EXTRACTED 1.00]
- **Locus Mutation System (split/merge/resize/undo)** — src_locimutationservice, src_locusupdater, src_localserver, src_locus, src_locusoutput [EXTRACTED 1.00]
- **Gene Model Hierarchy (Gene > Transcript > Exon)** — src_gene, src_gene_transcript, src_exon, src_gffparser [EXTRACTED 1.00]
- **Excel Export Subsystem** — export_excelexporter, export_exportregistry, export_xlsxwriter, export_snpcolumnprovider, export_locuscolumnprovider, export_columnspec, export_snpcontext, export_locuscontext [EXTRACTED 1.00]
- **Loci Identification Pipeline** — loci_lociidentifier, loci_lociprogress, concept_plink_clumping [EXTRACTED 1.00]
- **Project Lifecycle Test Suite** — src_testfingerprint, src_testmultiproject, src_testperprojectconfig, src_testprojectoutput, src_testserverendpoints, src_testintegration, src_testgwascolumnmapping [INFERRED 0.95]
- **rsID Recovery Pipeline Flow** — rsid_rsidpipeline, rsid_rsidrecovery, rsid_rsidmatcher, rsid_rsidapicompleter, rsid_rsidapicache, rsid_rsidprogress, rsid_dbsnprecord, rsid_matchresult [EXTRACTED 1.00]
- **External API Provider Strategy Pattern** — rsid_rsidapiprovider, rsid_gnomadprovider, rsid_ncbidbsnpprovider, rsid_ratelimiter, rsid_rsidapicompleter [EXTRACTED 1.00]
- **LYNXgwas Web Frontend** — lynxgwas_index, lynxgwas_viewer, assets_lynxgwas_icon [EXTRACTED 1.00]

## Communities (28 total, 2 thin omitted)

### Community 0 - "LD & RSID Cache"
Cohesion: 0.07
Nodes (35): Linkage Disequilibrium Computation, RsidApiCache, Config, File, Integer, List, Locus, Map (+27 more)

### Community 1 - "Allele Matching & DbSNP"
Cohesion: 0.07
Nodes (31): Forward/Reverse Allele Matching, DbSnpRecord, MatchResult, RsidMatcher, RsidPipeline, RsidRecovery, TestRsidRecovery, Integer (+23 more)

### Community 2 - "Excel Export Engine"
Cohesion: 0.08
Nodes (21): Closeable, ExcelExporter, ExportProgress, GenePos, ExportRegistry, Sheet, SheetData, XlsxWriter (+13 more)

### Community 3 - "RSID API Providers"
Cohesion: 0.06
Nodes (27): API-based rsID Completion, Token Bucket Rate Limiting, GnomadProvider, NcbiDbSnpProvider, RateLimiter, CompletionResult, Config, RsidApiCompleter (+19 more)

### Community 4 - "Local Server HTTP"
Cohesion: 0.18
Nodes (5): HttpExchange, InputStream, SplitRegion, String, LocalServer

### Community 5 - "Project Metadata"
Cohesion: 0.12
Nodes (14): Fingerprint-Based Staleness Detection, Config, File, List, Map, StaleReason, String, ProjectMetadata (+6 more)

### Community 6 - "GWAS Pipeline Config"
Cohesion: 0.08
Nodes (23): GWAS Analysis Pipeline, Multi-Project Orchestration with Staleness Detection, Global Configuration, GenomeSkyline, Config, String, File, ProjectState (+15 more)

### Community 7 - "Configuration Parsing"
Cohesion: 0.14
Nodes (9): Config, String, File, String, TestIntegration, File, String, TestPerProjectConfig (+1 more)

### Community 8 - "Pipeline Orchestration"
Cohesion: 0.11
Nodes (17): Config, GffParser, List, Locus, LocusOutput, Config, File, List (+9 more)

### Community 9 - "Gene Model & GFF"
Cohesion: 0.09
Nodes (16): Gene Annotation via GFF3 Streaming, Gene, Exon, String, Gene, String, Transcript (inner reference in Gene), GffParser (+8 more)

### Community 10 - "Web Frontend & Docs"
Cohesion: 0.13
Nodes (13): LYNXgwas Icon SVG, Loci Identification, Multi-Project GWAS Management, rsID Recovery Pipeline, LYNXgwas Home Page, LYNXgwas README, LYNXgwas Interactive Viewer, GlobalConfig (+5 more)

### Community 11 - "Locus Mutation Service"
Cohesion: 0.15
Nodes (13): Locus Mutation Operations (split/merge/resize/undo), Config, File, GffParser, Integer, List, Locus, LocusOutput (+5 more)

### Community 12 - "Loci Identification"
Cohesion: 0.14
Nodes (15): PLINK LD Clumping Pipeline, CandidateSnp, ClumpInterval, IdentifiedLocus, LociIdentifier, Params, LociProgress, LociProgress (+7 more)

### Community 13 - "SNP Annotation"
Cohesion: 0.13
Nodes (12): NCBI E-utilities rsID Lookup, String, Snp, Config, Integer, Map, Snp, String (+4 more)

### Community 14 - "JSON Export & LD Triangle"
Cohesion: 0.20
Nodes (13): Jb, LdTriangle, LocusRef, Runnable, Config, Exon, List, Locus (+5 more)

### Community 15 - "SNP Export Providers"
Cohesion: 0.12
Nodes (8): GwasPassthroughProvider, LdWithLeadProvider, LeadFlagsProvider, LocusIdentitySnpProvider, NearestGeneSnpProvider, RsidProvider, SnpColumnProvider, String

### Community 16 - "GWAS File Parser"
Cohesion: 0.20
Nodes (8): GwasParser, Config, List, Locus, String, Override, String, Locus

### Community 17 - "Server Endpoint Tests"
Cohesion: 0.29
Nodes (3): File, String, TestServerEndpoints

### Community 18 - "Export Column Spec"
Cohesion: 0.25
Nodes (6): ColType, Column Provider Plugin Pattern, ColumnSpec, LocusContext, SnpContext, String

### Community 19 - "Provider Value Methods"
Cohesion: 0.33
Nodes (4): ColumnSpec, LocusContext, Object, SnpContext

### Community 20 - "Multi-Project Tests"
Cohesion: 0.36
Nodes (4): File, Path, String, TestMultiProject

### Community 21 - "Locus Export Providers"
Cohesion: 0.18
Nodes (4): LocusCoreProvider, LocusLeadProvider, LocusNearestGeneProvider, LocusColumnProvider

### Community 22 - "Locus Column Provider"
Cohesion: 0.25
Nodes (6): LocusColumnProvider, ColumnSpec, List, LocusContext, Object, String

### Community 23 - "SNP Column Provider"
Cohesion: 0.25
Nodes (6): SnpColumnProvider, ColumnSpec, List, Object, SnpContext, String

### Community 24 - "RSID Detection"
Cohesion: 0.43
Nodes (4): DetectionResult, RsidDetector, List, String

### Community 25 - "Logo & Branding"
Cohesion: 0.53
Nodes (6): Data Visualization Motif, DNA Double Helix, Genomics Branding Concept, Stylized Lynx Head, LYNXgwas Logo, Network Graph Pattern

## Knowledge Gaps
- **104 isolated node(s):** `StringBuilder`, `String`, `String`, `Transcript`, `Config` (+99 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **2 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `ExcelExporter` connect `Excel Export Engine` to `Export Column Spec`?**
  _High betweenness centrality (0.057) - this node is a cross-community bridge._
- **Why does `LocalServer` connect `Local Server HTTP` to `LD & RSID Cache`, `Pipeline Orchestration`, `GWAS Pipeline Config`, `Allele Matching & DbSNP`?**
  _High betweenness centrality (0.051) - this node is a cross-community bridge._
- **Why does `GffParser` connect `Gene Model & GFF` to `Locus Mutation Service`, `GWAS Pipeline Config`?**
  _High betweenness centrality (0.036) - this node is a cross-community bridge._
- **What connects `StringBuilder`, `String`, `String` to the rest of the system?**
  _107 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `LD & RSID Cache` be split into smaller, more focused modules?**
  _Cohesion score 0.07242063492063493 - nodes in this community are weakly interconnected._
- **Should `Allele Matching & DbSNP` be split into smaller, more focused modules?**
  _Cohesion score 0.06646825396825397 - nodes in this community are weakly interconnected._
- **Should `Excel Export Engine` be split into smaller, more focused modules?**
  _Cohesion score 0.07987711213517665 - nodes in this community are weakly interconnected._