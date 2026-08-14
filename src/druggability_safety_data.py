"""Curated source data for the druggability + normal-tissue/selectivity review
phase (USP34, VEZF1, EML5, CITED2 -- frozen therapeutic shortlist, unchanged).

Data source and provenance: every field below was gathered by four
independent, source-verified research passes (one per candidate) against
UniProt, RCSB PDB, AlphaFold DB, ChEMBL, Pharos, Open Targets, GTEx Portal,
Human Protein Atlas, gnomAD, OMIM, ClinGen, MGI/IMPC, and PubMed/NCBI
eutils, each fact confirmed by directly fetching the source record (API or
page), not recalled from memory. Every accession, PMID, DOI, and numeric
value here is real and traceable to the source named in its `sources`
field; items the searches could not verify are recorded literally as
"NOT FOUND" rather than inferred or omitted. No claim of drug efficacy,
normal-tissue safety, or clinical translatability is made anywhere in this
module -- that synthesis (with its required hedges) lives in the report.

This module performs no CRISPR/RNA-seq/network/DepMap computation and
reads no upstream project data file -- it is pure curation, the
appropriate "source data" for a druggability/safety-literature phase
(the same role config/config.yaml or PREANALYSIS.md play elsewhere in this
project, and the same pattern used by
src/literature_mechanism_build_tables.py for the mechanism-review phase).
"""

from __future__ import annotations

CANDIDATES = ["USP34", "VEZF1", "EML5", "CITED2"]

# ---------------------------------------------------------------------------
# A. DRUGGABILITY
# ---------------------------------------------------------------------------
DRUGGABILITY_ROWS = [
    dict(
        candidate="USP34",
        protein_class="Cysteine-protease deubiquitinase (DUB), peptidase family C19. UniProt Q70CQ2, 3546 aa. Single annotated USP catalytic domain, residues 1894-2239; no DUSP/UBA/zinc-finger domain is UniProt-annotated for USP34 specifically (those occur in other USP-family members, not this one).",
        enzymatic_activity="Yes: cysteine-protease DUB, EC 3.4.19.12. Catalytic dyad Cys1903/His2164 (UniProt ACT_SITE); C1903S mutant = loss of function (PMID 21383061).",
        catalytic_or_binding_pocket="Catalytic cleft resolved by X-ray: apo form is conformationally INACTIVE with a misaligned catalytic histidine, realigned only upon ubiquitin engagement ('fingertips' closure mechanism) -- PMID 35588869. No dedicated pocket-druggability study (FTMap/SiteMap or similar) was found.",
        known_ligands_tools_probes_degraders="NOT FOUND. ChEMBL target search by accession returns zero USP34 entries; the one text hit for 'USP34' in ChEMBL is a mislabeled synonym on the USP35 target record (CHEMBL4630862), not real USP34 data. No chemical probe (Probes & Drugs portal unreachable/no entry), inhibitor, or degrader with a disclosed structure exists. A Progenra Inc. SBIR project (idiopathic pulmonary fibrosis indication) describes screening ~36,000 compounds against a USP34 isopeptidase assay, but discloses no hit structures or IC50s (secondary, search-level confirmation only).",
        structural_data="Experimental PDB: 7W3R (1.92 A, apo catalytic domain) and 7W3U (3.13 A, +ubiquitin-propargylamide activity-based probe), both from PMID 35588869 -- catalytic domain only, ~12% of the 3546-aa protein. No AlphaFold model exists for Q70CQ2 (AlphaFold API returns 404; no UniProt DR AlphaFoldDB cross-reference) -- an empirical absence, not confirmed as a length-limit artifact.",
        targetable_ppi_interfaces="AXIN1/AXIN2 binding (UniProt SUBUNIT comment, PMID 21383061); no AXIN-USP34 co-structure exists. The ubiquitin-binding 'fingertips' interface within the catalytic domain (PMID 35588869) is a documented substrate interface, not a separate protein partner interface.",
        direct_small_molecule_realistic="The only one of the four candidates with a real, solved catalytic-domain structure and a defined two-state (apo-inactive / substrate-bound-active) mechanism -- a plausible structure-based drug-discovery starting point in principle. No USP34-specific chemical matter exists yet, and DUB catalytic clefts are generally shallow and historically hard to drug selectively.",
        plausible_alternative_modalities="Allosteric/conformational-trap inhibitor exploiting the apo-inactive state (by loose analogy to USP14/IU1, not demonstrated for USP34); antisense/siRNA; CRISPR/epigenetic perturbation.",
        clinical_preclinical_precedent="None for USP34 itself (only the undisclosed-structure SBIR screen). Related DUBs, class-level context only: USP30 inhibitors (Mission Therapeutics/Vincere Biosciences) are in IND-enabling/early-clinical development for Parkinson's disease via mitophagy; several USP7 tool compounds (P5091, P22077, FT671, FT827, GNE-6640, GNE-6776) exist preclinically with none in clinical trials per a 2022 review.",
        patent_activity="NOT FOUND as USP34-specific. One candidate patent (US9752151B2) surfaced in search snippets but its full text contains zero occurrences of 'USP34' on direct verification and is excluded as not USP34-specific.",
        druggability_classification="POTENTIALLY_DRUGGABLE",
        sources="UniProt Q70CQ2; PMID 21383061; PMID 35588869 (PDB 7W3R, 7W3U); ChEMBL REST API; AlphaFold API (404, confirmed absent); SBIR.gov (secondary/search-level)",
    ),
    dict(
        candidate="VEZF1",
        protein_class="C2H2 zinc-finger transcription factor. UniProt Q14119, six tandem zinc fingers (ZnF1-6, residues 74-308) plus a C-terminal repeat region (394-485). No catalytic domain of any kind.",
        enzymatic_activity="No. UniProt lists no catalytic-activity comment and no EC number; function is DNA-binding transcriptional regulation (e.g. binds the IL-3 promoter CT/GC-rich region).",
        catalytic_or_binding_pocket="No small-molecule binding pocket is UniProt-annotated (only the six structural, zinc-coordinating fingers). One published approach targets the DNA-binding interface itself via a homology model (built on the Zif268/EGR1 DNA-bound zinc-finger structure), not a classical enzyme pocket (PMID 29970794).",
        known_ligands_tools_probes_degraders="One VEZF1-specific tool-compound series exists: He et al. 2018 (PMID 29970794) computationally docked the NCI Diversity Library against a VEZF1 homology model and identified three hits that block VEZF1-DNA binding by EMSA -- T4 (IC50=20 uM, abolished endothelial tube formation without affecting viability at <=IC50), T6 (IC50=100 uM), NSC1012 (IC50=500 uM). These are unoptimized, low-potency hit-level compounds, not chemical probes or leads. ChEMBL: zero target entries. No degrader/PROTAC found.",
        structural_data="No experimental PDB structure exists for human VEZF1 (RCSB full-text search returns zero entries, confirmed twice). AlphaFold model AF-Q14119-F1 exists but is LOW CONFIDENCE (global pLDDT 59.4; 47.8% of the model 'very low confidence'), consistent with a largely intrinsically-disordered protein outside the six zinc fingers.",
        targetable_ppi_interfaces="VEZF1 directly represses CITED2 transcription (ChIP-verified in mouse ESCs, PMID 29794136) -- a documented regulatory link between two of this project's own four candidates. Target-gene/co-regulatory interactions also documented for VEGFR2/NRP1 and TIMP3/MMP2 (human primary endothelial cells, PMID 24280686) and ETV2 (direct co-IP-verified protein interaction, mouse ESC/EB model, PMID 36923254). VEZF1 protein is a substrate of STUB1(CHIP)-mediated ubiquitination/degradation in an HCC model (PMID 36241701, lower-confidence sourcing -- primary text was paywalled). None of these constitute a structurally characterized, druggable PPI interface.",
        direct_small_molecule_realistic="Zinc-finger transcription factors have no catalytic active site; VEZF1 additionally has no experimental structure of any kind and only a low-confidence AlphaFold model. The only tool compounds are weak (20-500 uM), homology-model-derived DNA-binding-interface blockers with no medicinal-chemistry optimization reported. Pharos classifies VEZF1 as target-development-level 'Tbio' (known biology, no qualifying bioactive ligand).",
        plausible_alternative_modalities="Targeted degradation/PROTAC or molecular glue of the transcription factor (no VEZF1-specific example exists, but STUB1-mediated natural degradation demonstrates the protein is degradation-controllable); antisense/siRNA (miR-191 is reported to directly target VEZF1 mRNA in an ischemic-stroke angiogenesis model, PMID 31064890, demonstrating the transcript is a druggable RNA node by natural-regulation precedent); CRISPR/epigenetic perturbation; further medicinal-chemistry optimization of the He et al. 2018 DNA-binding-interface blockers.",
        clinical_preclinical_precedent="None for VEZF1 itself. Class-level zinc-finger-TF drugging precedent was not independently source-verified in this research pass and is not reported as a specific claim.",
        patent_activity="NOT FOUND. Search returned no VEZF1-specific patents.",
        druggability_classification="CURRENTLY_POORLY_DRUGGABLE",
        sources="UniProt Q14119; PMID 29970794 (tool compounds); RCSB Search API (0 hits); AlphaFold API (AF-Q14119-F1, pLDDT 59.4); PMID 29794136; PMID 24280686; PMID 36923254; PMID 36241701 (secondary); PMID 31064890; ChEMBL REST API",
    ),
    dict(
        candidate="EML5",
        protein_class="Non-catalytic, microtubule-associated WD40/HELP-domain scaffold protein (EML family). UniProt Q05BV3, 1969 aa. 29 annotated WD-repeat features plus HELP-domain and beta-propeller (EML-domain) InterPro/Pfam cross-references; no coiled-coil region is annotated in this UniProt entry.",
        enzymatic_activity="No. UniProt function comment (evidence: inferred by homology, not direct human data): 'May modify the assembly dynamics of microtubules.' No catalytic activity, no EC number.",
        catalytic_or_binding_pocket="NOT FOUND. No druggability assessment of any kind located. Pharos classifies EML5 as target-development-level 'Tdark' -- its lowest tier, denoting minimal functional annotation and no qualifying bioactive ligand.",
        known_ligands_tools_probes_degraders="NOT FOUND. ChEMBL target search returns zero entries (confirmed directly). No inhibitor, probe, degrader, or PROTAC found in PubMed/patent search.",
        structural_data="No experimental PDB structure (RCSB search returns zero entries). AlphaFold model AF-Q05BV3-F1 exists with a comparatively high mean pLDDT of 87.0 (60% of residues >90 pLDDT), consistent with a well-modeled WD40 beta-propeller core.",
        targetable_ppi_interfaces="No dedicated microtubule-binding-interface structural data found. UniProt/STRING/BioGRID list only low-confidence, largely genomic-neighbor-driven computational associations (e.g. PTPN21, TTC8, ZC3H14, SPATA7 -- all physically adjacent genes at the 14q31.3 locus per PMID 31406157) dominated by text-mining rather than experimental evidence. No validated, targetable interactor identified.",
        direct_small_molecule_realistic="No enzymatic active site, no known ligand-binding pocket, zero ChEMBL bioactivity, Pharos 'Tdark' status, and no experimental structure. The single most sparsely characterized of the four candidates by every druggability metric checked.",
        plausible_alternative_modalities="No validated modality-specific rationale identified for EML5 itself; in principle, PPI disruption of a validated interactor (none currently known), targeted degradation, or antisense/siRNA/CRISPR knockdown would be the generic non-small-molecule options for an uncharacterized scaffold protein.",
        clinical_preclinical_precedent="None for EML5. Family-level context only (NOT EML5-specific): the EML5 paralog EML4 is the fusion partner in EML4-ALK, the oncogenic fusion kinase targeted by ALK inhibitors (e.g. crizotinib) in non-small-cell lung cancer (PMID 17625570) -- the drugged moiety there is ALK's kinase domain, not any EML4/EML5 region.",
        patent_activity="NOT FOUND.",
        druggability_classification="CURRENTLY_POORLY_DRUGGABLE",
        sources="UniProt Q05BV3; Pharos (Tdark); ChEMBL REST API (0 hits); RCSB Search API (0 hits); AlphaFold API (AF-Q05BV3-F1, pLDDT 87.0); STRING API; PMID 31406157; PMID 17625570 (family-level context)",
    ),
    dict(
        candidate="CITED2",
        protein_class="Intrinsically disordered transcriptional co-regulator with NO DNA-binding domain. UniProt Q99967, 270 aa. Acts primarily -- though not exclusively -- through the C-terminal 'CR2' transactivation domain (~residues 216-269), which is disordered when free and folds only upon binding its partners; UniProt SUBUNIT additionally lists CITED2 interactions with SMAD2/3, TFAP2A/B/C, WT1 (by similarity), LHX2, and PPARA, none of which are independently co-structured (CR2-CBP/p300 is the only co-structured interface, see below).",
        enzymatic_activity="No. UniProt GO terms are limited to transcription coactivator/corepressor activity, chromatin binding, and histone-acetyltransferase binding -- no catalytic GO term, no EC number. Confirmed non-enzymatic PPI-based co-regulator.",
        catalytic_or_binding_pocket="No classical pocket (there is no catalytic domain), but the CITED2-CBP/p300 interface IS structurally characterized (unlike VEZF1 or EML5): NMR structures PDB 1P4Q (CITED2 216-259 + p300 CH1 domain, PMID 12778114) and 1R8U (CITED2 220-269 + CBP TAZ1 domain, PMID 14594809) show a coupled folding-and-binding interaction over an extended, non-pocket-like groove, the same TAZ1/CH1 surface competitively used by HIF-1alpha. Two later engineered-fusion co-structures (7LVS, 7QGS) further dissect this same interface's energetics.",
        known_ligands_tools_probes_degraders="NOT FOUND. ChEMBL: zero target entries. Open Targets 'chemicalProbes' field is an empty list; 'drugAndClinicalCandidates' count = 0. No compound intentionally targeting CITED2 or the CITED2-p300/CBP interaction was found. One indirect finding: p300/CBP-directed PROTACs (e.g. CBPD-409) cause secondary loss of CITED2 protein as a downstream consequence of degrading its binding partners -- a pharmacodynamic side effect of a different drug, not a CITED2-targeted therapeutic.",
        structural_data="PDB 1P4Q, 1R8U (NMR, native complex fragments, residues ~216-269 only), 7LVS, 7QGS (X-ray, engineered fusion/hybrid peptides). No full-length (1-270) structure exists, expected for an IDP where >60% of the sequence has no fixed structure outside the bound state. AlphaFold model AF-Q99967-F1 is LOW CONFIDENCE (mean pLDDT 51.5; 57% of residues 'very low confidence'), consistent with predominant intrinsic disorder.",
        targetable_ppi_interfaces="The CITED2-CBP/p300 (TAZ1/CH1) interface is the best-characterized target-relevant surface of any candidate in this review -- real native-complex NMR structures exist (1P4Q, 1R8U) showing an extended, non-globular, coupled-folding interaction (not a compact hydrophobic pocket), competitively shared with HIF-1alpha; the source papers describe this as a high-affinity, extensive interface, not a small/shallow one. UniProt SUBUNIT also lists CITED2 interactions with SMAD2/3, TFAP2A/B/C, WT1 (by similarity), LHX2, PPARA -- none independently co-structured.",
        direct_small_molecule_realistic="No catalytic active site exists by definition (no enzyme). Open Targets' formal tractability assessment returns FALSE across every small-molecule bucket checked (Approved Drug, Advanced Clinical, Phase 1, Structure-with-Ligand, High/Med-Quality Pocket, Druggable Family -- all false), consistent with the extended, disordered, IDP-typical binding mode described above being generally regarded as difficult for classical pocket-based small-molecule inhibitors.",
        plausible_alternative_modalities="PPI-interface disruptor (stapled or macrocyclic peptide mimicking the CITED2 activation-domain helix/extended element, given the real structural template at 1P4Q/1R8U); molecular glue; targeted degradation/PROTAC (would require a novel, currently unidentified small-molecule 'hook' on CITED2 -- none exists); antisense/siRNA/ASO knockdown; CRISPR/epigenetic perturbation of the CITED2 locus.",
        clinical_preclinical_precedent="None for CITED2 itself (0 ChEMBL entries, 0 Open Targets clinical candidates, 0 chemical probes). The shared binding hub, p300/CBP, IS a validated clinical target in its own right (bromodomain inhibitors e.g. CCS1477/inobrodib in clinical trials; PROTACs e.g. dCBP1, QC-182) -- but that is a distinct target, not CITED2.",
        patent_activity="NOT FOUND as a CITED2-directed composition-of-matter/therapeutic patent; only incidental mentions of CITED2 as a gene-expression biomarker in unrelated diagnostic patents.",
        druggability_classification="INDIRECT_OR_MODALITY_DEPENDENT",
        sources="UniProt Q99967; PDB 1P4Q (PMID 12778114); PDB 1R8U (PMID 14594809); PDB 7LVS, 7QGS; ChEMBL REST API; Open Targets Platform GraphQL API (tractability, chemicalProbes fields); AlphaFold API (AF-Q99967-F1, pLDDT 51.5)",
    ),
]

# ---------------------------------------------------------------------------
# B. NORMAL-TISSUE EXPRESSION (GTEx v8 + Human Protein Atlas)
# ---------------------------------------------------------------------------
NORMAL_TISSUE_ROWS = [
    dict(
        candidate="USP34",
        gtex_top_tissues="Cells-EBV-transformed lymphocytes 38.7 TPM; Cells-Cultured fibroblasts 35.7; Artery-Tibial 30.1; Ovary 28.9; Uterus 28.9; Brain-Cerebellum 27.6; Testis 26.1 (GTEx v8, all 54 tissues expressed, no near-zero tissue)",
        gtex_breast_tpm=18.61,
        gtex_blood_tpm=4.28,
        gtex_brain_range_tpm="6.6-27.6 (all 13 GTEx brain subregions)",
        gtex_liver_tpm=4.06,
        gtex_heart_tpm=5.44,
        gtex_kidney_tpm=6.71,
        gtex_gi_tract_tpm="12.5-24.0 (stomach to colon)",
        gtex_reproductive_tpm="25.1-28.9 (testis, ovary, uterus)",
        hpa_tissue_specificity_category="Low tissue specificity (Tau 0.14); Low cell type specificity (Tau 0.17); Low immune cell specificity",
        hpa_notes="HPA consensus nTPM broadly expressed 34-72 across most tissues; skeletal muscle notably high (54.1 nTPM, 3rd-highest of 50 HPA tissues); liver lowest GTEx tissue (4.06 TPM) but moderate in HPA (39.9 nTPM) -- a real cross-dataset normalization discrepancy, reported as-is.",
        expression_breadth_summary="Broadly and fairly uniformly expressed across essentially all sampled normal tissues (whole blood and liver are the comparative low points; all values remain well above noise). Not tissue-restricted.",
        sources="GTEx Portal API v2 (gtex_v8, ENSG00000115464); HPA proteinatlas.org/ENSG00000115464-USP34 (tissue + rna_tissue_consensus.tsv)",
    ),
    dict(
        candidate="VEZF1",
        gtex_top_tissues="Uterus 58.5 TPM; Cervix-Endocervix 54.5; Colon-Sigmoid 53.5; Fallopian tube 49.5; Cervix-Ectocervix 47.7; Esophagus-Muscularis 48.5 (GTEx v8, 50 tissues, range 6.7-58.5, <9-fold spread)",
        gtex_breast_tpm=30.5,
        gtex_blood_tpm=9.8,
        gtex_brain_range_tpm="7.9-23.9 (13 GTEx brain subregions; lowest major-organ group along with liver)",
        gtex_liver_tpm=6.7,
        gtex_heart_tpm=14.1,
        gtex_kidney_tpm=14.1,
        gtex_gi_tract_tpm="20.7-53.5 (stomach to colon)",
        gtex_reproductive_tpm="25.1-58.5 (testis to uterus, highest tissue in dataset)",
        hpa_tissue_specificity_category="Low tissue specificity; 'ubiquitous nuclear expression' (HPA consensus description)",
        hpa_notes="Highest single HPA value: brain corpus callosum (white matter), 95.6 nTPM. Cell-type enrichment noted in adrenal cortex cells, skin/splenic endothelial cells, testicular spermatogonia -- consistent with, but not exclusive to, VEZF1's namesake endothelial biology. Low immune cell specificity across all major leukocyte subsets.",
        expression_breadth_summary="Broad, low-tissue-specificity expression across nearly all sampled tissues; liver is the single lowest GTEx tissue. GTEx/HPA artery values reflect bulk vessel tissue, not purified endothelium.",
        sources="GTEx Portal API v2 (gtex_v8, ENSG00000136451); HPA proteinatlas.org/ENSG00000136451-VEZF1 (tissue + immune cell)",
    ),
    dict(
        candidate="EML5",
        gtex_top_tissues="Ovary 9.62 TPM (highest of 54 tissues); Brain-Cerebellum 5.55; Brain-Cerebellar Hemisphere 3.93 (other brain subregions much lower, 0.22-0.60); Testis 2.46; Pituitary 1.89",
        gtex_breast_tpm=0.273,
        gtex_blood_tpm=0.021,
        gtex_brain_range_tpm="0.22-5.55 (cerebellum-predominant; other regions near-zero)",
        gtex_liver_tpm=0.033,
        gtex_heart_tpm=0.0195,
        gtex_kidney_tpm=0.095,
        gtex_gi_tract_tpm="0.28-0.87 (stomach to esophagus)",
        gtex_reproductive_tpm="1.16-9.62 (prostate to ovary, highest tissue in dataset)",
        hpa_tissue_specificity_category="Tissue enhanced (ovary, retina) -- per UniProt/HPA cross-reference",
        hpa_notes="HPA consensus: ovary 10.4 nTPM (highest), retina 8.4; most other tissues <=0.7 nTPM. Notable mRNA/protein discordance: ovary immunohistochemistry flagged 'ND' (not detected at protein level) despite highest mRNA nTPM. Immune single-cell data: highest in plasmablasts (4.3 nTPM); neutrophils/monocytes/eosinophils = not detected.",
        expression_breadth_summary="Genuinely tissue-restricted (ovary/retina-enriched, cerebellum-predominant in brain), with very low or undetectable expression across most other normal tissues -- the most tissue-restricted expression pattern of the four candidates.",
        sources="GTEx Portal API v2 (gtex_v8, ENSG00000165521); HPA proteinatlas.org/ENSG00000165521-EML5 (tissue + single cell)",
    ),
    dict(
        candidate="CITED2",
        gtex_top_tissues="Cells-Cultured fibroblasts 492.1 TPM; Ovary 360.7; Fallopian tube 227.6; Nerve-Tibial 206.3; Thyroid 206.2 (GTEx v8, broad expression, no near-zero tissue)",
        gtex_breast_tpm=133.7,
        gtex_blood_tpm=26.3,
        gtex_brain_range_tpm="4.2-73.7 (13 GTEx brain subregions, widest regional spread of the four candidates)",
        gtex_liver_tpm=54.8,
        gtex_heart_tpm=30.0,
        gtex_kidney_tpm=43.1,
        gtex_gi_tract_tpm="39.6-139.6 (stomach to esophagus)",
        gtex_reproductive_tpm="31.2-360.7 (testis to ovary, highest tissue in dataset)",
        hpa_tissue_specificity_category="Low tissue specificity (Tau 0.26 tissue-level, 0.46 single-cell level); RNA cluster 'Ovary - Mixed function'",
        hpa_notes="Single-cell 'Cell type enhanced' in extravillous trophoblasts, breast lactating cells, Hofbauer cells, migrating cytotrophoblasts. Highest brain nTPM: cerebellar cortex 93.9. Not detected in blood plasma by immunoassay or MS proteomics (not a circulating protein at detectable levels).",
        expression_breadth_summary="Broadly and, in absolute terms, highly expressed across nearly all normal tissues sampled (among the higher absolute-TPM candidates of the four); reproductive tissues (ovary, fallopian tube) and cultured fibroblasts are highest; heart and whole blood are comparatively lower but clearly non-zero.",
        sources="GTEx Portal API v2 (gtex_v8, ENSG00000164442); HPA proteinatlas.org/ENSG00000164442-CITED2 (tissue + single cell)",
    ),
]

# ---------------------------------------------------------------------------
# C. GENETIC CONSTRAINT (gnomAD, OMIM, ClinGen, MGI/IMPC)
# ---------------------------------------------------------------------------
GENETIC_CONSTRAINT_ROWS = [
    dict(
        candidate="USP34",
        gnomad_version="v4.0/v4.1 constraint schema (cross-checked against gnomAD's own 2024-03 v4.0 gene-constraint post; live API query 2026-08-14)",
        loeuf=0.152,
        loeuf_90pct_ci="0.125-0.186",
        pli=1.0,
        oe_lof_observed_expected="70 observed / 460.3 expected",
        missense_z=1.65,
        omim_entry="*615295 (gene entry; phenotype-linkage status to the 2026 haploinsufficiency paper not directly confirmable, OMIM full-text fetch 403-blocked)",
        clingen_dosage_classification="No formal ClinGen 0-3 haploinsufficiency/triplosensitivity score on record. ClinGen gene-disease validity: Congenital Heart Disease classification 'Limited' (AD, evaluated 2024-06-10). DECIPHER haploinsufficiency index (older metric) = 12.35.",
        omim_clingen_discrepancy="No direct contradiction identified, but ClinGen's formal curation lags the most recent (2026) primary human-genetics finding (see mouse_ko/human_phenotype column).",
        mouse_ko_phenotype_summary="Conditional (MSC/pre-osteoblast-Cre) Usp34 knockout: reduced osteogenic differentiation and bone formation, impaired BMP2-driven bone regeneration (PMID 30181118, EMBO J 2018). No constitutive/germline knockout viability data located.",
        human_phenotype_summary="Wigoda et al. 2026 (PMID 42315110, Clin Genet): 6 individuals with heterozygous LoF USP34 variants (5 de novo) -- global developmental delay, craniofacial dysmorphism, speech impairment, variable ASD, AND distal limb anomalies, overlapping 2p15p16.1 microdeletion syndrome. Human, germline, developmental evidence -- not an adult/pharmacological-inhibition study.",
        impc_status="0/24 physiological systems tested; 0 significant phenotypes reported (not yet systematically phenotyped)",
        normal_cell_essentiality_data="NOT FOUND -- no independently verified normal (non-cancer) human-cell CRISPR essentiality data point located for USP34.",
        sources="gnomAD GraphQL API (gnomad.broadinstitute.org/api, ENSG00000115464); OMIM *615295 (403-blocked, unconfirmed phenotype link); ClinGen HGNC:20066; PMID 30181118; PMID 42315110; MGI:109473; IMPC mousephenotype.org",
    ),
    dict(
        candidate="VEZF1",
        gnomad_version="v4.1.1 (gnomad.broadinstitute.org/gene page self-labeled; API cross-checked, live query 2026-08-14)",
        loeuf=0.24,
        loeuf_90pct_ci="0.06-0.24",
        pli=1.0,
        oe_lof_observed_expected="5 observed / 44.1 expected",
        missense_z=4.78,
        omim_entry="*606747 (gene); disease link Cardiomyopathy, Dilated, 1OO (CMD1OO, MIM #620247, AD) -- single 4-generation family, exome-identified nonsense variant p.(Lys164*), functional MYH7/EDN1 transactivation-failure assay",
        clingen_dosage_classification="0 dosage-sensitivity classifications on record (no haploinsufficiency/triplosensitivity score).",
        omim_clingen_discrepancy="YES -- ClinGen's Dilated Cardiomyopathy Gene Curation Expert Panel formally classifies VEZF1-CMD1OO as 'No Known Disease Relationship' (evaluated 2026-03-04), directly superseding/contradicting the single-family OMIM disease listing. Reported as an unresolved discrepancy, not silently resolved.",
        mouse_ko_phenotype_summary="Homozygous Vezf1-null: embryonic lethal at midgestation (angiogenic remodeling defects, hemorrhaging, defective EC adhesion/tight junctions) -- PMID 15882861. Heterozygous: ~20%-penetrant lymphatic hypervascularization with hemorrhaging/edema (cystic-hygroma-like). MGI phenotype categories: cardiovascular, embryonic, growth/size, immune, mortality -- no skeletal category flagged.",
        human_phenotype_summary="Only the single-family, ClinGen-disputed dilated-cardiomyopathy report (see omim_entry). No other confirmed human Mendelian phenotype located.",
        impc_status="0 significant phenotypes reported; viability status not populated (not systematically phenotyped)",
        normal_cell_essentiality_data="NOT FOUND -- not searched in this pass beyond the explicit DepMap exclusion.",
        sources="gnomAD browser + API (ENSG00000136451, v4.1.1); OMIM *606747 (search-indexed, primary full text 403-blocked); ClinGen HGNC:12949 (direct fetch); PMID 15882861; MGI:1313291; IMPC mousephenotype.org",
    ),
    dict(
        candidate="EML5",
        gnomad_version="v4 series (GraphQL API query 2026-08-14; ClinGen cross-check consistent)",
        loeuf=0.558,
        loeuf_90pct_ci="0.414-0.558",
        pli=9.18e-11,
        oe_lof_observed_expected="120 observed / 250.0 expected",
        missense_z=3.73,
        omim_entry="*618119 (gene-only entry, asterisk prefix = no associated phenotype MIM number, i.e. no established human disease phenotype). Unverified GeneCards-attributed disease names ('Leber Plus Disease', 'Mucopolysaccharidosis Type IVb') were flagged as likely spurious and are NOT reported as established facts (GeneCards page itself 403-blocked, no corroborating literature).",
        clingen_dosage_classification="Haploinsufficiency and Triplosensitivity both 'Awaiting Review' (not yet curated, distinct from 'no evidence found'). DECIPHER %HI = 30.61%.",
        omim_clingen_discrepancy="None identified (neither source asserts a phenotype).",
        mouse_ko_phenotype_summary="MGI:2442513 explicitly: '0 phenotypes from 0 alleles in 0 genetic backgrounds' despite 6 alleles existing.",
        human_phenotype_summary="No human Mendelian phenotype established (gene-only OMIM entry; ClinGen dosage review pending). One non-human Mendelian-type association exists: OMIA:002904-9913, EML5-related male subfertility in cattle (PMID 35478957) -- explicitly bovine, not human.",
        impc_status="0/24 physiological systems tested; 0 significant phenotypes; 0 associated diseases (not yet systematically phenotyped)",
        normal_cell_essentiality_data="NOT FOUND as a normal-tissue-specific data point; UniProt cross-references BioGRID-ORCS ('9 hits in 1151 CRISPR screens') but this aggregates predominantly cancer-cell-line screens and is excluded from normal-tissue inference per project scope.",
        sources="gnomAD GraphQL API (ENSG00000165521); OMIM *618119 (title only, 403-blocked); ClinGen HGNC:18197 (direct fetch); MGI:2442513 (direct fetch); IMPC mousephenotype.org; PMID 35478957",
    ),
    dict(
        candidate="CITED2",
        gnomad_version="v4.1.1 (807,162 samples, GRCh38, MANE Select transcript ENST00000367651.4 -- confirmed via both live browser page and GraphQL API, values match)",
        loeuf=0.17,
        loeuf_90pct_ci="0.07-0.54",
        pli=0.94,
        oe_lof_observed_expected="2 observed / 11.8 expected",
        missense_z=-1.01,
        omim_entry="#614431 Ventricular Septal Defect 2 (VSD2) and #614433 Atrial Septal Defect 8 (ASD8), both linked to specific in-frame CITED2 deletions/insertions retaining 50-75% of wild-type transactivation/repression activity (PMID 16287139). An initially-surfaced association with OMIM #614980 (CHTD2) was directly verified to be FALSE -- CHTD2 is caused by TAB2, not CITED2 -- and is excluded.",
        clingen_dosage_classification="Haploinsufficiency score 0, 'No Evidence for Haploinsufficiency'; Triplosensitivity score 0, 'No Evidence for Triplosensitivity' (evaluated 2023-11-29). ClinGen's stated rationale: reported variants only marginally affect repressive activity, some are inherited from unaffected parents (incomplete penetrance), and were found in both patients and controls.",
        omim_clingen_discrepancy="YES -- OMIM lists two named CITED2-linked congenital heart defect phenotypes (VSD2, ASD8) based on rare hypomorphic variants, while ClinGen's current, systematic, expert-panel review finds insufficient evidence for haploinsufficiency. Reported as an unresolved discrepancy, not silently resolved.",
        mouse_ko_phenotype_summary="Constitutive Cited2-null: embryonic/perinatal lethal with cardiac malformations (VSD, ASD, double-outlet right ventricle, persistent truncus arteriosus), neural tube defects (exencephaly), adrenal agenesis, cranial neural crest and left-right patterning defects (PMID 11694877, 12149478, 15750185). None of these three papers report a distinct skeletal phenotype (absence of report, not a documented negative finding).",
        human_phenotype_summary="See omim_entry; disputed by ClinGen (see omim_clingen_discrepancy).",
        impc_status="Dozens of MGI/IMPC-curated Mammalian Phenotype terms (via Open Targets aggregation), overwhelmingly cardiac/neural-tube/laterality/adrenal/lethality-related; no explicit skeletal/bone Mammalian Phenotype term appears in the returned list (reflects what has been curated, not proof of formal exclusion).",
        normal_cell_essentiality_data="NOT FOUND within available tools; Open Targets 'isEssential' = false (a DepMap-cancer-cell-line-derived field, excluded from normal-tissue inference per project scope).",
        sources="gnomAD browser + API (ENSG00000164442, v4.1.1); OMIM #614431, #614433 (UniProt DISEASE CC cross-refs); PMID 16287139; ClinGen CCID:006871 (direct fetch); PMID 11694877; PMID 12149478; PMID 15750185; Open Targets Platform GraphQL API",
    ),
]

# ---------------------------------------------------------------------------
# D. BONE / MUSCULOSKELETAL SYSTEM
# ---------------------------------------------------------------------------
BONE_MUSCULOSKELETAL_ROWS = [
    dict(
        candidate="USP34",
        bone_marrow_expression_hpa_ntpm=35.5,
        skeletal_muscle_expression_hpa_ntpm=54.1,
        cell_type_resolved_data_available="No osteoblast/osteocyte/osteoclast/chondrocyte entry exists in HPA's single-cell atlas for any candidate in this review (154 cell types profiled, none bone-lineage). Closest musculoskeletal-adjacent single-cell entries for USP34: 'myonuclei' (skeletal-muscle nuclei) at 766.5 nCPM (2nd-highest of 154 HPA cell types) and hematopoietic stem cells at 395.0 nCPM.",
        published_bone_role_summary="YES, direct and causal: Guo et al. 2018 (PMID 30181118, EMBO J) -- USP34 stabilizes Smad1 and RUNX2 (opposing SMURF1-mediated degradation) to promote osteogenic differentiation; siRNA/shRNA depletion in human MSCs in vitro inhibits osteogenic differentiation; conditional (MSC- or pre-osteoblast-Cre) knockout mice show low bone mass and impaired BMP2-driven bone regeneration in vivo, rescued by Smurf1 co-depletion.",
        evidence_species="Human (MSCs, in vitro) + Mouse (conditional knockout, in vivo)",
        evidence_type="DOCUMENTED_CAUSAL_POSTNATAL (mouse conditional knockout using MSC/pre-osteoblast-Cre, active during skeletal development/differentiation -- NOT an adult-onset deletion in already-mature bone, unlike CITED2's Mx1-Cre studies below) plus DOCUMENTED_DEVELOPMENTAL_HUMAN (2026 human germline haploinsufficiency causing distal limb anomalies as part of a neurodevelopmental syndrome, PMID 42315110 -- see genetic-constraint table)",
        bone_concern_category="DOCUMENTED_CAUSAL_POSTNATAL",
        sources="HPA rna_tissue_consensus.tsv, rna_single_cell_type.tsv; PMID 30181118; PMID 42315110",
    ),
    dict(
        candidate="VEZF1",
        bone_marrow_expression_hpa_ntpm=53.5,
        skeletal_muscle_expression_hpa_ntpm=33.2,
        cell_type_resolved_data_available="No bone-lineage single-cell data located. Bone-marrow tissue-level IHC is discordant between antibodies (HPA027520 'not detected' vs HPA048315 'medium' in the same hematopoietic-cell population); MS proteomics 'not detected' in all 3 replicates -- a genuine cross-method discrepancy, not a clean positive protein-level finding.",
        published_bone_role_summary="NOT FOUND directly. Closest indirectly-relevant paper: Das et al. 2023 (PMID 36923254) shows Vezf1-KO mouse embryoid bodies have reduced hematoendothelial marker genes (Flt1, Gata2, Flk1) -- developmental, embryonic hematoendothelial lineage-specification evidence, upstream of (not the same as) adult bone-marrow niche biology. Must not be over-read as direct bone-marrow-niche or osteogenesis evidence.",
        evidence_species="Mouse (embryonic hematoendothelial lineage only; no direct bone/marrow-niche study in any species)",
        evidence_type="INFERRED_ONLY (no direct bone/skeletal paper exists; only upstream developmental-lineage inference plus generic vascularization-dependent reasoning about endochondral ossification, which this review explicitly declines to extrapolate as evidence)",
        bone_concern_category="INFERRED_ONLY",
        sources="HPA proteinatlas.org/ENSG00000136451-VEZF1/tissue/bone+marrow; PMID 36923254; PMID 15882861 (no skeletal phenotype reported); MGI:1313291 (no skeletal phenotype category flagged)",
    ),
    dict(
        candidate="EML5",
        bone_marrow_expression_hpa_ntpm=0.3,
        skeletal_muscle_expression_hpa_ntpm=0.0,
        cell_type_resolved_data_available="No bone-lineage single-cell data located; no cell-type breakdown available for the low bone-marrow bulk-tissue signal.",
        published_bone_role_summary="NOT FOUND. PubMed esearch 'EML5 AND (bone OR osteoblast OR osteoclast OR skeletal)' returns 0 results; EML5's entire PubMed literature is only 17 papers total, none concerning bone biology.",
        evidence_species="None available",
        evidence_type="NONE_IDENTIFIED (MGI: 0 phenotypes from 0 alleles despite 6 alleles existing; IMPC: 0/24 systems tested -- absence of testing, not a negative finding, but no positive bone evidence of any kind exists to weigh)",
        bone_concern_category="NONE_IDENTIFIED",
        sources="HPA rna_tissue_consensus.tsv; NCBI eutils esearch (EML5 AND bone/osteoblast/osteoclast/skeletal, 0 hits); MGI:2442513; IMPC mousephenotype.org",
    ),
    dict(
        candidate="CITED2",
        bone_marrow_expression_hpa_ntpm=None,
        skeletal_muscle_expression_hpa_ntpm=None,
        cell_type_resolved_data_available="No dedicated osteoblast single-cell entry in HPA (HPA's single-cell 'cell type enhanced' calls for CITED2 are trophoblast/lactating-breast/Hofbauer-cell/eosinophil categories, not bone-lineage). No single cross-cell-type ranking of osteoblasts vs osteoclasts vs HSCs vs MSCs by CITED2 expression was located.",
        published_bone_role_summary="YES, the most extensive and mechanistically direct bone/marrow literature of the four candidates, spanning four distinct axes: (1) HSC/bone-marrow maintenance -- adult conditional (Mx1-Cre) Cited2 knockout causes HSC loss, loss of quiescence, and multilineage bone-marrow failure in mice (PMID 19951693, 22308296, 34715054), with human confirmatory data in primary CD34+ cells/AML xenografts (PMID 25184385); (2) osteoclast differentiation -- Cited2 is described as 'the molecular switch triggering terminal differentiation of osteoclasts', with in vivo lineage-specific deletion in osteoclast precursors causing failure to commit to the osteoclast fate (PMID 33288951, Nat Metab 2020); (3) fracture healing -- CITED2 identified as a negative regulator of fracture healing via MMP-2/-3/-9/-13 suppression in a rat mandibular-osteotomy model (PMID 19607804); (4) chondrocyte/cartilage mechanotransduction -- CITED2 mediates MMP-1/MMP-13 repression under physiological joint loading in human chondrocytes and rat in vivo joint-loading models, chondroprotective when moderately loaded (PMID 12960175, 20826544). Also implicated in adipogenic-vs-osteogenic lineage commitment of human MSCs (siRNA knockdown redirects hMSCs toward adipogenic fate even under osteogenic-inducing conditions, PMID 33691475).",
        evidence_species="Mouse (adult conditional knockout, postnatal) + Human (primary CD34+ cells, AML xenografts, primary chondrocytes) + Rat (fracture healing, joint loading)",
        evidence_type="DOCUMENTED_ADULT_CAUSAL: the HSC-maintenance studies (PMID 19951693, 22308296) used Mx1-Cre, genuinely induced in already-mature adult mice (poly(I:C)-inducible), not a developmental promoter -- the only candidate in this review with a conditional deletion performed in already-formed adult tissue. Osteoclastogenesis (PMID 33288951), fracture healing (rat, PMID 19607804), and cartilage mechanotransduction (rat/human, PMID 12960175, 20826544) are additional, independent postnatal/adult-tissue axes. Separately, constitutive knockout is embryonic lethal (cardiac/neural-tube, PMID 11694877/12149478/15750185) without a distinctly reported skeletal phenotype in those specific papers.",
        bone_concern_category="DOCUMENTED_ADULT_CAUSAL",
        sources="PMID 19951693; PMID 22308296; PMID 34715054; PMID 25184385; PMID 33288951; PMID 19607804; PMID 12960175; PMID 20826544; PMID 33691475",
    ),
]

# ---------------------------------------------------------------------------
# E. VERIFIED REFERENCES -- every PMID cited by letter (candidate) throughout
# this module. PMID/journal/year identity was confirmed by the research
# agents' direct source fetches (NCBI eutils/EuropePMC/PMC record) for
# every row. Titles were independently re-verified against Codex's own
# PubMed check for a sample of load-bearing rows during the phase-review
# pass (2026-08-14); rows not in that sample carry an abbreviated/working
# title (the paper's real title paraphrased from the fetched abstract,
# not retyped character-for-character from PubMed) -- this is flagged
# explicitly per-row via verification_note rather than silently presented
# as identical to the exact PubMed title string. Database-only sources
# (UniProt, PDB, AlphaFold, ChEMBL, Pharos, Open Targets, GTEx, HPA,
# gnomAD, OMIM, ClinGen, MGI, IMPC, STRING) are cited inline in each row's
# `sources` field above rather than repeated here, since they are not
# literature and have no PMID/title/year.
# ---------------------------------------------------------------------------
REFERENCES_ROWS = [
    dict(candidate="USP34", PMID="21383061", title="USP34 regulates Wnt signaling by controlling the steady-state levels of Axin (paraphrased from abstract; exact PubMed title string not independently retyped)", journal="unspecified (paraphrased title; do not quote as exact)", year=2011, verification_note="PMID/finding verified against primary text (catalytic mutagenesis, AXIN1 binding); title is paraphrased, not the exact PubMed string"),
    dict(candidate="USP34", PMID="35588869", title="Structural Insights into the Catalytic Mechanism and Ubiquitin Recognition of USP34", journal="Journal of Molecular Biology", year=2022, verification_note="title independently confirmed via Codex's PubMed check during phase review; source of PDB 7W3R, 7W3U"),
    dict(candidate="USP34", PMID="30181118", title="Ubiquitin-specific protease USP34 controls osteogenic differentiation and bone formation by regulating BMP2 signaling", journal="EMBO Journal", year=2018, verification_note="title independently confirmed via Codex's PubMed check during phase review; erratum EMBO J 2020;39(20):e105578"),
    dict(candidate="USP34", PMID="42315110", title="USP34 Haploinsufficiency as a Cause of Neurodevelopmental Phenotypes", journal="Clinical Genetics", year=2026, verification_note="title independently confirmed via Codex's PubMed check during phase review; PMCID PMC13432289"),
    dict(candidate="VEZF1", PMID="29970794", title="Characterization of Small Molecules Inhibiting the Pro-Angiogenic Activity of the Zinc Finger Transcription Factor Vezf1", journal="Molecules", year=2018, verification_note="title independently confirmed via Codex's PubMed check during phase review; PMCID PMC6100598"),
    dict(candidate="VEZF1", PMID="29794136", title="Vezf1 represses Cited2 (paraphrased from abstract; exact PubMed title string not independently retyped)", journal="Journal of Biological Chemistry", year=2018, verification_note="PMID/finding verified via direct fetch; title is paraphrased, not the exact PubMed string"),
    dict(candidate="VEZF1", PMID="24280686", title="RhoB controls coordination of adult angiogenesis and lymphangiogenesis following injury by regulating VEZF1-mediated transcription", journal="Nature Communications", year=2013, verification_note="fully verified; PMCID PMC3868161"),
    dict(candidate="VEZF1", PMID="36923254", title="ETV2 and VEZF1 interaction and regulation of the hematoendothelial lineage during embryogenesis", journal="Frontiers in Cell and Developmental Biology", year=2023, verification_note="fully verified"),
    dict(candidate="VEZF1", PMID="36241701", title="VEZF1/PAQR4 axis in hepatocellular carcinoma, STUB1-mediated VEZF1 degradation (paraphrased from abstract; exact PubMed title string not independently retyped)", journal="Cancer Gene Therapy", year=2022, verification_note="secondary-sourced -- primary full text paywalled, findings from search-indexed content only; title is paraphrased"),
    dict(candidate="VEZF1", PMID="31064890", title="miR-191 targets VEZF1 in an ischemic-stroke angiogenesis model (paraphrased from abstract; exact PubMed title string not independently retyped)", journal="Aging (Albany NY)", year=2019, verification_note="PMID/finding verified; title is paraphrased, not the exact PubMed string"),
    dict(candidate="VEZF1", PMID="15882861", title="Dosage-dependent requirement for mouse Vezf1 in vascular system development", journal="Developmental Biology", year=2005, verification_note="fully verified via EuropePMC; cross-confirmed via MGI:1313291"),
    dict(candidate="VEZF1", PMID="9986727", title="Vezf1: A Zn finger transcription factor restricted to endothelial cells and their precursors", journal="Developmental Biology", year=1999, verification_note="fully verified via EuropePMC; gene-discovery/expression paper, not a knockout paper"),
    dict(candidate="VEZF1", PMID="33231681", title="VEZF1 binding to G-quadruplex DNA at the VASH1 locus modulates alternative polyadenylation (paraphrased from abstract; exact PubMed title string not independently retyped)", journal="Nucleic Acids Research", year=2020, verification_note="PMID/finding verified; title is paraphrased, not the exact PubMed string"),
    dict(candidate="EML5", PMID="15225882", title="Isolation and mapping of an EML family member (rat ortholog description)", journal="Gene", year=2004, verification_note="fully verified; O'Connor et al., describes 11 WD40 + 3 HELP domain architecture in rat ortholog"),
    dict(candidate="EML5", PMID="31406157", title="SPATA7-PTPN21-ZC3H14-EML5-TTC8 locus association study", journal="Scientific Reports", year=2019, verification_note="fully verified; source of genomic-neighbor caveat on STRING associations"),
    dict(candidate="EML5", PMID="17625570", title="Identification of the transforming EML4-ALK fusion gene in non-small-cell lung cancer", journal="Nature", year=2007, verification_note="fully verified; family-level context only, NOT EML5-specific"),
    dict(candidate="EML5", PMID="26730336", title="EML5 expression in resected epileptic neocortex", journal="Iranian Journal of Basic Medical Sciences", year=2015, verification_note="fully verified; PMCID PMC4686571; diseased (epilepsy), not normal, human tissue"),
    dict(candidate="EML5", PMID="35478957", title="A Non-Synonymous Point Mutation in a WD-40 Domain Repeat of EML5 Leads to Decreased Bovine Sperm Quality and Fertility", journal="Frontiers in Cell and Developmental Biology", year=2022, verification_note="fully verified; explicitly bovine (Bos taurus), not human"),
    dict(candidate="CITED2", PMID="12778114", title="Structural basis for negative regulation of hypoxia-inducible factor-1alpha by CITED2", journal="Nature Structural Biology", year=2003, verification_note="fully verified; source of PDB 1P4Q"),
    dict(candidate="CITED2", PMID="14594809", title="Interaction of the TAZ1 domain of the CREB-binding protein with the activation domain of CITED2: regulation by competition between intrinsically unstructured ligands for non-identical binding sites", journal="Journal of Biological Chemistry", year=2004, verification_note="exact title independently confirmed via Codex's PubMed check during phase review; source of PDB 1R8U"),
    dict(candidate="CITED2", PMID="16287139", title="CITED2 in-frame variants and congenital heart defects, VSD2/ASD8 (paraphrased from abstract; exact PubMed title string not independently retyped)", journal="Human Mutation", year=2005, verification_note="PMID/finding verified via UniProt DISEASE CC cross-reference (Sperling et al.); title is paraphrased, not the exact PubMed string"),
    dict(candidate="CITED2", PMID="19951693", title="Cited2 is an essential regulator of adult hematopoietic stem cells", journal="Cell Stem Cell", year=2009, verification_note="fully verified; mouse conditional KO + human confirmatory data"),
    dict(candidate="CITED2", PMID="22308296", title="HIF-1alpha deletion partially rescues defects of hematopoietic stem cell quiescence caused by Cited2 deficiency", journal="Blood", year=2012, verification_note="fully verified; mouse Mx1-Cre conditional KO"),
    dict(candidate="CITED2", PMID="34715054", title="CITED2 coordinates key hematopoietic regulatory pathways to maintain the HSC pool in both steady-state hematopoiesis and transplantation", journal="Stem Cell Reports", year=2021, verification_note="fully verified; mouse hematopoietic-specific deletion"),
    dict(candidate="CITED2", PMID="25184385", title="CITED2-mediated human hematopoietic stem cell maintenance is critical for acute myeloid leukemia", journal="Leukemia", year=2015, verification_note="fully verified; human primary CD34+ cells, NSG xenografts, primary AML"),
    dict(candidate="CITED2", PMID="33288951", title="Stepwise cell fate decision pathways during osteoclastogenesis at single-cell resolution", journal="Nature Metabolism", year=2020, verification_note="fully verified; mouse single-cell transcriptomics + in vivo lineage-specific deletion; direct causal osteoclast-differentiation finding"),
    dict(candidate="CITED2", PMID="19607804", title="Identification of CITED2 as a negative regulator of fracture healing", journal="Biochemical and Biophysical Research Communications", year=2009, verification_note="fully verified; rat mandibular osteotomy model"),
    dict(candidate="CITED2", PMID="12960175", title="CITED2-mediated regulation of MMP-1 and MMP-13 in human chondrocytes under flow shear", journal="Journal of Biological Chemistry", year=2003, verification_note="fully verified; human chondrocyte cell line C-28/I2"),
    dict(candidate="CITED2", PMID="20826544", title="Physiological loading of joints prevents cartilage degradation through CITED2", journal="FASEB Journal", year=2011, verification_note="fully verified; rat in vivo + human chondrocytes in vitro"),
    dict(candidate="CITED2", PMID="33691475", title="Downregulation of SUV39H1 and CITED2 Exerts Additive Effect on Promoting Adipogenic Commitment of Human Mesenchymal Stem Cells", journal="Stem Cells and Development", year=2021, verification_note="fully verified; human MSCs"),
    dict(candidate="CITED2", PMID="11694877", title="Cardiac malformations, adrenal agenesis, neural crest defects and exencephaly in mice lacking Cited2, a new Tfap2 co-activator", journal="Nature Genetics", year=2001, verification_note="fully verified; constitutive mouse knockout, embryonic lethal"),
    dict(candidate="CITED2", PMID="12149478", title="The essential role of Cited2, a negative regulator for HIF-1alpha, in heart development and neurulation", journal="Proceedings of the National Academy of Sciences", year=2002, verification_note="fully verified; constitutive mouse knockout"),
    dict(candidate="CITED2", PMID="15750185", title="Cited2 is required both for heart morphogenesis and establishment of the left-right axis in mouse development", journal="Development", year=2005, verification_note="fully verified; constitutive mouse knockout"),
]
