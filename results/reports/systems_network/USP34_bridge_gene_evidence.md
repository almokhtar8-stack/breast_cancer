# USP34 bridge-gene evidence audit

Focused audit of the four genes that form USP34's shortest (2-edge) network
paths to CTNNB1, PTEN, EP300, and SOX2
(`results/tables/systems_network/USP34_shortest_paths.tsv`): **USP9X,
RPS27A, UBC, UBB**.

Source data: only already-frozen systems-network outputs were read --
`results/networks/systems_network/cytoscape/network_nodes.tsv`, the frozen
Phase 2 per-dataset ranked-gene tables
(`results/tables/systems_network/{gse118713,gse240112,gse111151,gse245601}_ranked_genes.tsv`),
and `data/processed/labels.parquet` (the frozen CRISPR label source). Nothing
was rerun; nothing under `results/networks/systems_network/` was modified.
Full numbers: `results/tables/systems_network/USP34_bridge_gene_evidence.tsv`.
Figure: `results/figures/systems_network/final_review/08_USP34_bridge_gene_evidence.png`.

GSE118713's raw DE table carries three contrasts (TAMR_vs_MCF7,
FASR_vs_MCF7, TAMR_vs_FASR); the frozen ranked-genes table already selects
TAMR_vs_MCF7 only (the resistance contrast used throughout this project),
and that is the value reported below.

CRISPR sign convention: negative `effect_size` = sensitising_KO, positive =
tolerance_associated_KO. No value below is FDR<0.05 unless marked.
**GSE245601 is the acute 12h response and is reported as a fifth layer for
completeness only -- it is never counted toward "resistance evidence" and
never counted toward the classification below.**

## 1. Per-gene evidence table

| Gene | CRISPR effect (p / FDR) | GSE118713 log2FC (p / FDR) | GSE240112 log2FC (p / FDR) | GSE111151 log2FC (p / FDR) | GSE245601 acute log2FC (p / FDR) |
|---|---|---|---|---|---|
| **USP9X** | -0.188 (p=0.531, FDR=0.914) | +0.543 (p=0.042*, FDR=0.106) | -0.174 (p=0.416, FDR=0.559) | -0.119 (p=0.490, FDR=0.823) | -0.176 (p=0.055, FDR=0.371) |
| **UBB** | +1.170 (p=0.015*, FDR=0.585) | +0.105 (p=0.536, FDR=0.674) | -0.624 (p=0.093, FDR=0.192) | -0.177 (p=0.356, FDR=0.745) | +0.328 (p=0.036*, FDR=0.313) |
| **RPS27A** | -0.745 (p=0.285, FDR=0.827) | +0.324 (p=0.063, FDR=0.144) | +0.645 (p=0.182, FDR=0.310) | -0.090 (p=0.491, FDR=0.823) | +0.006 (p=0.958, FDR=0.987) |
| **UBC** | +0.583 (p=0.106, FDR=0.762) | +0.238 (p=0.095, FDR=0.197) | +0.351 (p=0.370, FDR=0.514) | +0.066 (p=0.698, FDR=0.911) | +0.332 (p=0.011*, FDR=0.208) |

\* nominal p<0.05, not FDR<0.05. **No gene reaches FDR<0.05 in any layer.**

Resistance-direction consistency across the 3 resistance datasets (GSE118713 / GSE240112 / GSE111151):
- USP9X: up / down / down -- **mixed**
- UBB: up / down / down -- **mixed**
- RPS27A: up / up / down -- **mixed**
- UBC: up / up / up -- **consistent direction, but not statistically significant in any of the three**

## 2. Network facts (from USP34_shortest_paths.tsv)

| Bridge gene | Bridges USP34 to | Edge type(s) to USP34 |
|---|---|---|
| USP9X | CTNNB1, SOX2 | functional_association (STRING, 0.793) |
| RPS27A | CTNNB1, PTEN | physical_PPI (STRING, 0.961) |
| UBC | PTEN, EP300 | physical_PPI (STRING, 0.982) |
| UBB | CTNNB1 | physical_PPI (STRING, 0.954) |

## 3. Classification (conservative; A requires FDR<0.05 in CRISPR and/or a resistance dataset; B requires >=1 nominal p<0.05 hit in CRISPR or a resistance dataset; C requires neither)

| Gene | Classification | Rationale |
|---|---|---|
| **USP9X** | **B. PARTIAL SUPPORT** | One nominal hit: GSE118713 up, p=0.042 (not FDR-significant), directly contradicted (non-significantly) by down-trending GSE240112/GSE111151. CRISPR effect near null (p=0.53). No consistent resistance-direction signal. |
| **UBB** | **B. PARTIAL SUPPORT** | One nominal CRISPR hit (p=0.015, tolerance-associated direction) plus a nominal acute-only hit (p=0.036, not counted toward resistance). No nominal or FDR hit in any resistance dataset; resistance direction mixed. |
| **RPS27A** | **C. NETWORK-ONLY / GENERIC BRIDGE** | No nominal (p<0.05) or FDR hit in any layer, resistance or CRISPR. Flattest profile of the four genes. |
| **UBC** | **C. NETWORK-ONLY / GENERIC BRIDGE** | No nominal or FDR hit in CRISPR or any resistance dataset (closest is GSE118713 p=0.095). The only nominal hit is in the acute-only GSE245601 layer, which by design does not count as resistance or CRISPR evidence. Direction is at least consistently up across all three resistance datasets, but none reach even nominal significance. |

No tier above is evidence of mechanism, and no expression association here is evidence of causality -- these are descriptive classifications of how much (if any) of our own CRISPR/resistance data independently backs a STRING-derived bridge.

## 4. Answers to the specific questions

**Which bridge has the strongest independent support from our own data?**
None reach FDR<0.05 anywhere, so "strongest" only distinguishes among weak signals. UBB has the single lowest p-value across all cells for any of the four genes (CRISPR p=0.015) plus a second independent nominal hit (acute p=0.036, same up/tolerance-leaning direction) -- the most internally consistent (if still non-FDR-significant) weak signal. USP9X's one nominal hit (GSE118713 p=0.042) is undercut by being contradicted, non-significantly, by the other two resistance datasets. Neither should be described as data-supported; UBB is marginally the least weak of the four.

**Is USP9X genuinely supported, or mainly a STRING-network hypothesis?**
Mainly a STRING-network hypothesis. USP9X's only signal above noise is a single nominal (not FDR-significant) hit in one of three resistance datasets, contradicted in direction by the other two, and its CRISPR result is essentially null (p=0.53, FDR=0.91). Being a deubiquitinase paralog of USP34 makes it the most *biologically plausible* of the four bridges, but plausibility is not the same as support, and our own data does not currently back it beyond the STRING functional_association edge itself.

**Do UBC/UBB/RPS27A look like generic ubiquitin hubs?**
RPS27A and UBC show no significant or even nominal support in CRISPR or resistance expression -- consistent with them functioning as generic, highly-connected ubiquitin-pathway hub genes in STRING (ubiquitin/ribosomal-fusion and polyubiquitin precursor genes are expected to have dense, largely non-specific interaction degree) rather than tamoxifen-resistance-specific nodes. UBB is slightly different: it has one nominal CRISPR hit, so it is not quite as flat as RPS27A/UBC, but its resistance-dataset evidence is just as absent and directionally mixed -- it is closer to "generic hub with one added functional flag" than to a resistance-specific signal.

**Does any bridge deserve follow-up before literature/druggability analysis?**
Given the conservative classification, no bridge currently qualifies as data-supported (no A-tier gene). If a follow-up priority must be chosen: USP9X is the most biologically motivated candidate to keep on a watch-list (deubiquitinase paralog of USP34, plausible functional redundancy/substrate overlap) despite the current lack of independent data support -- any future look at it should be explicit that support today is network-topology-only. UBB is the second candidate to watch given its nominal CRISPR hit, with the same caveat. RPS27A and UBC show no evidence pattern that distinguishes them from generic ubiquitin-pathway connectivity and do not warrant prioritized follow-up on the current evidence.
