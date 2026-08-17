# Structure / Druggability Figure v1 -- Short Data Note

**Question:** can these four tamoxifen-sensitising candidate
vulnerabilities realistically be targeted? The four are NOT equally
tractable, and that difference is the result.

## 1-2. Evidence source and PDB IDs

All factual claims come from the project's already-audited
`results/tables/post_audit_sensitivity/06b_structural_tractability_audit.tsv`
(columns A-G, `precise_summary`, `sources`), read unmodified. The PDB IDs
displayed are **recovered by scanning that table's own text** at build
time (`pdb_ids_from_audit()`), not hand-asserted; the figure asserts that
each rendered ID appears in the audit text for that gene.

| Gene | PDB shown | Experimental? | Bound species shown |
|---|---|---|---|
| KDM1A | **6NQU** | Yes, X-ray | **GSK2879552** — a selective small-molecule inhibitor (ligand KWM) |
| TLK2 | **5O0Y** | Yes, X-ray 2.86 Å | **ATP-γ-S** (AGS), a non-hydrolysable ATP/substrate analog — **not** an inhibitor or drug-like ligand |
| USP34 | **7W3U** | Yes, X-ray | **Ubiquitin-propargylamide covalent activity-based probe** (warhead AYE, covalently linked to Cys1903, LINK record 1.59 Å) — **not** a drug |
| VEZF1 | none | **No experimental structure exists** | n/a |

Other IDs present in the audit text but deliberately **not** rendered:
KDM1A 2Z5U (tranylcypromine-class covalent inhibitor), USP34 7W3R (apo
catalytic domain), and VEZF1's 1AAY — which is the **unrelated Zif268
zinc finger** used as a homology-model template in the published VEZF1
inhibitor study, *not* a VEZF1 structure. Structures came from the
already-downloaded local files with recorded provenance
(`usp34_structures/`, `kdm1a_tlk2_structures/`); nothing was
re-downloaded or substituted.

## 3-5. Ligand / probe / inhibitor distinctions (kept strictly separate)

- **Crystallographic ligand or substrate analog** (TLK2's ATP-γ-S): proves
  the pocket is occupiable by a nucleotide mimetic. It is *not* an
  inhibitor and implies no therapeutic chemistry.
- **Covalent activity-based probe** (USP34's UbPA): experimentally proves
  the catalytic cysteine (Cys1903) is reactive and covalently
  addressable. This is a **tractability hypothesis**, not a drug, and not
  proof of druggability.
- **Validated selective inhibitor** (KDM1A's GSK2879552): a real
  selective small molecule with a co-crystal structure.
- **Clinical-stage pharmacology** (KDM1A only): the audit records 9 LSD1
  inhibitors reported to have reached clinical stage (including
  iadademstat/ORY-1001); **none is breast-cancer-approved**.

## 6. Evidence supporting each pharmacology statement

Derived programmatically from the audit's own columns
(`evidence_track()`), never typed in:

| Gene | Experimental structure | Direct ligand/probe | Selective or clinical pharmacology |
|---|---|---|---|
| KDM1A | established | established | established |
| TLK2 | limited (audit: PARTIAL, kinase domain only) | limited (ATP analog, not a drug) | none |
| USP34 | limited (audit: PARTIAL, catalytic domain ≈12% of protein) | limited (probe, not a drug) | none |
| VEZF1 | none | none | none |

Filled = established, half = limited/not-a-drug, open = none. TLK2 and
USP34 deliberately occupy the **same** levels rather than being ranked
against each other: they differ in *kind*, not degree (TLK2 = strong
structural class with no selective chemistry; USP34 = direct covalent
catalytic-cysteine evidence with no selective chemistry). No numerical
maturity score was invented — the project has no validated scoring
framework for this, so a 3-row evidence track was used instead of a
single continuum.

VEZF1's "none" for pharmacology is not "nothing exists": the audit
records a real but weak, homology-model-derived screening hit
(T4/503-1-83, IC50 ≈ 20 µM, DNA-binding inhibition), shown as the
panel's small secondary annotation. It is far from validated or
selective.

## 7. Why VEZF1 has no structure panel

The audit records `A_experimental_human_structure_exists = False`: no
solved zinc-finger-domain structure exists, and the published inhibitor
study built a homology model on the unrelated Zif268 finger (1AAY) for
exactly that reason. A homology or AlphaFold model would be visually
indistinguishable from the three experimental panels while being a
categorically weaker kind of evidence, so none is shown. Instead the
panel carries a restrained six-zinc-finger domain schematic (from the
protein's annotated architecture) plus the explicit statement "No
experimental structure" — the absence is the finding.

## 8. No docking performed

No docking, pose prediction, binding-affinity estimation, energy
minimisation, or pocket re-detection was performed in this figure or its
render script. The panels are static renders of deposited experimental
coordinates. No AlphaFold model is used anywhere.

## 9-10. Scope and limits

No candidate ranking, CRISPR value, or frozen scientific result was
changed; this is one translational evidence layer, not a final candidate
ranking. **Structural and pharmacological tractability does not imply
efficacy** in ER+ breast cancer or against tamoxifen resistance: an
existing selective or clinical-stage inhibitor (KDM1A) means the target
class is chemically addressable, not that the drug works in this
indication — the audit itself notes none of the LSD1 clinical compounds
is breast-cancer-approved, and that KDM1A's novelty as a tamoxifen-response
modifier is LOW because its ER-corepressor role is established prior
literature.
