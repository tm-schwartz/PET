# I am validating a manually created composite MNI-space reference mask for FDG-PE

**Question:** I am validating a manually created composite MNI-space reference mask for FDG-PET (pons + cerebellar vermis) against published volumetry. The mask was built from public atlases: (a) pons from the NextBrain histological MNI152 1mm atlas (Casamitjana et al., Nature 2025; label 580/10580 'pons', measured 11.8 cm3 total, lateral extent x = ±19 mm, z from -51 to -15), and (b) vermis from the Diedrichsen/SUIT 2009 probabilistic cerebellar atlas deterministic MNI labels 6,9,12,15,18,21,24,27 (Vermis VI, CrusI, CrusII, VIIb, VIIIa, VIIIb, IX, X; measured 7.2 cm3 on a 1mm MNI grid). The intended provenance equivalence is FreeSurfer/FastSurfer: FastSurfer cerebnet (cerebellar subnuclei incl. vermis) and FIRST/Bayesian brainstem pons label.

Please answer with cited evidence:
1. Does the FreeSurfer 'pons' brainstem subfield label (Iglesias et al. 2015, NeuroImage 113:184-195, and the FIRST brainstem extension) INCLUDE the middle cerebellar peduncle (MCP/brachium pontis)? What anatomical boundaries does its delineation protocol use (superior/inferior/lateral limits)?
2. Published volume norms from that exact pipeline: mean pons, midbrain, medulla, and whole-brainstem volumes in healthy adults (UK Biobank Elvsåshagen 2020, Bergsland, Fujimoto 2023, Müller 2023). I have pons ≈ 14.5–15.8 cm3 in native subject space — confirm or correct with specific numbers.
3. Total human cerebellar VERMIS volume (all vermal lobules combined) in vivo: what do segmentation studies using the SUIT/Diedrichsen atlas report (cm3)? E.g., Bocchetta 2016 (genetic FTD, Diedrichsen atlas label fusion) control vermis volume; Bernard 2013; any pediatric/adult morphometry papers. What is the expected volume of the vermis in the MNI152 template itself?
4. NextBrain atlas (Casamitjana et al. Nature 2025): exact DOI/title, what segmentation it provides for pons (does the pons label include MCP or pontine tegmentum separately), and any volume statistics reported.
5. Is 11.8 cm3 plausible for a pons label on an MNI template vs 14.5-15.8 cm3 in native MRI space — i.e., do MNI-template structures systematically appear smaller than in vivo individual brains (MNI152 nonlinear average shrinkage of the brainstem)? Any measurements comparing brainstem/pons volume in MNI152 template vs subjects?
Report exact numbers with citations (DOI + quoted table values where possible).

# Your **MNI-space mask** is plausible, but published volumetry does not establish **pipeline equivalence** or systematic template shrinkage.

The main validation issue is **anatomical definition, not volume alone**: FreeSurfer’s pons has explicit peduncular cutoffs, SUIT’s named vermis labels do not establish coverage of all vermal lobules, and NextBrain’s whole-region validation does not establish equivalence to FreeSurfer brainstem subfields [[1.1]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/) [[2.1]](https://consensus.app/papers/cerebnet-a-fast-and-reliable-deeplearning-pipeline-for-faber-kügler/a60b5a9f54235f8eb88503a55027580d/) [[3.1]](https://consensus.app/papers/a-probabilistic-histological-atlas-of-the-human-brain-for-casamitjana-mancini/1430631409aa5b51b0461df2a6e2e5a5/).

Some requested numerical tables are absent from the supplied excerpts. The numbers below distinguish directly recoverable measurements from benchmarks that cannot be verified here.

## 1. FreeSurfer Pons Boundaries

**FreeSurfer’s protocol includes some MCP-region tissue within a bounded pontine segmentation; it does not establish inclusion of the entire MCP/brachium pontis.** The strongest evidence is Appendix A of Iglesias et al., *Bayesian segmentation of brainstem structures in MRI*, DOI **10.1016/j.neuroimage.2015.02.065**, which explicitly instructs tracing around “the bright pons and middle cerebellar peduncle regions” [[1.2]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/).

Crucially, that instruction operates alongside a posterior cutoff: where the middle and superior cerebellar peduncles contact the pons, the protocol draws a straight line between their superior and inferior branching points [[1.1]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/). Therefore, neither “MCP completely excluded” nor “whole MCP included” accurately summarizes the published delineation.

| Boundary | Published FreeSurfer delineation |
|---|---|
| Superior | Superior pontine notch to inferior quadrigeminal-plate edge, with lateral anatomical adjustments [[1.3]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/) |
| Inferior | Straight sagittal boundary, subsequently linked to the inferior pontine notch and axial posterior boundary [[1.4]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/) |
| Posterior, medial | Tissue–CSF interface along the fourth ventricle [[1.5]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/) |
| Posterior, peduncular | Straight line joining superior and inferior peduncular branching points [[1.6]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/) |
| Anterior/lateral surface | Tissue–CSF interface, excluding exiting vessels and nerves [[1.7]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/) |

**Figure 1:** Anatomical landmarks defining the FreeSurfer pontine segmentation boundaries

The protocol does **not** prescribe a fixed MNI lateral limit such as ±19 mm. It uses anatomical landmarks after rigid alignment to FreeSurfer reference space, then returns annotations to original space [[1.8]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/) [[1.9]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/). Its inferior delineation also follows the fourth-ventricular tissue boundary past the vestibular-nuclear notch, rather than segmenting only the ventral pontine bulge [[1.10]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/).

**FIRST should not be treated as synonymous with this Bayesian subfield pipeline.** Iglesias distinguishes FSL’s whole-brainstem segmentation from the newly introduced four-structure method [[1.11]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/) [[1.12]](https://consensus.app/papers/bayesian-segmentation-of-brainstem-structures-in-mri-iglesias-leemput/71adc12849cd52818b4cc64a87d912f6/). An intermethod comparison found substantial FIRST differences, including fourth-ventricle inclusion, and poor agreement with other segmentations [[4.1]](https://consensus.app/papers/brainstem-enlargement-in-preschool-children-with-autism-bosco-giuliano/d8c6198df0c85196ae989aa8751f3747/) [[4.2]](https://consensus.app/papers/brainstem-enlargement-in-preschool-children-with-autism-bosco-giuliano/d8c6198df0c85196ae989aa8751f3747/). No identifiable “FIRST brainstem extension” protocol in the supplied evidence resolves its specific MCP boundary.

This ambiguity is not unique to FreeSurfer: Nigro’s separate LABS method explicitly reports no shared anatomical criterion for separating pons from MCP and introduces its own landmark [[5.1]](https://consensus.app/papers/fully-automated-segmentation-of-the-pons-and-midbrain-nigro-cerasa/60a59bd71d11541eaf616fc10ec11cb1/) [[5.2]](https://consensus.app/papers/fully-automated-segmentation-of-the-pons-and-midbrain-nigro-cerasa/60a59bd71d11541eaf616fc10ec11cb1/).

## 2. Native-Space Brainstem Volumes

**Your 14.5–15.8 cm³ pons benchmark is compatible with one recoverable healthy-control result, but cannot be confirmed as the range of means across all the named studies.** Rowe et al. report a healthy-control mean of **15.45714 cm³**, with substantial between-subject variability and only seven controls [[6.1]](https://consensus.app/papers/brainstem-substructure-atrophy-in-lateonset-rowe-rangaprakash/dfea6d85a172548bb35fd2ccd39ab896/) [[6.2]](https://consensus.app/papers/brainstem-substructure-atrophy-in-lateonset-rowe-rangaprakash/dfea6d85a172548bb35fd2ccd39ab896/).

The following are the exact absolute measurements recoverable from the supplied evidence; disease-cohort measurements are explicitly separated from healthy norms.

| Source and population | Structure | Quoted measurement, mm³ | Converted mean ± SD, cm³ |
|---|---|---:|---:|
| Rowe 2025, seven controls | Pons | “15,457.14 ± 2,748.41” | **15.457 ± 2.748** [[6.1]](https://consensus.app/papers/brainstem-substructure-atrophy-in-lateonset-rowe-rangaprakash/dfea6d85a172548bb35fd2ccd39ab896/) |
| Sander 2019, 80 MS patients, FreeSurfer | Pons | “14,005.89”; SD “2,084.54” | 14.006 ± 2.085 [[7.1]](https://consensus.app/papers/accurate-rapid-and-reliable-fully-automated-mri-brainstem-sander-pezold/9be7f82b50e55452afa780487534a87d/) |
| Same MS cohort | Midbrain | “5,803.48”; SD “603.40” | 5.803 ± 0.603 [[7.1]](https://consensus.app/papers/accurate-rapid-and-reliable-fully-automated-mri-brainstem-sander-pezold/9be7f82b50e55452afa780487534a87d/) |
| Same MS cohort | Medulla | “4,437.50”; SD “541.35” | 4.438 ± 0.541 [[7.1]](https://consensus.app/papers/accurate-rapid-and-reliable-fully-automated-mri-brainstem-sander-pezold/9be7f82b50e55452afa780487534a87d/) |
| Same MS cohort | Brainstem | “24,246.86”; SD “3,041.96” | 24.247 ± 3.042 [[7.1]](https://consensus.app/papers/accurate-rapid-and-reliable-fully-automated-mri-brainstem-sander-pezold/9be7f82b50e55452afa780487534a87d/) |

**Figure 2:** Recoverable absolute brainstem measurements with healthy and disease cohorts distinguished

Rowe’s DOI is **10.1007/s12311-025-01803-4**; Sander’s is **10.1002/hbm.24687**. Sander’s table is useful for checking pipeline magnitude, **not** for establishing healthy-adult norms. Another FreeSurfer study provides a healthy female-control medulla measurement of **3,236 ± 411 mm³**, illustrating why population and processing details must accompany quoted values [[8.1]](https://consensus.app/papers/volumetric-alteration-of-brainstem-in-female-migraineurs-he-kis-jakab/40fa4e8403b358518726c274aa98ed80/).

Elvsåshagen’s *The genetic architecture of human brainstem structures and their involvement in common brain disorders*, DOI **10.1038/s41467-020-17376-1**, confirms FreeSurfer 6.0 plus Bayesian brainstem segmentation, with extensive visual QC [[9.1]](https://consensus.app/papers/the-genetic-architecture-of-human-brainstem-structures-elvsåshagen-bahrami/945ce3ba48495fa08a548c521df58b31/) [[9.2]](https://consensus.app/papers/the-genetic-architecture-of-human-brainstem-structures-elvsåshagen-bahrami/945ce3ba48495fa08a548c521df58b31/). However, its supplied excerpts do **not** contain the requested absolute mean-volume table. The relevant Bergsland, Fujimoto 2023, and Müller 2023 numerical tables likewise cannot be identified and quoted from the supplied indexed evidence.

Whole-brainstem definitions also require checking: a later large genetic study explicitly distinguishes Bayesian total brainstem **including SCP** from standard subcortical brainstem **excluding SCP** [[10.1]](https://consensus.app/papers/the-genetic-architecture-of-brainstem-structures-xue-fu/12c9f84e9cee5097bcbd6eb07e52fddd/) [[10.2]](https://consensus.app/papers/the-genetic-architecture-of-brainstem-structures-xue-fu/12c9f84e9cee5097bcbd6eb07e52fddd/). Consequently, apparently matching “whole-brainstem” measurements need not measure identical anatomy.

## 3. Vermis Volume and Coverage

**The clearest Diedrichsen-based control measurement is 5.034 ± 0.375 cm³—not 7.2 cm³.** Bocchetta et al., *Patterns of regional cerebellar atrophy in genetic frontotemporal dementia*, DOI **10.1016/j.nicl.2016.02.008**, report Table 1 control vermis values as **“5034 375” mm³**, meaning mean and SD, in 18 controls [[11.1]](https://consensus.app/papers/patterns-of-regional-cerebellar-atrophy-in-genetic-bocchetta-cardoso/cf1c34962680544ebf04d4269a827f7d/).

They used atlas propagation and label fusion of the Diedrichsen cerebellar atlas [[11.2]](https://consensus.app/papers/patterns-of-regional-cerebellar-atrophy-in-genetic-bocchetta-cardoso/cf1c34962680544ebf04d4269a827f7d/). However, all volumes were corrected for total intracranial volume; the table should therefore not be presented without qualification as an unadjusted native-space anatomical norm [[11.3]](https://consensus.app/papers/patterns-of-regional-cerebellar-atrophy-in-genetic-bocchetta-cardoso/cf1c34962680544ebf04d4269a827f7d/).

| Measurement or definition | What the evidence establishes |
|---|---|
| Bocchetta control vermis | **5.034 ± 0.375 cm³**, with TIV correction [[11.1]](https://consensus.app/papers/patterns-of-regional-cerebellar-atrophy-in-genetic-bocchetta-cardoso/cf1c34962680544ebf04d4269a827f7d/) |
| Your eight-label union | Named VI–X labels do not demonstrate inclusion of anterior vermis I–V |
| Complete anterior vermis | Requires an operational midline/paramedian definition because separation is indistinct [[12.1]](https://consensus.app/papers/cerebellum-development-during-childhood-and-adolescence-tiemeier-lenroot/eb644698e2f653b0b7788471dff11bd1/) |
| SUIT versus CerebNet | Protocols differ in vermal aggregation and white-matter handling [[2.1]](https://consensus.app/papers/cerebnet-a-fast-and-reliable-deeplearning-pipeline-for-faber-kügler/a60b5a9f54235f8eb88503a55027580d/) |
| MNI maximum-probability ROI | An ROI derived from normalized probability maps, not a reported population mean [[13.1]](https://consensus.app/papers/a-probabilistic-mr-atlas-of-the-human-cerebellum-diedrichsen-balsters/c62234fca94f524bbb0dffaed2e96433/) |

**Figure 3:** Differences between atlas-defined vermis masks and anatomical vermis measurements

Your stated label list contains no separately designated vermis I–V. That omission matters for a claim of **“all vermal lobules combined”**: longitudinal morphometry explicitly notes that anterior vermis must be operationally defined as midline and paramedian anterior-lobe sectors [[12.1]](https://consensus.app/papers/cerebellum-development-during-childhood-and-adolescence-tiemeier-lenroot/eb644698e2f653b0b7788471dff11bd1/). It is more defensible to call your mask the **union of the specified SUIT vermal labels** than the complete anatomical vermis.

CerebNet is also not an exact replacement for that union. It separates gray and white matter and adds detailed white-matter strands that earlier protocols did not capture [[2.2]](https://consensus.app/papers/cerebnet-a-fast-and-reliable-deeplearning-pipeline-for-faber-kügler/a60b5a9f54235f8eb88503a55027580d/) [[2.3]](https://consensus.app/papers/cerebnet-a-fast-and-reliable-deeplearning-pipeline-for-faber-kügler/a60b5a9f54235f8eb88503a55027580d/). Its SUIT-derived pretraining labels required additional FreeSurfer cerebellar white matter because SUIT did not supply that segmentation [[2.4]](https://consensus.app/papers/cerebnet-a-fast-and-reliable-deeplearning-pipeline-for-faber-kügler/a60b5a9f54235f8eb88503a55027580d/). “Cerebellar lobular/vermal subsegments” is therefore more accurate here than “subnuclei.”

The pediatric/adult studies reinforce the definition problem rather than supply interchangeable norms. Hodgdon measures I–V and VI–VII, including gray and white matter, but not the entire I–X vermis [[14.1]](https://consensus.app/papers/mri-morphometry-of-the-anterior-and-posterior-cerebellar-hodgdon-anderson/e34e842d25c55748bdcee66dbd2cdcd3/) [[14.2]](https://consensus.app/papers/mri-morphometry-of-the-anterior-and-posterior-cerebellar-hodgdon-anderson/e34e842d25c55748bdcee66dbd2cdcd3/). Tiemeier measures three vermal regions longitudinally across ages 5–24 [[12.2]](https://consensus.app/papers/cerebellum-development-during-childhood-and-adolescence-tiemeier-lenroot/eb644698e2f653b0b7788471dff11bd1/). The supplied Bernard connectivity study reports coordinates and networks, not an absolute total-vermis norm [[15.1]](https://consensus.app/papers/resting-state-corticocerebellar-functional-bernard-seidler/1ee7afd4d0055d97921dcffb94d301f8/).

No supplied excerpt establishes an expected **MNI152-template vermis volume of 7.2 cm³**. Diedrichsen’s atlas provides probabilistic and maximum-probability spatial assignments from 20 young adults, not a directly interchangeable template/native volumetric standard [[13.2]](https://consensus.app/papers/a-probabilistic-mr-atlas-of-the-human-cerebellum-diedrichsen-balsters/c62234fca94f524bbb0dffaed2e96433/) [[13.3]](https://consensus.app/papers/a-probabilistic-mr-atlas-of-the-human-cerebellum-diedrichsen-balsters/c62234fca94f524bbb0dffaed2e96433/).

## 4. NextBrain Provenance

The exact publication is **Casamitjana et al., “A probabilistic histological atlas of the human brain for MRI segmentation,” Nature, DOI 10.1038/s41586-025-09708-2**. It provides a probabilistic atlas built from five hemispheres, with 333 ROIs and a companion Bayesian MRI segmentation tool [[3.2]](https://consensus.app/papers/a-probabilistic-histological-atlas-of-the-human-brain-for-casamitjana-mancini/1430631409aa5b51b0461df2a6e2e5a5/) [[3.3]](https://consensus.app/papers/a-probabilistic-histological-atlas-of-the-human-brain-for-casamitjana-mancini/1430631409aa5b51b0461df2a6e2e5a5/).

| Provenance question | Verified answer |
|---|---|
| Published atlas resolution | Probabilities on a **0.2-mm isotropic grid** [[3.4]](https://consensus.app/papers/a-probabilistic-histological-atlas-of-the-human-brain-for-casamitjana-mancini/1430631409aa5b51b0461df2a6e2e5a5/) |
| Relationship to Allen MNI labels | Explicitly distinct, with 333 versus 138 ROIs [[3.5]](https://consensus.app/papers/a-probabilistic-histological-atlas-of-the-human-brain-for-casamitjana-mancini/1430631409aa5b51b0461df2a6e2e5a5/) |
| FreeSurfer comparison | Whole-region, silver-standard validation—not pons-subfield equivalence [[3.6]](https://consensus.app/papers/a-probabilistic-histological-atlas-of-the-human-brain-for-casamitjana-mancini/1430631409aa5b51b0461df2a6e2e5a5/) |
| Pons-related terminology | Preprint excerpt includes “Pons-nuc,” without a volumetric interpretation [[16.1]](https://consensus.app/papers/a-nextgeneration-histological-atlas-of-the-human-brain-and-casamitjana-mancini/fed4983563dc54f5b0bfacf842465247/) |
| Exact 580/10580 boundaries and volume | Not documented in the supplied excerpts |

**Figure 4:** Verified NextBrain provenance and unresolved pontine label correspondence questions

Thus, the **1-mm MNI file you used needs its own release-level provenance check**. The supplied publication does not verify that label IDs 580/10580 represent the FreeSurfer-equivalent complete pons, nor whether MCP and pontine tegmentum are separately assigned. The presence of a pons-related name is insufficient to establish those anatomical relationships.

The accelerated NextBrain work explicitly acknowledges protocol differences when matching atlas ROIs and describes merging subregions to **approximate** FreeSurfer whole-region definitions [[17.1]](https://consensus.app/papers/fast-segmentation-with-the-nextbrain-histological-atlas-puonti-nolan/61715fb65ba95612aedb50f5dd4e303f/) [[17.2]](https://consensus.app/papers/fast-segmentation-with-the-nextbrain-histological-atlas-puonti-nolan/61715fb65ba95612aedb50f5dd4e303f/). Its supplied volume comparisons concern hippocampal, thalamic, and amygdalar substructures—not a pontine normative table [[17.3]](https://consensus.app/papers/fast-segmentation-with-the-nextbrain-histological-atlas-puonti-nolan/61715fb65ba95612aedb50f5dd4e303f/).

## 5. Template Size and Validation

**An 11.8-cm³ template label is not inherently implausible, but systematic MNI brainstem shrinkage is not an evidence-supported explanation.** Relative to your proposed native means of 14.5–15.8 cm³, it is approximately **18.6–25.3% smaller**; these percentages are arithmetic comparisons, not measured normalization effects.

| Relevant comparison | Finding and limitation |
|---|---|
| Hammers, native versus MNI normalization | Total regional volume increased **46%**, not decreased; not a pons-specific result [[18.1]](https://consensus.app/papers/three‐dimensional-maximum-probability-atlas-of-the-human-hammers-allom/3f5daa1554c1562a999eff1a3e1ff985/) |
| Giff, native versus MNI152 | Several subcortical measurements increased after transformation; not brainstem volumetry [[19.1]](https://consensus.app/papers/spatial-normalization-discrepancies-between-native-and-giff-noren/2e0a280bc1025fe49da78099c1cbda7b/) |
| SUIT versus whole-brain normalization | Cerebellar fissure spatial spread decreased **60%**; alignment, not shrinkage [[20.1]](https://consensus.app/papers/a-spatially-unbiased-atlas-template-of-the-human-diedrichsen/33e599a6e9a9596f905e9a271c9afc6a/) |
| Elderly versus MNI152 template | Brainstem overlap differences were minimal; no pontine volume scaling factor [[21.1]](https://consensus.app/papers/developing-and-validating-an-elderly-brain-template-a-ota-nakazato/c7916a686a845ed6847ba69ff03f96f3/) |

**Figure 5:** Normalization evidence does not establish systematic pontine template shrinkage

Hammers’ DOI is **10.1002/hbm.10123**; Giff’s is **10.1016/j.pscychresns.2023.111595**. Neither provides the requested matched comparison between an MNI152 pontine label and native FreeSurfer pons volumes. Their findings reject a general presumption of shrinkage but do not determine the direction of local deformation in your particular mask.

Boundary differences remain a concrete alternative explanation. For example, Luft’s older volumetry explicitly grouped **pons, pontine tegmentum, and cerebellar peduncles**, whereas Lambert’s brainstem mask included only parts of the middle and superior peduncles [[22.1]](https://consensus.app/papers/patterns-of-agerelated-shrinkage-in-cerebellum-and-luft-skalej/b405b5a5c4905158a5143ecb6a9528e2/) [[23.1]](https://consensus.app/papers/multiparametric-brainstem-segmentation-using-a-modified-lambert-lutti/9f2744c85ca1566fb8d917a43c635c93/). Conversely, later brainstem methods acknowledge sufficiently different delineation protocols that direct FreeSurfer overlap comparisons were not performed [[24.1]](https://consensus.app/papers/automated-brainstem-parcellation-using-multiatlas-magnusson-love/bea5946d09c753b9a07595c422d7ae44/). Linear measurements and planimetric norms also cannot substitute for volumes: the lower-ranked Metwally study measures diameters and ratios, not cm³ [[25.1]](https://consensus.app/papers/neuroanatomical-mri-study-reference-values-for-the-metwally-basha/776b04e760135541be075b9718380c57/).

For your validation, the defensible distinction is between **volume plausibility** and **demonstrated anatomical correspondence**. The most important shift is from matching headline volumes to matching label definitions; the biggest open question is whether your exact NextBrain-derived pons and SUIT vermis union cover the intended FreeSurfer/CerebNet anatomy. Neither 11.8 cm³ nor 7.2 cm³ alone resolves that question.
 
_These search results were found and analyzed using Consensus, an AI-powered search engine for research. Try it at https://consensus.app. © 2026 Consensus NLP, Inc. Personal, non-commercial use only; redistribution requires copyright holders’ consent._
 
## References
 
`1.` Iglesias J, Van Leemput K, Bhatt P, et al. Bayesian segmentation of brainstem structures in MRI. *NeuroImage.* 2015;113:184 - 195. doi:10.1016/j.neuroimage.2015.02.065 · [Access through your institution](https://libkey.io/libraries/519/articles/54230606/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [1.1] Once the middle and superior cere- bellar peduncles make contact with the pons, draw a straight line from the superior point where the pedunclefirst branches from the pons to the most inferior point where the peduncle branches from the pons. (Methods)
 
> [1.2] In this one coronal slice, trace around the curvature of the bright pons and middle cerebellar peduncle regions andfill in the region. (Methods)
 
> [1.3] Tracing of the superior boundary: in the midsagittal plane,first trace the line passing through the superior pontine notch and the inferior edge of the quadriminal plate, then the line to the quadriminal plate and, after all regions have been traced, erase the extraneous portions of the lines. (Methods)
 
> [1.4] Then, in sagittal slices, trace the inferior boundary as the straight diagonal line that extends anteriorly from the inferior pontine notch to the posterior voxel created by the axial-defined boundary. (Methods)
 
> [1.5] 3. Tracing of the posterior boundary: in sagittal view,first trace the line along the tissue–CSF boundary. (Methods)
 
> [1.6] Then, repeat this step in subsequent sagittal slices. (Methods)
 
> [1.7] Then, in lat- eral slices, trace the inferior boundary as defined by the tissue–CSF interface, without including the blood vessels and nerves that extend from the middle of the pons. (Methods)
 
> [1.8] Delineation protocol for brainstem dataset Rather than delineating the brainstem structures in the native space of the scans directly, these scans werefirst rigidly registered to the FreeSurfer reference space (“fsaverage”). (Methods)
 
> [1.9] The manual annotations were made on the registered scans, which helps reduce the variability in the annotations, and then warped back to the original space using the inverse transform and nearest neighbor interpolation. (Methods)
 
> [1.10] Trace along the CSF–tissue boundary just past the vestibular nuclei (CN VIII), which can be visualized by the notch of the fourth ventricle bound- ary, which becomes a vertical line. (Methods)
 
> [1.11] The segmentation algorithm (Patenaude et al., 2011) in FSL (Smith et al., 2004) also includes the brainstem in its parcellation, which is based on active shape and appearance models. (Introduction)
 
> [1.12] In this paper, we present a supervised segmentation method for the midbrain, pons, medulla and superior cerebellar peduncle (SCP). (Introduction)
 
`2.` Faber J, Kügler D, Bahrami E, et al. CerebNet: A fast and reliable deep-learning pipeline for detailed cerebellum sub-segmentation. *Neuroimage.* 2022;264. doi:10.1016/j.neuroimage.2022.119703 · [Access through your institution](https://libkey.io/libraries/519/articles/537342226/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [2.1] The protocols for segmenting the cerebellum on MR images differ in the level of detail at which single anterior lobules and vermal subsegments are distinguished or aggregated ( Bogovic et al., 2013b; Diedrichsen, 2006; Park et al., 2014 ). (Introduction)
 
> [2.2] To allow deeper analysis of the GM/WM boundary in the cerebellum, we extend our protocol by a fine-grained segmentation of CWM strands projecting into the cerebellar cortex. (Introduction)
 
> [2.3] This work introduces CerebNet, a fully automated, extensively validated, deep learning method for the lobular segmentation of the cerebellum, including the separation of gray and white matter. (Abstract)
 
> [2.4] Since SUIT does not provide segmentations for CWM, we additionally process all 160 images with FreeSurfer (FS) ( Fischl et al., 2002 ) and merge FSgenerated CWM with the cerebellar sub-regions labels from SUIT. (Methods)
 
`3.` Casamitjana A, Mancini M, Robinson ED, et al. A probabilistic histological atlas of the human brain for MRI segmentation. *Nature.* 2025;648:678 - 685. doi:10.1038/s41586-025-09708-2 · [Access through your institution](https://libkey.io/libraries/519/articles/687708211/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [3.1] Although these scores cannot assess segmentation accuracy at the subregion level, they do enable evaluation on a much larger multisite cohort, as well as comparison with the Allen MNI template—the only competing histological (or rather, histology-inspired) atlas that can segment the whole brain in vivo. (Results)
 
> [3.2] We developed artificial intelligence-enabled methods to align roughly 10,000 histological sections from five whole brain hemispheres into three-dimensional volumes and to produce delineations for 333 ROIs on these sections. (Abstract)
 
> [3.3] We also created a companion Bayesian tool for automatic segmentation of these ROIs in magnetic resonance imaging (MRI) scans. (Abstract)
 
> [3.4] Figure 3 shows close-ups of orthogonal slices of the atlas, which model voxel-wide probabilities for the 333 ROIs on a 0.2-mm isotropic grid. (Results)
 
> [3.5] For this reason, it only has 138 ROIs—whereas NextBrain has 333. (Results)
 
> [3.6] Using FreeSurfer 7.0 as a silver standard, we computed Dice scores for our segmentations at the level of whole regions: that is, the level of granularity provided by FreeSurfer. (Results)
 
`4.` Bosco P, Giuliano A, Delafield-Butt J, Muratori F, Calderoni S, Retico A. Brainstem enlargement in preschool children with autism: Results from an intermethod agreement study of segmentation algorithms. *Human Brain Mapping.* 2018;40:19 - 7. doi:10.1002/hbm.24351 · [Access through your institution](https://libkey.io/libraries/519/articles/230069069/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [4.1] The differences (see figure 3) are particularly striking with FSL-FIRST method, where the estimate of the volume is in general lower, except for the systematic higher estimate of the brainstem structure in the posterior part corresponding to the fourth ventricle which is included in FSL-FIRST brainstem segmentation protocol. (Results)
 
> [4.2] Both Pearson correlations and Dice indices, showed poor agreement between FSL-FIRST and the other methods (ANTs and FS), which by contrast, yielded Pearson correlations greater than 0.93 and average Dice indices greater than 0.76 when compared with each other. (Abstract)
 
`5.` Nigro S, Cerasa A, Zito G, et al. Fully Automated Segmentation of the Pons and Midbrain Using Human T1 MR Brain Images. *PLoS ONE.* 2014;9. doi:10.1371/journal.pone.0085618 · [Access through your institution](https://libkey.io/libraries/519/articles/49221515/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [5.1] In the neuroradiological and neuroimaging community there is no consensus on which anatomical landmark might be useful to improve this separation. (Methods)
 
> [5.2] To do that, the lateral margins of the superior cerebellar peduncle were used for demarcating the boundaries between the middle cerebellar peduncles and pons. (Methods)
 
`6.` Rowe O, Rangaprakash D, Eichler FS, Schmahmann J, Barry RL, Stephen CD. Brainstem Substructure Atrophy in Late-Onset GM2-Gangliosidosis Imaging Using Automated Segmentation. *Cerebellum (London, England).* 2025;24:50 - 50. doi:10.1007/s12311-025-01803-4 · [Access through your institution](https://libkey.io/libraries/519/articles/652208042/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [6.1] There were differences between LOGG and controls in the pons (12,785.06 ± 1,603.84 vs. 15,457.14 ± 2,748.41 mm3, p = 0.0069) and SCP (196.93 ± 31.20 vs. 293.57 ± .70.16 mm3, p = 0.0003). (Abstract)
 
> [6.2] Ten genetically defined LOGG patients (7 LOTS, 3 LOSD) and 7 age-matched healthy controls were consented. (Methods)
 
`7.` Sander L, Pezold S, Andermatt S, et al. Accurate, rapid and reliable, fully automated MRI brainstem segmentation for application in multiple sclerosis and neurodegenerative diseases. *Human Brain Mapping.* 2019;40:4091 - 4104. doi:10.1002/hbm.24687 · [Access through your institution](https://libkey.io/libraries/519/articles/325751005/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [7.1] Accuracy | |-------------------|---------------------|--------------|---------------------------------|---------------------------------|---------------------|----------| | | Mean volume [mm 3 ] | SD | Mean volume [mm 3 ] | SD | Mean volume [mm 3 ] | SD | | Brainstem | 24,246.86 | 3,041.96 | 24,250.25 | 3,046.03 | 24,418.10 | 3,164.43 | | Mesencephalon | 5,803.48 | 603.40 | 5,799.21 | 609.89 | 5,910.43 | 613.95 | | Pons | 14,005.89 | 2,084.54 | 13,994.06 | 2,085.17 | 14,081.53 | 2,213.21 | | Medulla oblongata | 4,437.50 | 541.35 | 4,456.98 | 531.86 | 4,426.15 | 536.38 | Note: Mean brainstem and brainstem substructure volumes, obtained by the different segmentation methods with FreeSurfer (Iglesias et al., 2015), manually edited segmentations and MD-GRU. (Results)
 
`8.` He M, Kis-Jakab G, Komáromy H, et al. Volumetric alteration of brainstem in female migraineurs with and without aura.. *Clinical neurology and neurosurgery.* 2023;236:108089. doi:10.1016/j.clineuro.2023.108089 · [Access through your institution](https://libkey.io/libraries/519/articles/598708932/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [8.1] Results: Migraineurs had greater medulla volume (MwoA 3552 ± 459 mm3, MwA 3424 ± 448 mm3) than healthy controls (3236 ± 411 mm3). (Abstract)
 
`9.` Elvsåshagen T, Bahrami S, Van Der Meer D, et al. The genetic architecture of human brainstem structures and their involvement in common brain disorders. *Nature Communications.* 2019;11. doi:10.1038/s41467-020-17376-1 · [Access through your institution](https://libkey.io/libraries/519/articles/404250311/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [9.1] The whole brainstem, midbrain, pons, SCP, and medulla oblongata were then delineated using Freesurfer 6.0 22 and Bayesian brainstem segmentation 20 . (Methods)
 
> [9.2] We then manually assessed the delineations in all MRI data sets ( n = 57,298) by visually inspecting twelve sagittal view fi gures of the segmentations for each participant, as shown in Supplementary Fig. 1. (Methods)
 
`10.` Xue H, Fu J, Geng Z, et al. The genetic architecture of brainstem structures. *Nature Communications.* 2025;17. doi:10.1038/s41467-025-67221-6 · [Access through your institution](https://libkey.io/libraries/519/articles/691534142/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [10.1] The total brainstem volume (including SCP) was included as an additional covariate to perform GWAS for the relative volumes of the medulla, pons, midbrain, and SCP 12 . (Methods)
 
> [10.2] As the prior GWAS 14 applied the automatic subcortical segmentation method to obtain the whole brainstem (excluding SCP) volume, we also employed the same approach to calculate the whole brainstem volume for each participant from CHIMGEN, ABCD, or UKB, which was used in GWAS for the trait. (Methods)
 
`11.` Bocchetta M, Cardoso JM, Cash D, Ourselin S, Warren J, Rohrer J. Patterns of regional cerebellar atrophy in genetic frontotemporal dementia. *NeuroImage : Clinical.* 2016;11:287 - 290. doi:10.1016/j.nicl.2016.02.008 · [Access through your institution](https://libkey.io/libraries/519/articles/232589166/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [11.1] Controls (18) C9orf72 (17) MAPT (20) GRN (7) C9orf72 vs Controls MAPT vs Controls GRN vs Controls C9orf72 MAPT GRN GRN vs C9orf72 GRN vs MAPT C9orf72 vs MAPT Structure Mean SD Mean SD Mean SD Mean SD p-Value p-Value p-Value % difference vs controls % difference vs controls % difference vs controls p-Value p-Value p-Value Lobule I–IV 6236 670 6363 861 6250 433 7053 674 0.660 0.534 0.012 −2.0 −0.2 −13.1 0.099 0.013 0.775 Lobule V 8911 897 8623 905 8922 742 9371 737 0.503 0.553 0.125 3.2 −0.1 −5.2 0.034 0.145 0.311 Lobule VI 15639 1452 14587 2000 15277 1390 16590 1930 0.067 0.534 0.244 6.7 2.3 −6.1 0.034 0.104 0.133 Lobule VIIa-Crus I 21091 2467 17885 1984 20360 3303 19298 2237 b0.0005 0.675 0.125 15.2 3.5 8.5 0.099 0.498 0.013 Lobule VIIa-Crus II 17610 2311 15807 2863 16015 1822 17666 752 0.067 0.033 0.534 10.2 9.1 −0.3 0.147 0.036 0.845 Lobule VIIb 7565 802 7330 976 7488 781 7460 611 0.483 0.515 0.534 3.1 1.0 1.4 0.951 0.850 0.729 Lobule VIIIa 7799 704 7724 927 7957 558 8033 756 0.909 0.361 0.495 1.0 −2.0 −3.0 0.534 0.766 0.517 Lobule VIIIb 7451 900 6926 894 7318 1077 7233 834 0.245 0.806 0.790 7.0 1.8 2.9 0.288 0.685 0.270 Lobule IX 6200 828 5563 799 5520 709 5659 864 0.013 0.010 0.158 10.3 11.0 8.7 0.318 0.533 0.988 Lobule X 860 88 850 108 909 124 882 146 0.883 0.228 1.000 1.1 −5.7 −2.6 0.757 0.607 0.149 Vermis 5034 375 4624 584 4554 412 4862 408 0.025 0.001 0.297 8.1 9.5 3.4 0.260 0.145 0.940 Dentate nuclei 3291 426 3303 406 3214 422 3325 524 1.000 0.361 0.836 −0.4 2.3 −1.0 0.951 0.766 0.341 Interposed nuclei 392 44 351 31 357 31 394 54 0.004 0.007 0.836 10.4 8.9 −0.4 0.099 0.104 0.821 Fastigial nuclei 57 10 51 6 51 6 56 10 0.062 0.033 0.836 10.5 11.1 0.9 0.147 0.104 0.988 288 M. Bocchetta et al. (Results)
 
> [11.2] We then used an atlas propagation and label fusion strategy of the Diedrichsen cerebellar atlas to automatically extract subregions of the cerebellum: the cerebellar lobules (I-IV, V, VI, VIIa-Crus I, VIIa-Crus II, VIIb, VIIIa, VIIIb, IX and X), the vermis and the deep nuclei (dentate, in- terposed and fastigial) (Cardoso et al., 2015; Diedrichsen et al., 2009, 2011). (Methods)
 
> [11.3] All vol- umes were corrected for TIV, which was calculated using the Statistical Parametric Mapping (SPM) 12 software, version 6470 (www.fil.ion.ucl. ac.uk/spm), running under Matlab R2014b (Math Works, Natick, MA, USA) (Malone et al., 2015). (Methods)
 
`12.` Tiemeier H, Lenroot R, Greenstein D, Tran L, Pierson R, Giedd J. Cerebellum development during childhood and adolescence: a longitudinal morphometric MRI study. *NeuroImage.* 2009;49:63 - 70. doi:10.1016/j.neuroimage.2009.08.016 · [Access through your institution](https://libkey.io/libraries/519/articles/44223868/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [12.1] As there is not a clearly separated vermis within the anterior lobe, this term is applied to the midline and paramedian sectors of the anterior lobe as an extension of the structures visible in the posterior and inferior aspects of the cerebellum (Schmahmann, 2000). (Methods)
 
> [12.2] In this study we apply an innovative image analysis technique to quantify total cerebellar volume and 11 subdivisions (i.e. anterior, superior posterior, and inferior posterior lobes, corpus medullare, and three vermal regions) from anatomic brain MRI scans from 25 healthy females and 25 healthy males aged 5-24 years, each of whom was scanned at least three times at approximately 2-year intervals. (Abstract)
 
`13.` Diedrichsen J, Balsters J, Flavell J, Cussans E, Ramnani N. A probabilistic MR atlas of the human cerebellum. *NeuroImage.* 2009;46 1:39-46. doi:10.1016/j.neuroimage.2009.01.045 · [Access through your institution](https://libkey.io/libraries/519/articles/57296644/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [13.1] Furthermore, maximum probability maps derived from these atlases can be used to define regions of interest (ROIs) in functional neuroimaging and neuroanatomical research. (Abstract)
 
> [13.2] We separately masked the lobules on T1-weighted MRI scans (1 mm isotropic resolution) of 20 healthy young participants (10 male, 10 female, average age 23.7 yrs). (Abstract)
 
> [13.3] These cerebella were then aligned to the standard or non-linear version of the whole-brain MNI152 template using a number of commonly used normalization algorithms, or to a previously published cerebellum-only template (Diedrichsen, J., 2006. A spatially unbiased atlas template of the human cerebellum. NeuroImage 33, 127-138.). (Abstract)
 
`14.` Hodgdon EA, Anderson R, Al-Azzawi HM, et al. MRI morphometry of the anterior and posterior cerebellar vermis and its relationship to sensorimotor and cognitive functions in children. *Developmental Cognitive Neuroscience.* 2024;67. doi:10.1016/j.dcn.2024.101385 · [Access through your institution](https://libkey.io/libraries/519/articles/614019523/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [14.1] Rigorous morphometric measurement of the cerebellar vermis lobules I-V and VI-VII was performed using a manual segmentation protocol assisted by computerized formula of FreeSurfer (Fischl, Sereno, Tootell, et al., 1999). (Methods)
 
> [14.2] The obtained cerebellar volume included gray and white matter. (Methods)
 
`15.` Bernard JA, Seidler R, Hassevoort KM, et al. Resting state cortico-cerebellar functional connectivity networks: a comparison of anatomical and self-organizing map approaches. *Frontiers in Neuroanatomy.* 2012;6. doi:10.3389/fnana.2012.00031 · [Access through your institution](https://libkey.io/libraries/519/articles/51384714/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [15.1] Timecourses were extracted from the lobules of the right hemisphere and vermis. (Abstract)
 
`16.` Casamitjana A, Mancini M, Robinson ED, et al. A next-generation, histological atlas of the human brain and its application to automated brain MRI segmentation. *bioRxiv.* 2024;648:678 - 685. doi:10.1101/2024.02.05.579016
 
> [16.1] | ctx-parsorbitalis | 0.43236 | 0.42366 | | Pons-nuc | 0.65229 | 0.5514 | | ctx-rh-posteriorcingulate | 0.69731 | 0.56382 | | ctx-rostralanteriorcingulate | 0.67935 | 0.41073 | | ctx-pericalcarine | 0.5253 | 0.38703 | | ctx-entorhinal | 0.72093 | 0.65273 | | ctx-temporalpole | 0.54659 | 0.49421 | | ctx-parahippocampal | 0.73123 | 0.66232 | | Left-PuL | 0.76937 | 0.78917 | | Left-external-pallidum | 0.78074 | 0.79407 | | Left-MDl | 0.87316 | 0.86494 | | ctx-frontalpole | 0.11589 | 0.3493 | | ctx-transversetemporal | 0.63184 | 0.56123 | | Left-VA | 0.6994 | 0.73112 | | CA1 | 0.78967 | 0.60729 | | Fornix | 0.30751 | 0.29593 | | Claustrum | 0.48364 | 0.3618 | (Body)
 
`17.` Puonti O, Nolan J, Dicamillo R, et al. Fast segmentation with the NextBrain histological atlas. *bioRxiv.* 2025;4. doi:10.1101/2025.09.22.673638
 
> [17.1] Because these atlases were constructed from manual delineations based on anatomical protocols that differ from those of NextBrain, we mapped each label in the atlases to the closest corresponding NextBrain ROI or set of ROIs. (Results)
 
> [17.2] Instead, we assess performance at the whole-ROI level by merging the automatically segmented subregions to approximate the FreeSurfer protocol (Fischl et al., 2002). (Methods)
 
> [17.3] As no reference segmentations are available for 20 out of the 21 subjects, we instead report the volumes of the substructures of the hippocampus, thalamus and amygdala in tables in Tables 1, 2, and 3, and compare them with those from other existing ex vivo atlases (Iglesias et al., 2015, 2018; Saygin et al., 2017). (Results)
 
`18.` Hammers A, Allom R, Koepp M, et al. Three‐dimensional maximum probability atlas of the human brain, with particular reference to the temporal lobe. *Human Brain Mapping.* 2003;19. doi:10.1002/hbm.10123 · [Access through your institution](https://libkey.io/libraries/519/articles/26141009/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [18.1] The total regional volume was 46% bigger, which, assuming equal differences in each direction, corresponds to an average 13.5% increase in each dimension (Table I). (Discussion)
 
`19.` Giff A, Noren G, Magnotti J, et al. Spatial Normalization Discrepancies between Native and MNI152 Brain Template Scans in Gamma Ventral Capsulotomy Patients. *Psychiatry research. Neuroimaging.* 2023;329:111595 - 111595. doi:10.1016/j.pscychresns.2023.111595 · [Access through your institution](https://libkey.io/libraries/519/articles/543822548/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [19.1] MNI measurements were significantly larger than native measurements across several structures in both coronal and axial slices, and we found that MNI transformation increases the size of many subcortical structures in a significant and proportional way for both females and males. (Abstract)
 
`20.` Diedrichsen J. A spatially unbiased atlas template of the human cerebellum. *NeuroImage.* 2006;33 1:127-38. doi:10.1016/j.neuroimage.2006.05.056 · [Access through your institution](https://libkey.io/libraries/519/articles/57081272/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [20.1] Compared to normalization to the MNI whole-brain template, the new method strongly improves the alignment of individual fissures, reducing their spatial spread by 60%, and improves the overlap of the deep cerebellar nuclei. (Abstract)
 
`21.` Ota K, Nakazato Y, Oyama G. Developing and validating an elderly brain template: A comprehensive comparison with MNI152 for age-specific neuroimaging analyses. *NeuroImage.* 2025. doi:10.1016/j.neuroimage.2025.121473 · [Access through your institution](https://libkey.io/libraries/519/articles/675969805/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [21.1] The largest interaction effects (≈1-4 % gain) occurred in the caudate, thalamus, hippocampus, and amygdala, whereas the brainstem, pallidum, and putamen showed minimal differences between templates. (Abstract)
 
`22.` Luft AR, Skalej M, Schulz J, et al. Patterns of age-related shrinkage in cerebellum and brainstem observed in vivo using three-dimensional MRI volumetry.. *Cerebral cortex.* 1999;9 7:712-21. doi:10.1093/cercor/9.7.712 · [Access through your institution](https://libkey.io/libraries/519/articles/27074734/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [22.1] The brainstem was divided into midbrain, pons/tegmentum pontis/ cerebellar peduncles and medulla. (Methods)
 
`23.` Lambert C, Lutti A, Helms G, Frackowiak RSJ, Ashburner J. Multiparametric brainstem segmentation using a modified multivariate mixture of Gaussians. *NeuroImage : Clinical.* 2013;2:684 - 694. doi:10.1016/j.nicl.2013.04.017 · [Access through your institution](https://libkey.io/libraries/519/articles/203500722/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [23.1] The following anatomical landmarks were used: Inferiorly, the level of the foramen magnum; posteriorly, a vertical line following the most posterior aspect of the medulla, which included part of the middle and superior cerebellar peduncles. (Methods)
 
`24.` Magnusson MI, Love A, Ellingsen LM. Automated brainstem parcellation using multi-atlas segmentation and deep neural network. 2021;11596:115962E - 115962E-6. doi:10.1117/12.2581129
 
> [24.1] We did not compute the DSC for FreeSurfer based on our manual masks, since the two methods are based on slightly different delineation protocols. (Results)
 
`25.` Metwally MI, Basha M, AbdelHamid GA, et al. Neuroanatomical MRI study: reference values for the measurements of brainstem, cerebellar vermis, and peduncles.. *The British journal of radiology.* 2021. doi:10.1259/bjr.20201353 · [Access through your institution](https://libkey.io/libraries/519/articles/441236824/full-text-file?utm_source=api_3508&allow_speedbump=true)
 
> [25.1] On midsagittal T 1 weighted and axial T 2 weighted images specific linear diameters and ratios of the brainstem, cerebellar vermis, and peduncles were attained. (Methods)
