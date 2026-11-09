# pet-fdg container

Packaged `pet_py` FDG SUVR pipeline + real FSL 6.0.7.9 (CPU) + uv venv
(python 3.11, `deepbet==0.0.2`, `pydicom<3`, CPU torch, itk-elastix,
dcm2niix wheel). The image carries **no data**; everything dataset-specific
is bind-mounted.

## Build

```console
cd PET                      # repo root
docker build -f container/Dockerfile -t pet-fdg:dev .
```

## Run

```console
docker run --rm -u $(id -u):$(id -g) \
  -v $INPUT:/input:ro -v $OUTPUT:/output -v $TEMPLATES:/templates \
  pet-fdg:dev run /input --nifti_path /input/sub-<MRN>_..._PT.nii.gz ...
```

Mounts: `/input` (ro), `/output` (rw), `/templates` (rw when seeding masks).
For clusters: `apptainer build pet-fdg.sif docker-daemon:pet-fdg:dev`, run
with `apptainer run -B ... pet-fdg.sif run ...` (no `--user` needed; SIF
runs as caller).

## Templates dir

Needs the Darkner FDG-PET template (download via **browser** — the KU portal
CDN bot-blocks curl):
<https://researchprofiles.ku.dk/en/publications/fdg-pet-template-mni152-1mm/>

Then seed once:

```console
docker run --rm -v $TEMPLATES:/templates -v $DIR_WITH_NII:/seed:ro \
  pet-fdg:dev masks /seed/MNI152_PET_1mm.nii --outdir /templates
```

Produces: smoothed + stripped fixed images, all `*-mask.nii.gz` ROI masks
(MNI structural + Harvard-Oxford via fslpy atlases), `param27.txt`.

## ponsvermis reference mask

The legacy `ponsvermis_1mm_bin` mask (FastSurfer cerebnet + pons, manually
created) has no surviving output on storage, so we reconstruct it from public
atlases onto the Darkner fixed-image grid (`build_ponsvermis.py` pattern):

- **pons** = NextBrain histological MNI152 1mm atlas labels 580/10580
  (Casamitjana et al., *Nature* 2025, doi:10.1038/s41586-025-09708-2);
  11.8 cc, midline centroid (0, -24, -32).
- **vermis** = Diedrichsen/SUIT 2009 probabilistic cerebellar atlas
  (doi:10.1016/j.neuroimage.2009.01.045) deterministic MNI labels
  6,9,12,15,18,21,24,27 (Vermis VI/CrusI/CrusII/VIIb/VIIIa/VIIIb/IX/X);
  7.2 cc, centroid (0, -66, -31); zero overlap with pons (4th ventricle
  between — anatomically correct).

Citation-verified caveats (2026-11, full report in
`/tmp/opencode/ponsvermis_citation_verification.md`):

1. Pons volume is 19-25% below FreeSurfer/FIRST native-space control means
   (14.5-15.8 cc); plausible for a deterministic template-space label.
   FreeSurfer's protocol (Iglesias 2015, App. A) includes bounded MCP tissue
   via a posterior cutoff between MCP/SCP branching points, so NextBrain's
   MCP-spanning lateral extent (x ± 19 mm) is equivalent-enough — do NOT
   trim the peduncles.
2. Vermis is ~40% above the same-family in-vivo benchmark (Bocchetta 2016,
   label fusion + TIV correction: 5.0 cc); deterministic max-prob labels are
   inherently more inclusive. Our union is **posterior vermis only** — SUIT
   places lobules I-IV (largely anterior vermis) in hemispheric labels 1/2,
   which cerebnet would count as vermal. Optional refinement, not needed for
   reference-region use.

Composite (MetaROI) mask: pass Berkeley `Composite.nii` with
`masks --composite`; it is resampled nearest-neighbour onto the fixed grid
as `Composite-mask.nii.gz`.

Still missing (not reproducible from public data): the exact legacy masks —
obtain from LONI IDA "UC Berkeley – FDG MetaROI NPDKA Analysis" if the
reconstruction above must be replaced 1:1.

## Commands

- `pet-fdg info` — version report
- `pet-fdg run <data_path> [--nifti_path|--zip_file ...]` — pipeline (tests/run.py)
- `pet-fdg masks <template.nii> [--outdir] [--placeholder-ponsvermis]` — provision templates
- `pet-fdg batch <bids_root> [--pattern GLOB] [-n N]` — one series per subject
  into `/output/batch_summary.csv` (status per subject; failures never fatal)
- `pet-fdg shell` — debugging

## Notes

- `deepbet` downloads model weights on first use into `$HOME`; set
  `HOME=/tmp` (or a mounted rw path) when running `--user`.
- FSL is only used for `robustfov` and fslpy resampling; nipype sees it via
  `FSLDIR=/opt/fsl`.
