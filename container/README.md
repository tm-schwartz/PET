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
Still missing (not reproducible from public data): the hand-drawn
`ponsvermis_1mm_bin-mask.nii.gz` / `Composite-mask.nii.gz` — obtain from
LONI IDA "UC Berkeley – FDG MetaROI NPDKA Analysis" or reconstruct.

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
