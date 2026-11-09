"""Provision the external templates directory for the pet-fdg pipeline.

Run inside the container where the (writable) templates dir is mounted:

    pet-fdg masks /path/to/MNI152_PET_1mm.nii.gz --outdir /templates

Inputs: Darkner's FDG-PET template (MNI152 1mm),
https://researchprofiles.ku.dk/en/publications/fdg-pet-template-mni152-1mm/
(cite the dataset when used).

Outputs in --outdir:
    MNI152_PET_1mm_raw.nii.gz                verbatim copy of the input
    MNI152_PET_1mm_coreg_smoothed.nii.gz     FWHM-smoothed fixed image (affine step)
    masked-MNI-PET-STRIPPED-SMOOTHED.nii.gz  smoothed x MNI152 brain mask (elastix step)
    *-mask.nii.gz                            MNI + Harvard-Oxford ROI masks
    param27.txt                              elastix Bspline parameters (from repo)
"""

from argparse import ArgumentParser
from os import environ
from pathlib import Path

import nibabel as nib
import numpy as np
from nilearn.image import load_img, math_img, resample_to_img, smooth_img

from pet_py.util import generate_masks

RAW = "MNI152_PET_1mm_raw.nii.gz"
SMOOTHED = "MNI152_PET_1mm_coreg_smoothed.nii.gz"
STRIPPED = "masked-MNI-PET-STRIPPED-SMOOTHED.nii.gz"
PONSVERMIS = "ponsvermis_1mm_bin-mask.nii.gz"


def placeholder_ponsvermis(outdir: Path) -> None:
    """Stand-in reference region: pons band of HO Brain-Stem + cerebellar
    midline column, on the fixed-image grid. For pipeline smoke tests only."""
    brainstem = load_img(
        next(outdir.glob("Brain-Stem_harvardoxford-subcortical-mask.nii.gz"))
    )
    cerebellum = load_img(next(outdir.glob("Cerebellum_mni-mask.nii.gz")))

    def z_indices(img, lo: float, hi: float):
        zz = np.arange(img.shape[2]) * img.affine[2, 2] + img.affine[2, 3]
        return ((zz >= lo) & (zz <= hi)).nonzero()[0]

    bs = brainstem.get_fdata(dtype="float32")
    vox = np.zeros_like(bs)
    zi = z_indices(brainstem, -35, -20)  # pons
    vox[:, :, zi] = bs[:, :, zi]

    cb = cerebellum.get_fdata(dtype="float32")
    zj = z_indices(cerebellum, -55, -15)
    xx = np.arange(cb.shape[0]) * cerebellum.affine[0, 0] + cerebellum.affine[0, 3]
    vm = np.zeros_like(cb)
    vm[:, :, zj] = cb[:, :, zj] * (np.abs(xx) < 5.0)[:, None, None]

    ref = np.clip((vox > 0) + (vm > 0), 0, 1).astype("float32")
    nib.Nifti1Image(ref, brainstem.affine).to_filename(outdir / PONSVERMIS)
    print(f"placeholder ponsvermis written ({int(ref.sum())} voxels)")


def main() -> None:
    parser = ArgumentParser(__doc__)
    parser.add_argument("template_in", help="Darkner MNI152_PET_1mm.nii(.gz)")
    parser.add_argument("--outdir", default="/templates")
    parser.add_argument("--fwhm", type=float, default=6.0, help="smoothing mm")
    parser.add_argument(
        "--composite",
        help="Berkeley MetaROIs Composite.nii (2mm MNI152Lin); resampled "
        "to the fixed-image grid as Composite-mask.nii.gz",
    )
    parser.add_argument(
        "--placeholder-ponsvermis",
        action="store_true",
        help="Build a stand-in ponsvermis_1mm_bin-mask.nii.gz from the "
        "HO Brain-Stem (pons z-band) + MNI Cerebellum (midline column) "
        "masks until the Berkeley original is available. NOT for production.",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    src = load_img(args.template_in, dtype="f4")
    nib.save(src, outdir / RAW)

    smoothed = smooth_img(src, fwhm=args.fwhm)
    smoothed.to_filename(outdir / SMOOTHED)

    fsldir = Path(environ.get("FSLDIR", "/opt/fsl"))
    brain_mask = resample_to_img(
        load_img(fsldir / "data/standard/MNI152_T1_1mm_brain_mask.nii.gz"),
        smoothed,
        interpolation="nearest",
    )
    math_img("img * m", img=smoothed, m=brain_mask).to_filename(outdir / STRIPPED)

    param = Path(__file__).parents[1].joinpath("pet-py", "resources", "param27.txt")
    (outdir / "param27.txt").write_bytes(param.read_bytes())

    generate_masks(outdir / SMOOTHED)

    if args.composite:
        comp = load_img(args.composite)
        resample_to_img(comp, load_img(outdir / SMOOTHED), interpolation="nearest"
                        ).to_filename(outdir / "Composite-mask.nii.gz")
        print("Composite-mask.nii.gz resampled onto fixed grid")

    if args.placeholder_ponsvermis:
        placeholder_ponsvermis(outdir)

    n_masks = len(list(outdir.glob("*-mask.nii.gz")))
    print(f"templates ready in {outdir}: {n_masks} ROI masks")


if __name__ == "__main__":
    main()
