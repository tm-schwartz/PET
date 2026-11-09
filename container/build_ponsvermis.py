"""Reconstruct ponsvermis_1mm_bin-mask.nii.gz on the fixed-image grid.

Stands in for the legacy FastSurfer cerebnet+pons mask (no surviving
output on storage). Provenance and citation-verified caveats:
container/README.md.

    pons   = NextBrain histological MNI152 1mm atlas labels 580/10580
             (doi:10.1038/s41586-025-09708-2)
    vermis = Diedrichsen/SUIT 2009 MNI labels 6,9,12,15,18,21,24,27
             (Vermis VI..X; doi:10.1016/j.neuroimage.2009.01.045)

Inputs (download; see README links):
    nextbrain_wholebrain_MNI152_1mm.nii.gz  from compneurobilbao/nextbrain-mni-atlas
    atl-Anatom_space-MNI_dseg.nii           from DiedrichsenLab/cerebellar_atlases
                                            (branch master)

    python build_ponsvermis.py NEXTBRAIN.nii.gz SUIT_DSEG.nii \
        FIXED.nii.gz OUT-mask.nii.gz
"""

from argparse import ArgumentParser
from pathlib import Path

import nibabel as nib
import numpy as np
from nilearn.image import load_img, resample_to_img

PONS_LABELS = [580, 10580]
VERMIS_LABELS = [6, 9, 12, 15, 18, 21, 24, 27]


def main() -> None:
    ap = ArgumentParser(__doc__)
    ap.add_argument("nextbrain")
    ap.add_argument("suit_dseg")
    ap.add_argument("fixed", help="MNI152_PET_1mm_coreg_smoothed.nii.gz")
    ap.add_argument("out")
    args = ap.parse_args()

    nb = load_img(args.nextbrain)
    pons = np.isin(np.asanyarray(nb.dataobj), PONS_LABELS)

    suit = load_img(args.suit_dseg)
    vermis_2d = np.isin(np.asanyarray(suit.dataobj), VERMIS_LABELS)
    vermis = resample_to_img(
        nib.Nifti1Image(vermis_2d.astype("uint8"), suit.affine), nb,
        interpolation="nearest",
    ).get_fdata() > 0

    union = (pons | vermis).astype("float32")
    for name, m in (("pons", pons), ("vermis", vermis)):
        c = np.argwhere(m).mean(axis=0)
        w = nb.affine[:3, :3] @ c + nb.affine[:3, 3]
        print(f"{name}: {int(m.sum())} vox, centroid {np.round(w, 1)}")
    print("overlap:", int((pons & vermis).sum()))

    out = resample_to_img(
        nib.Nifti1Image(union, nb.affine), load_img(args.fixed),
        interpolation="nearest",
    )
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_filename(args.out)
    print("wrote", args.out, int(out.get_fdata().sum()), "voxels")


if __name__ == "__main__":
    main()
