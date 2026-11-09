# pet-fdg: FDG PET SUVR pipeline (pet_py) packaged for generalized dataset use.
#
# Data-free image by design. Bind at run time:
#   /input      ro   source dataset (BIDS-like study dirs)
#   /output     rw   derivatives + CSVs
#   /templates  ro   MNI-PET fixed images + *-mask.nii.gz ROIs (see README)
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
         curl ca-certificates bzip2 tar xz-utils git libgomp1 procps \
    && rm -rf /var/lib/apt/lists/*

# Real FSL, CPU-only variant, pinned for reproducibility.
# Provides robustfov (nipype) and libflimage (fslpy resampleToPixdims).
RUN curl -Ls https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/releases/getfsl.sh \
    | sh -s -- /opt/fsl -V 6.0.7.9 -r -n \
    && /opt/fsl/bin/robustfov -help 2>&1 | head -1

ENV FSLDIR=/opt/fsl \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLMULTIFILEQUIT=TRUE \
    LD_LIBRARY_PATH=/opt/fsl/lib

# uv-managed venv (AGENTS.md: uv only). CPU torch via uv's torch backend.
RUN curl -LsSf https://astral.sh/uv/0.7.8/install.sh | env UV_INSTALL_DIR=/usr/local/bin sh

# uv's managed interpreter must live outside /root (mode 700) or `python` in
# the venv is unreachable when the container runs with --user.
ENV UV_PYTHON_INSTALL_DIR=/opt/python

COPY pet-py /opt/src/pet-py
COPY container/make_templates.py container/batch_run.py container/build_ponsvermis.py /opt/src/container/

RUN uv venv --python 3.11 /opt/venv \
    && uv pip install --python /opt/venv/bin/python --torch-backend=cpu \
         /opt/src/pet-py "pydicom<3" "deepbet==0.0.2" dcm2niix \
    && /opt/venv/bin/python -c "import pet_py, itk, dcm2niix, fsl.data.image; print('imports ok')"

ENV PATH=/opt/venv/bin:/opt/fsl/bin:$PATH \
    VIRTUAL_ENV=/opt/venv \
    TEMPLATE_DIR=/templates \
    PARAM_FILE=/templates/param27.txt

COPY container/pet-fdg /usr/local/bin/pet-fdg

ENTRYPOINT ["pet-fdg"]
CMD ["--help"]
