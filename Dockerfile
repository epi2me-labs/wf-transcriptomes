ARG BASEIMAGE=epi2melabs/base-workflow-image:latest
FROM $BASEIMAGE

# Minimal install for example purposes 
RUN \
    . $CONDA_DIR/etc/profile.d/mamba.sh \
    && micromamba activate \ 
    && micromamba install -y \
            pysam \
        -c anaconda -c conda-forge -c bioconda -q -y \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME

USER $WF_UID
WORKDIR $HOME
