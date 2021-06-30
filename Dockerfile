ARG BASEIMAGE=ontresearch/base-workflow-image:v0.1.0
FROM $BASEIMAGE
ARG ENVFILE=environment.yaml

COPY $ENVFILE $HOME/environment.yaml
RUN \
    . $CONDA_DIR/etc/profile.d/mamba.sh \
    && micromamba activate \
    && micromamba install -n base --file $HOME/environment.yaml \
    && micromamba clean --all --yes \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME

USER $WF_UID
WORKDIR $HOME
