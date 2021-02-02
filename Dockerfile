ARG BASEIMAGE=epi2melabs/base-workflow-image:latest
FROM $BASEIMAGE

# Minimal install for example purposes
COPY environment.yaml $HOME/environment.yaml 
RUN \
    . $CONDA_DIR/etc/profile.d/mamba.sh \
    && micromamba activate \
    && micromamba install --file $HOME/environment.yaml \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME

USER $WF_UID
WORKDIR $HOME
