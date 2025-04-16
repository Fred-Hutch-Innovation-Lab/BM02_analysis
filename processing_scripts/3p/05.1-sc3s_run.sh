#!/bin/bash
# ml Apptainer
# ml JupyterLab/4.0.3-GCCcore-12.2.0
apptainer exec \
    --bind /fh/fast/_IRC/FHIL/grp/BM_paper/analysis/:/home/dgratz/workspace \
    --bind /app/software/JupyterLab/4.0.3-GCCcore-12.2.0/share/jupyter/lab \
    /fh/fast/_IRC/FHIL/grp/BM_paper/analysis/config/sc3s_jupyterlab.sif \
    jupyter lab --ip=$(hostname) --port=$(fhfreeport) --no-browser --allow-root