#!/bin/bash
#
# engaging_launch_Rstudio.sh
#
# Launch an RStudio Server session in the project's Apptainer container on the
# MIT ORCD / Engaging cluster.
#
#   Usage:
#       sbatch engaging_launch_Rstudio.sh
#       squeue -u $USER                      # wait for state = R
#       cat rstudio.<node>.<jobid>.out       # tunnel command + password
#
# The container is a pre-built, read-only, maintainer-managed SIF. This script
# deliberately does NOT build or pull it:
#   * Pulling docker:// at job start runs mksquashfs inside the job and gets
#     OOM-killed on default memory allocations.
#   * A fixed SIF guarantees every figure in the manuscript came from the same
#     software stack.
# If the image is missing or unreadable, contact the BCC (see README) -- do not
# build your own into the shared directory.
#
# Repo: KochInstitute-Bioinformatics/bh4-manuscript-publicData-analysis
# ---------------------------------------------------------------------------

#SBATCH --job-name=rstudio
#SBATCH --output=rstudio.%N.%j.out
#SBATCH --partition=mit_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-user=charliew@mit.edu # edit this with your own email address
#SBATCH --mail-type=BEGIN,END,FAIL

set -uo pipefail

# ===========================================================================
# CONFIG -- override by exporting before sbatch, e.g.
#   RSTUDIO_SIF=/path/to/alt.sif sbatch engaging_launch_Rstudio.sh
# ===========================================================================
DSMI=${DSMI:-/orcd/data/ki/003/core/bcc/DSMI_Resources}

# Shared read-only container image (managed by the BCC; do not modify)
RSTUDIO_SIF=${RSTUDIO_SIF:-$DSMI/containers/bulk_r451_v1.sif}

# Persistent container HOME: R library, .Rprofile, RStudio prefs. Per-user.
RSTUDIO_HOME=${RSTUDIO_HOME:-$DSMI/rstudio_home/$USER}

# Per-job scratch (/run, /tmp, rserver sqlite DB). Real disk, not tmpfs.
SCRATCH_BASE=${SCRATCH_BASE:-$DSMI/tmp/$USER}

# Login node used in the printed tunnel instructions
LOGIN_NODE=${LOGIN_NODE:-orcd-login.mit.edu}
# ===========================================================================

module load apptainer
module load miniforge          # python used for the free-port pick

# --- preflight: the image must already exist ---------------------------------
if [[ ! -r "$RSTUDIO_SIF" ]]; then
  cat 1>&2 <<ERR

======================================================================
ERROR: cannot read the analysis container image

    $RSTUDIO_SIF

This image is built and maintained centrally by the KI Bioinformatics
Core. Please do NOT build your own copy into the shared directory --
contact the BCC (see README) and we will fix or restore it.

Diagnostics to include in your message:

    namei -l $RSTUDIO_SIF
    groups
    id -un

Most common cause: your account is not in the BCC storage group
(orcd_rg_hstor003_pg_ki_bcc), so you cannot traverse to the image.
======================================================================

ERR
  exit 1
fi

NCPUS=${SLURM_CPUS_PER_TASK:-1}   # SLURM_JOB_CPUS_PER_NODE can be "8(x1)"

mkdir -p "$RSTUDIO_HOME" "$SCRATCH_BASE" || {
  echo "ERROR: cannot create $RSTUDIO_HOME / $SCRATCH_BASE" 1>&2; exit 1; }

# --- per-job working directory ----------------------------------------------
workdir=$(mktemp -d "${SCRATCH_BASE}/rstudio-${SLURM_JOB_ID}.XXXXXX") || exit 1
cleanup() { rm -rf "$workdir"; }
trap cleanup EXIT INT TERM

mkdir -p -m 700 "${workdir}/run" \
                "${workdir}/tmp" \
                "${workdir}/var/lib/rstudio-server"

cat > "${workdir}/database.conf" <<'CONF'
provider=sqlite
directory=/var/lib/rstudio-server
CONF

cat > "${workdir}/rsession.sh" <<RSESSION
#!/bin/sh
# Confine BLAS/OpenMP to the cores Slurm actually allocated.
export OMP_NUM_THREADS=${NCPUS}
export OPENBLAS_NUM_THREADS=${NCPUS}
export MKL_NUM_THREADS=${NCPUS}
export R_PARALLELLY_AVAILABLECORES_FALLBACK=${NCPUS}
export TMPDIR=/tmp
exec /usr/lib/rstudio-server/bin/rsession "\${@}"
RSESSION
chmod +x "${workdir}/rsession.sh"

# --- bind mounts -------------------------------------------------------------
binds="${workdir}/run:/run"
binds+=",${workdir}/tmp:/tmp"
binds+=",${workdir}/database.conf:/etc/rstudio/database.conf"
binds+=",${workdir}/rsession.sh:/etc/rstudio/rsession.sh"
binds+=",${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server"
# Whole ORCD tree at real paths so project data resolves as it does on the host
binds+=",/orcd:/orcd"
# Convenience aliases used by scripts in this repo
binds+=",/orcd/data/ki/003/core/bcc/IGB_Resources/annotation_files:/annotationFiles"
binds+=",/orcd/data/ki/003/core/bcc/IGB_Resources/scripts:/scripts"
binds+=",/orcd/data/ki/003/core/bcc/Genomes:/Genomes"
binds+=",${DSMI}:/DSMI_Resources"
export APPTAINER_BIND="$binds"

# --- apptainer / session environment ----------------------------------------
# Keep any incidental apptainer scratch inside the job dir, never $HOME.
export APPTAINER_CACHEDIR="${workdir}/apptainer-cache"
export APPTAINER_TMPDIR="${workdir}/tmp"

export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0
export APPTAINERENV_USER=$(id -un)
export APPTAINERENV_PASSWORD=$(openssl rand -base64 18 | tr -dc 'A-Za-z0-9' | head -c 16)

PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
readonly PORT

# --- connection instructions --------------------------------------------------
NODE=$(hostname -s)
INFO="rstudio-${SLURM_JOB_ID}.connect.txt"

cat > "$INFO" <<INFOEOF
==============================================================================
 RStudio Server   job ${SLURM_JOB_ID}   node ${NODE}   port ${PORT}
 image : ${RSTUDIO_SIF}
 cores : ${NCPUS}      mem: ${SLURM_MEM_PER_NODE:-?} MB
==============================================================================

OPTION A -- VS Code Remote
  Ports panel > Forward a Port > enter:   ${NODE}:${PORT}
  Then open:  http://localhost:${PORT}

OPTION B -- SSH tunnel from your laptop (one hop; try this first)
  ssh -N -L ${PORT}:${NODE}:${PORT} ${USER}@${LOGIN_NODE}
  Then open:  http://localhost:${PORT}

OPTION C -- two-hop tunnel (if the login node cannot route to ${NODE})
  ssh -t -L ${PORT}:localhost:${PORT} ${USER}@${LOGIN_NODE} \\
      ssh -t ${NODE} -L ${PORT}:localhost:${PORT}
  Then open:  http://localhost:${PORT}

LOGIN
  user:     ${APPTAINERENV_USER}
  password: ${APPTAINERENV_PASSWORD}

NOTES
  * Container HOME (persistent R library / prefs): ${RSTUDIO_HOME}
  * Project data is visible at its real /orcd/... paths.
  * Job scratch (auto-deleted on exit): ${workdir}
  * Image provenance: see ${RSTUDIO_SIF}.provenance.txt

WHEN FINISHED
  1. In RStudio click the power button (top right) to quit the session.
  2. On a login node:   scancel -f ${SLURM_JOB_ID}
==============================================================================
INFOEOF

chmod 600 "$INFO"
cat "$INFO" 1>&2

# --- record what we ran, for reproducibility ---------------------------------
{
  echo "--- provenance ---"
  echo "sif        : $RSTUDIO_SIF"
  echo "sif mtime  : $(stat -c %y "$RSTUDIO_SIF" 2>/dev/null)"
  echo "sif sha256 : $(sha256sum "$RSTUDIO_SIF" 2>/dev/null | cut -d' ' -f1)"
  echo "apptainer  : $(apptainer --version)"
  echo "node       : $NODE"
  echo "git commit : $(git -C "$SLURM_SUBMIT_DIR" rev-parse --short HEAD 2>/dev/null || echo 'n/a')"
} 1>&2

# --- launch -------------------------------------------------------------------
apptainer exec --cleanenv \
  -H "${RSTUDIO_HOME}:/home/rstudio" \
  "$RSTUDIO_SIF" \
  /usr/lib/rstudio-server/bin/rserver \
    --server-user "${USER}" \
    --www-port "${PORT}" \
    --auth-none=0 \
    --auth-pam-helper-path=pam-helper \
    --auth-stay-signed-in-days=30 \
    --auth-timeout-minutes=0 \
    --rsession-path=/etc/rstudio/rsession.sh

rc=$?
printf 'rserver exited with status %d\n' "$rc" 1>&2
rm -f "$INFO"
exit $rc