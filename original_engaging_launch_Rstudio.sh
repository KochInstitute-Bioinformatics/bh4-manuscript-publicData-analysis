#!/bin/bash

#SBATCH -n 8
#SBATCH --job-name=Rstudio       # Assign an short name to your job
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --partition=mit_normal
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-user=charliew@mit.edu # The user to email
#SBATCH --mail-type=ALL

module load apptainer
module load miniforge

workdir=$(python -c 'import tempfile; print(tempfile.mkdtemp())')

mkdir -p -m 700 ${workdir}/run ${workdir}/tmp ${workdir}/var/lib/rstudio-server
cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

cat > ${workdir}/rsession.sh <<END
#!/bin/sh
export OMP_NUM_THREADS=${SLURM_JOB_CPUS_PER_NODE}
exec /usr/lib/rstudio-server/bin/rsession "\${@}"
END

chmod +x ${workdir}/rsession.sh

export APPTAINER_BIND="${workdir}/run:/run,${workdir}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server"
export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0
export APPTAINERENV_USER=$(id -un)
export APPTAINERENV_PASSWORD="koch76"
#export APPTAINERENV_PASSWORD=$(echo $RANDOM | base64 | head -c 20)

readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')

cat 1>&2 <<END

Connect to RStudio using either of these options:

1. Forward port in VSCode.

   Ports > Forward a Port > Set "Port" to: ${HOSTNAME}:${PORT}

   Navigate to http://localhost:${PORT} in your web browser.
   Login with the following credentials:

   user: ${APPTAINERENV_USER}
   password: ${APPTAINERENV_PASSWORD}

2. Manually forward ports with SSH.

   From your MacOS Terminal / Windows Powershell, run:
   
   ssh -N -L ${PORT}:${HOSTNAME}:${PORT} ${USER}@orcd-login.mit.edu
   
   Navigate to http://localhost:${PORT} in your web browser.
   Login with the following credentials:

   user: ${APPTAINERENV_USER}
   password: ${APPTAINERENV_PASSWORD}

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f ${SLURM_JOB_ID}
END

apptainer exec --cleanenv \
               -B /orcd:/orcd \
               docker://docker.io/bumproo/r4_5_3_singlecell_bulk_rnaseq:latest /usr/lib/rstudio-server/bin/rserver \
                  --server-user ${USER} --www-port ${PORT} \
                  --auth-none=0 \
                  --auth-pam-helper-path=pam-helper \
                  --auth-stay-signed-in-days=30 \
                  --auth-timeout-minutes=0 \
                  --rsession-path=/etc/rstudio/rsession.sh