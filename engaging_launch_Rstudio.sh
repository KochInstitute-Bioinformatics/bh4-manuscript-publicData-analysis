#!/bin/bash

#SBATCH -n 1
#SBATCH --job-name=Rstudio           # Assign an short name to your job
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --mail-user=example@mit.edu  # The user to email
#SBATCH --mail-type=ALL

module load apptainer
module load miniforge

workdir=$(python -c 'import tempfile; print(tempfile.mkdtemp())')
echo "workdir: " $workdir

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

binds="${workdir}/run:/run"
binds+=",${workdir}/tmp:/tmp"
binds+=",${workdir}/database.conf:/etc/rstudio/database.conf"
binds+=",${workdir}/rsession.sh:/etc/rstudio/rsession.sh"
binds+=",${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server"
binds+=",/orcd/data/ki/003/core/bcc/IGB_Resources/annotation_files:/annotationFiles"
binds+=",/orcd/data/ki/003/core/bcc/IGB_Resources/scripts:/scripts"
binds+=",/orcd/data/ki/003/core/bcc/Genomes:/Genomes"

export APPTAINER_BIND="$binds"

export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0
export APPTAINERENV_USER=$(id -un)
export APPTAINERENV_PASSWORD="koch76"
#export APPTAINERENV_PASSWORD=$(echo $RANDOM | base64 | head -c 20)

readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')

cat 1>&2 <<END

1. SSH tunnel from your workstation using the following command:

   SOCKET:${PORT}

   From local host run:
    
    ssh -t -L PORT:localhost:SOCKET kerberosID@orcd-login004.mit.edu ssh -t node -L SOCKET:localhost:SOCKET
  
2. Once connected, open your web browser and go to:

      http://localhost:PORT

3. log in to RStudio Server using the following credentials:

   user: ${APPTAINERENV_USER}
   password: ${APPTAINERENV_PASSWORD}

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f ${SLURM_JOB_ID}
END

apptainer exec --cleanenv -H $PWD:/home/rstudio docker://bumproo/bulk_r451:v1 /usr/lib/rstudio-server/bin/rserver \
            --server-user ${USER} --www-port ${PORT} \
            --auth-none=0 \
            --auth-pam-helper-path=pam-helper \
            --auth-stay-signed-in-days=30 \
            --auth-timeout-minutes=0 \
            --rsession-path=/etc/rstudio/rsession.sh 
printf 'rserver exited' 1>&2