#!/bin/bash
# Script to launch flair-FLUKA jobs on a cluster using OpenPBS
# Author: Vasilis.Vlachoudis@cern.ch
# Date: Nov-2006

NAME="flair${!#}"
JOBFILE=${NAME}.job
NAME=`echo ${NAME} | cut -c 1-15`

cat > ${JOBFILE} << EOF
#!/bin/bash
#PBS -q normal
#PBS -N ${NAME}
export FLUKA=$FLUKA
export FLUPRO=$FLUPRO
cd $PWD
$*
EOF

echo "Submitting job ${JOBFILE}"
chmod +x ${JOBFILE}
qsub ${JOBFILE}
