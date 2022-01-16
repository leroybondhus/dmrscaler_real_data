#### dmrscaler_real_data_setup.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l highp,h_rt=2:00:00,h_data=4G
## Modify the parallel environment
## and the number of cores as needed:
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea



# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh

## Edit the line below as needed:
module load bedtools/2.30.0
module load anaconda3/2021.11
module load R/4.1.0

source $CONDA_DIR/etc/profile.d/conda.sh
conda activate myconda


Rscript real_data_setup.R


# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### dmrscaler_real_data_setup.sh  STOP ####
