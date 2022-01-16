#### dmrscaler_real_data.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.$TASK_ID
#$ -j y
## Edit the line below as needed:
#$ -l highp,h_rt=2:00:00,h_data=6G
## Modify the parallel environment
## and the number of cores as needed:
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
# Job array indexes

NUM_DATA_SETS=`wc -l < data_set_table.csv`
NUM_METHOD_SETS=`wc -l < method_table.csv`
UPPER_LIM=$(expr $NUM_SIMUL_SETS \* $NUM_METHOD_SETS )
#$ -t 1-${UPPER_LIM}:1


# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "


# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load R/4.1.0

DATA_SET_ID=$(expr 1+ $(expr $SGE_TASK_ID % $NUM_DATA_SETS ))
METHOD_SET_ID=$(expr 1 + $(expr $SGE_TASK_ID / $NUM_DATA_SETS ))


Rscript real_data_individual_run.R $SIMUL_SET_ID $METHOD_SET_ID


# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "

#### dmrscaler_real_data.sh STOP ####
