#!/bin/bash
path="codes/main_scripts/main_gam_designs.R"
log_path="log_files/designs/GAM/log_"

from=1
to=4

for indN in 1
do
  for indB in 2
  do
    for indE in {1..3}
    do
      for indU in {1..5}
      do
        job_u=$(qsub -N design_gam_N${indN}_B${indB}_U${indU}_E${indE} -v indexN=$indN,indexB=$indB,indexU=$indU,indexE=$indE,path=$path,log_path=$log_path -J $from-$to sub_gam_designs.sub)
        echo $job_u
      done
    done
  done
done

from=1
to=5

for indN in 2
do
  for indB in 2
  do
    for indE in {1..3}
    do
      for indU in {1..5}
      do
        job_u=$(qsub -N design_gam_N${indN}_B${indB}_U${indU}_E${indE} -v indexN=$indN,indexB=$indB,indexU=$indU,indexE=$indE,path=$path,log_path=$log_path -J $from-$to sub_gam_designs.sub)
        echo $job_u
      done
    done
  done
done
