#!/bin/bash
path="codes/main_scripts/main_util_corollary.R"
log_path="log_files/util_coroll/log"

from=1
to=4

indN=1
for indB in 2
do
  for indE in {1..3}
  do
    for indU in {1..5}
    do
      job_u=$(qsub -N util_coll_N${indN}_B${indB}_U${indU}_E${indE} -v indexN=$indN,indexB=$indB,indexU=$indU,indexE=$indE,path=$path,log_path=$log_path -J $from-$to sub_gam_designs.sub)
      echo $job_u
    done
  done
done

from=1
to=5

indN=2
for indB in 2
do
  for indE in {1..3}
  do
    for indU in {1..5}
    do
      job_u=$(qsub -N util_coll_N${indN}_B${indB}_U${indU}_E${indE} -v indexN=$indN,indexB=$indB,indexU=$indU,indexE=$indE,path=$path,log_path=$log_path -J $from-$to sub_gam_designs.sub)
      echo $job_u
    done
  done
done