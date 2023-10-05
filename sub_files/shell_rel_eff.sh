#!/bin/bash
path="codes/main_scripts/main_relative_efficiency.R"
log="log_files/rel_eff/log"

for indN in {1..2}
do
  job1=$(qsub -N rel_eff_N${indN} -v indexN=$indN,path=$path,log=$log -J 1-3 sub_rel_eff.sub)
  echo $job1
done