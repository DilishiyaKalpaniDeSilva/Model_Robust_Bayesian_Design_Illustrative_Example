#!/bin/bash
path="codes/main_scripts/main_polynomial_designs.R"
logp="log_files/designs/polynomial/log"

log1="fixed"
log2="all"

for indN in {1..2}
do
  job1=$(qsub -N design_polynomialF_N${indN} -v indexN=$indN,path=$path,logp=$logp,log=$log1 -J 1-3 sub_polynomial.sub)
  echo $job1
  job2=$(qsub -N design_polynomial_N${indN} -v indexN=$indN,path=$path,logp=$logp,log=$log2 -J 1-3 sub_polynomial.sub)
  echo $job2
done
