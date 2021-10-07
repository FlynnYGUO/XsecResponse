#!/bin/bash
 
#for entry in $(ls ../m3_dune/inputs/mtuples/GenieCAF/HaddedFiles/*); 
#do
  #infile=${entry##*/}
  #outfile=${entry%.root}".txt"
  #com="../nusystematics/makeSelectedMtuples ${infile} > ${outfile}" 
  #eval $com
#done

#echo "Finished with splines"; 

#for entry in $(ls ../m3_dune/inputs/mtuples/GenieSampleCuts/*); 
#do
  #infile=${entry##*/}
  #outfile2=${entry%.root}".txt"
   #outfile3="./SplineOutputs/"${infile%.root}"_wgt.txt"
  #com2="../nusystematics/xsec_weight_gen --infile ${infile} --fhicl DUNETDRv1.ParamHeaders.noprolog.fcl &> ${outfile3} & "
  #eval $com2
  #sleep 5m; 
#done

declare -i i=0

for entry in $(ls ../m3_dune/inputs/mtuples/GenieSampleCuts/*); 
do
  declare -i j=0
  for entry2 in $(ls ../m3_dune/inputs/mtuples/WeightFiles/*); 
  do
    infile=${entry##*/}
    outfile2=${infile%.root}"_splines.root"
    outfile3="./SplineOutputs/"${infile%.root}"_splines.txt"
    com3="./bin/make_xsec_response_1d -w ${entry2} -m ${entry} -o ${outfile2} -selec numu &> ${outfile3} &"
    if [ $i == $j ]
    then
      eval ${com3}
      sleep 1m; 
    fi
    ((j++))
  done
  ((i++))
done
