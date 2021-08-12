#!/bin/bash

#weightfolder=/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/genWeights_mtuples_19b_procv1_v4
#/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/genweight_neut_5.4.0_v6_cc1pi_berpaweights_flux13av2_Q2_params
#/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/genWeights_mtuples_19b_procv1_v2
#weightfolder=/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/asg2019oa/SKmtuples/SKweights
#weightfolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/weights_for_splines_090820
#weightfolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/weights_19b_15122020
weightfolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/weights_19b_15122020_no_mom_cuts
#
#
#/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/genweight_neut_5.4.0_first_check_28122019
weightprefix=""
#outfolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/splines_first_try_090820
#outfolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/splines_19b_production_24112020
#outfolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/splines_19b_production_noCoulombCorrections_24112020
#outfolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/splines_19b_15122020
outfolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/splines_19b_15122020_no_mom_cuts
#/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/splines_fitqun_cc1pi_berpaweighs_flux13av2_v3
#/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/splines_mtuples_19b_procv1_v3
##
#/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/splines2019_v1 #!!CHANGE THESE TWO EACH TIME
splineprefix=2019splinesfitqunPre #!!CHANGE THESE TWO EACH TIME
#mtuplefolder=/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/mtuples_19b_procv4
#mtuplefolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/old_production_yoshidasan_selection
#mtuplefolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/19b_production_24112020
#mtuplefolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/19b_production_15122020
mtuplefolder=/vols/t2k/users/ea2817/MaCh3_OA_numuCC1pi/19b_production_15122020_no_mom_cuts
#/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/mtuples_fitqun_cc1pi_berpaweights_flux13av2
#/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/mtuples_19b_procv1
#
#/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/mtuples_19b_procv1
#/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/sk_mtuple_run1_9dflux_incnormhists
# For older splines production
#/QMULZone1/home/asg/asg2017oa/SK/mtuples_fitqun_cc1pi_berpaweights_flux13av2

#weightfolder=/afs/cern.ch/work/b/bradics/T2K/NEUT/T2KReWeight-build/mtuples/Kevin/weights-final
#weightprefix=2018weightsfitqunv2
#outfolder=/afs/cern.ch/work/b/bradics/T2K/NEUT/T2KReWeight-build/mtuples/output/v2 #!!CHANGE THESE TWO EACH TIME
#splineprefix=2018splinesfitqunv2 #!!CHANGE THESE TWO EACH TIME
#mtuplefolder=/afs/cern.ch/work/b/bradics/T2K/NEUT/T2KReWeight-build/mtuples/

mtupleprefix=t2ksk19b.fqv4r0
#t2ksk.14a.neut5.3.2.13a_tuned_v1r1
#
#
#t2ksk.14a.neut5.3.2.13a_tuned_v1r1
isfitqun="true"
XSEC_RESPONSE="/vols/t2k/users/ea2817/build/XsecResponse/XsecResponse"
qsub_opt="-q hep.q -l h_rt=2:00:0 -l h_vmem=4G"
STD_ERR="${outfolder}/err"
STD_OUT="${outfolder}/out"

mkdir $outfolder
mkdir $STD_ERR
mkdir $STD_OUT
for horn in fhc rhc
#for horn in 250ka.fine m250ka.fine
do
  for flavour1 in numu nue
  do
	for flavour2 in nue numu
	do
	  if [[ "$flavour1" == 'nue' ]]
	  then
		if [[ "$flavour2" == 'numu' || "$flavour2" == 'numubar' || "$flavour2" == 'nuebar' ]]
		then
		  continue
		fi
	  elif [[ "$flavour1" == 'nuebar' ]]
	  then
		if [[ "$flavour2" == 'numubar' || "$flavour2" == 'numu' || "$flavour2" == 'nue' ]]
		then
		  continue
		fi
	  elif [[ "$flavour1" == 'numu' ]]
	  then
		if [[ "$flavour2" == 'numubar' || "$flavour2" == 'nuebar' ]]
		then
		  continue
		fi
	  elif [[ "$flavour1" == 'numubar' ]]
	  then
		if [[ "$flavour2" == 'numu' || "$flavour2" == 'nue' ]]
		then
		  continue
		fi
	  fi
	  for sel in nue numu cc1pi numucc1pi1de numucc1pi2de
	  #for sel in nue numu cc1pi
	  do
		if [[ -f ${weightfolder}/${weightprefix}${mtupleprefix}.${horn}.${flavour1}_x_${flavour2}_${sel}selec_spline_weights.root ]]
		then
		  script_base_nu="${outfolder}/spline_gen_${horn}_${flavour1}_x_${flavour2}_${sel}sel"
		  script_base_nubar="${outfolder}/spline_gen_${horn}_${flavour1}bar_x_${flavour2}bar_${sel}sel"

		  script_name_nu_2D="${script_base_nu}_erectheta.sh"
		  script_name_nubar_2D="${script_base_nubar}_erectheta.sh"
		  script_name_nu_1D="${script_base_nu}.sh"
		  script_name_nubar_1D="${script_base_nubar}.sh"

		  echo "cd $XSEC_RESPONSE" > ${script_name_nu_2D}
		  #echo "source setup_ed.sh" >> ${script_name_nu_2D}

		  echo "cd $XSEC_RESPONSE" > ${script_name_nubar_2D}
		  #echo "source setup_ed.sh" >> ${script_name_nubar_2D}

		  echo "cd $XSEC_RESPONSE" > ${script_name_nu_1D}
		  #echo "source setup_ed.sh" >> ${script_name_nu_1D}

		  echo "cd $XSEC_RESPONSE" > ${script_name_nubar_1D}
		  #echo "source setup_ed.sh" >> ${script_name_nubar_1D}


		  #2D
		  #nu
		  echo "./bin/make_xsec_response_sk_2019_2d -w ${weightfolder}/${weightprefix}${mtupleprefix}.${horn}.${flavour1}_x_${flavour2}_${sel}selec_spline_weights.root -m ${mtuplefolder}/${mtupleprefix}.${horn}.${flavour1}_x_${flavour2}_${sel}selec.root -o $outfolder/spline${splineprefix}erectheta_${mtupleprefix}.${horn}.${flavour1}_x_${flavour2}_${sel}selec.root -selec $sel -f $isfitqun" >> ${script_name_nu_2D}
		  #nubar
		  echo "./bin/make_xsec_response_sk_2019_2d -w ${weightfolder}/${weightprefix}${mtupleprefix}.${horn}.${flavour1}bar_x_${flavour2}bar_${sel}selec_spline_weights.root -m ${mtuplefolder}/${mtupleprefix}.${horn}.${flavour1}bar_x_${flavour2}bar_${sel}selec.root -o $outfolder/spline${splineprefix}erectheta_${mtupleprefix}.${horn}.${flavour1}bar_x_${flavour2}bar_${sel}selec.root -selec $sel -f $isfitqun" >> ${script_name_nubar_2D}

		  #1D
		  #nu
		  echo "./bin/make_xsec_response_sk_2019 -w ${weightfolder}/${weightprefix}${mtupleprefix}.${horn}.${flavour1}_x_${flavour2}_${sel}selec_spline_weights.root -m ${mtuplefolder}/${mtupleprefix}.${horn}.${flavour1}_x_${flavour2}_${sel}selec.root -o $outfolder/spline${splineprefix}_${mtupleprefix}.${horn}.${flavour1}_x_${flavour2}_${sel}selec.root -selec $sel" -f $isfitqun >> ${script_name_nu_1D}
		  #nubar
		  echo "./bin/make_xsec_response_sk_2019 -w ${weightfolder}/${weightprefix}${mtupleprefix}.${horn}.${flavour1}bar_x_${flavour2}bar_${sel}selec_spline_weights.root -m ${mtuplefolder}/${mtupleprefix}.${horn}.${flavour1}bar_x_${flavour2}bar_${sel}selec.root -o $outfolder/spline${splineprefix}_${mtupleprefix}.${horn}.${flavour1}bar_x_${flavour2}bar_${sel}selec.root -selec $sel" -f $isfitqun >> ${script_name_nubar_1D}


		  com_nu_2D="qsub -e ${outfolder}/err/${script_name_nu_2D}.err -o ${outfolder}/out/${script_name_nu_2D}.out ${qsub_opt} ${script_name_nu_2D}"
		  com_nubar_2D="qsub -e ${outfolder}/err/${script_name_nubar_2D}.err -o ${outfolder}/out/${script_name_nubar_2D}.out ${qsub_opt} ${script_name_nubar_2D}"
		  com_nu_1D="qsub -e ${outfolder}/err/${script_name_nu_1D}.err -o ${outfolder}/out/${script_name_nu_1D}.out ${qsub_opt} ${script_name_nu_1D}"
		  com_nubar_1D="qsub -e ${outfolder}/err/${script_name_nubar_1D}.err -o ${outfolder}/out/${script_name_nubar_1D}.out ${qsub_opt} ${script_name_nubar_1D}"

		  #eval $com_nu_2D
		  echo "submitted $flavour1 x $flavour2 2D"
		  sleep 1
		  #eval $com_nubar_2D
		  echo "submitted ${flavour1}bar x ${flavour2}bar 2D"
		  sleep 1
		  #eval $com_nu_1D
		  echo "submitted ${flavour1} x ${flavour2} 1D"
		  sleep 1
		  #eval $com_nubar_1D
		  echo "submitted ${flavour1}bar x ${flavour2}bar 1D"
		  sleep 1
		  
		else
		  echo "Couldn't find file ${weightfolder}/${weightprefix}${mtupleprefix}.${horn}.${flavour1}_x_${flavour2}_${sel}selec.root"

		fi
	  done
	done
  done
done

