// Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_bin/make_xsec_response_sk_2019.cc
//
// Usage:
//  ./bin/make_xsec_response_1d -w wtfile.root -m mtuple.root -o outfile.root -selec {nue, numu, cc1pi}
//
#include <stdlib.h>
#include <cstdlib>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <string.h>

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TAxis.h"

#include "XsecVary2D_NDGAr.h"

using std::cout;
using std::cerr;
using std::endl;

bool nue = false;
bool NHC = true;
bool byEvent_splines = false;
//std::string fTempFilenue =  "inputs/DUNE_numu_templates_v0.root";
//std::string fTempFilenumu =  "inputs/DUNE_numu_templates_v0.root";
//std::string fTempFilenue_erec =  "inputs/DUNE_numu_templates_v0_erec.root";
//std::string fTempFilenumu_erec =  "inputs/DUNE_numu_templates_v0_erec.root";
std::string fTempFile_yrec =  "inputs/DUNE_ND_templates_tdr_yrec.root";
std::string fTempFile_erec =  "inputs/DUNE_ND_templates_tdr_erec.root";
std::string fTempFile_etru =  "inputs/DUNE_ND_templates_tdr_erec.root";
//std::string fTempFile_etru =  "inputs/DUNE_numu_templates_v0_tdr_etru.root";

std::string fCovFile = "inputs/flux_covariance_13av1_prelim.root";
std::string wtfile = "";
std::string mtuple = "";
std::string outfile = "";
std::string thisSyst = "";
bool inputSyst = false;
std::string selection = "";
void Usage();
void ParseArgs(int argc, char **argv);
void ReadTemplates(TFile *ftemp, std::vector<TH1F*> &recoe_temps);


int main(int argc, char *argv[])
{
    //Get the arguments
    ParseArgs(argc, argv);

	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	std::cout << "Starting spline generation using " << wtfile << " and selection " << selection << std::endl;
    // Check name of file to see if it's numode or anti-numode, and whether nue or numu selection is required //!!PD NEEDS UPDATING FOR DUNE
    std::string mc_version(mtuple);

    //Change File paths
    //outfile = "/vols/dune/ljw20/DUNE_2023_ND_splines_split/" + outfile;
    std::cout << "THE OUTPUT FILE IS: " << outfile << std::endl; 

    //Open template files and covariance files
    TFile *ftemp_etru = new TFile(fTempFile_etru.c_str());
    TFile *ftemp_erec = new TFile(fTempFile_erec.c_str());
    TFile *ftemp_yrec = new TFile(fTempFile_yrec.c_str());

    //Load the templates
    std::vector<TH1F*> recoe_temps_etru;
    std::vector<TH1F*> recoe_temps_erec;
    std::vector<TH1F*> recoe_temps_yrec;
    ReadTemplates(ftemp_etru, recoe_temps_etru);
    ReadTemplates(ftemp_erec, recoe_temps_erec);
    ReadTemplates(ftemp_yrec, recoe_temps_yrec);
    std::cout << "Creating splines with ";
    std::cout << "Found " << recoe_temps_etru.size() << " etrue binning templates" << std::endl;
    std::cout << "Found " << recoe_temps_erec.size() << " erec binning templates" << std::endl;
    std::cout << "Found " << recoe_temps_yrec.size() << " yrec  binning templates" << std::endl;


    //Load the xsec variations, generate splines and graphs for xsec variations

    //For the nue sample ('0' is code for nue, '1' is code for numu). Actually this doesn't do anything in XsecVary at the moment, but could be used to put in the selections


    // ---------------------------------------------------------------------------------- //
    //                                                                                    //
    //               Start block for specifying systematics treatments.                   //
    //                                                                                    //
    // ---------------------------------------------------------------------------------- //

    // Vectors containing interaction mode codes that will eventually be associated with a
    // particular systematic. These are the modes this systematic affects.
    
    // CC+NC Mode lists
    std::vector<int> qemodes; qemodes.push_back(1); qemodes.push_back(15);
    std::vector<int> mecmodes; mecmodes.push_back(2); mecmodes.push_back(22);
    std::vector<int> resmodes; resmodes.push_back(4); resmodes.push_back(17);
    std::vector<int> dismodes; dismodes.push_back(3); dismodes.push_back(16);
    std::vector<int> disresmodes; disresmodes.push_back(3); disresmodes.push_back(4); disresmodes.push_back(16); disresmodes.push_back(17);
    std::vector<int> mostmodes; mostmodes.push_back(1); mostmodes.push_back(3); mostmodes.push_back(4);mostmodes.push_back(15); mostmodes.push_back(16); mostmodes.push_back(17);
    
    // CC modes only lists
    std::vector<int> ccqemodes; ccqemodes.push_back(1);
    std::vector<int> ccmecmodes; ccmecmodes.push_back(2);
    std::vector<int> ccresmodes; ccresmodes.push_back(4);
    std::vector<int> ccdismodes; ccdismodes.push_back(3);
    std::vector<int> ccdisresmodes; ccdisresmodes.push_back(3); ccdisresmodes.push_back(4);
    std::vector<int> ccmostmodes; ccmostmodes.push_back(1); ccmostmodes.push_back(3); ccmostmodes.push_back(4);
    
    //NC modes only list
    std::vector<int> ncqemodes; ncqemodes.push_back(15);
    std::vector<int> ncmecmodes; ncmecmodes.push_back(22);
    std::vector<int> ncresmodes; ncresmodes.push_back(17);
    std::vector<int> ncdismodes; ncdismodes.push_back(16);
    std::vector<int> ncdisresmodes; ncdisresmodes.push_back(16); disresmodes.push_back(17);
    std::vector<int> ncmostmodes; ncmostmodes.push_back(15); mostmodes.push_back(16); mostmodes.push_back(17);
    
    std::vector<int> ccallmodes; 
    std::vector<int> allmodes; 
    for(int i = 1; i <= 14; i++) 
     ccallmodes.push_back(i);

    for(int i = 1; i <= 26; i++) 
     allmodes.push_back(i);

    // Specify for which systematics we would like splines produced.
    // Constructor requires 4 inputs, with an optional 5th. In order:
    // (1) Name of systematic, formatted in correspondes with the input weight file
    //     tree name for that systematic, with "_tree" removed
    // (2) A short name for this systematic that will be featured in the name of
    //     output splines and graphs, along mode and bin information
    // (3) List of interaction modes affected by this parameter, built above
    // (4) The position of the nominal knot within the arrays contained in the
    //     supplied weight file
    // (5) [optional] Should extra knots be made to extend the range of systematic
    //      values in the splines, creating flat(ish) edges? (boolean)
   
   // DUNE SYSTEMATICS
 
//    SystematicProperties2D_NDGAr maqe_sysprop("MaCCQE","maqe",ccqemodes,3);  //1
//    SystematicProperties2D_NDGAr vecffqeshape_sysprop("VecFFCCQEshape","vecffqeshape", ccqemodes,0); //1
    //SystematicProperties2D_NDGAr mancel_sysprop("MaNCEL","mancel", allmodes,3); //? **
    //SystematicProperties2D_NDGAr etancel_sysprop("EtaNCEL","etancel", allmodes,3);//? **
	//ETA - this feels like a bug.... we're not giving any different modes to CC and NC events....
//    SystematicProperties2D_NDGAr mares_sysprop("MaCCRES","mares", resmodes,3);//4
//    SystematicProperties2D_NDGAr mvres_sysprop("MvCCRES","mvres", resmodes,3);//4
//    SystematicProperties2D_NDGAr mancres_sysprop("MaNCRES","mancres", resmodes,3);//4
//    SystematicProperties2D_NDGAr mvncres_sysprop("MvNCRES","mvncres", resmodes,3);//4
    //SystematicProperties2D_NDGAr rdecbr1gamma_sysprop("RDecBR1gamma", "rdecbr1gamma", resmodes,3);//4 **
    //SystematicProperties2D_NDGAr rdecbr1eta_sysprop("RDecBR1eta", "rdecbr1eta", resmodes,3);//4 **
    SystematicProperties2D_NDGAr thetadelta_sysprop("Theta_Delta2Npi", "thetadelta", resmodes,0); //4
    SystematicProperties2D_NDGAr ahtby_sysprop("AhtBY", "ahtby", dismodes,3); //3
    SystematicProperties2D_NDGAr bhtby_sysprop("BhtBY", "bhtby", dismodes,3); //3
    SystematicProperties2D_NDGAr cv1uby_sysprop("CV1uBY", "cv1uby", dismodes,3); //3
    SystematicProperties2D_NDGAr cv2uby_sysprop("CV2uBY", "cv2uby", dismodes,3); //3
    //SystematicProperties2D_NDGAr formzone_sysprop("FormZone", "formzone", dismodes,3); //3 **
    //SystematicProperties2D_NDGAr mfppi_sysprop("MFP_pi", "mfppi", dismodes,3);//3 **
	/*
    SystematicProperties2D_NDGAr frcexpi_sysprop("FrCEx_pi", "frcexpi", disresmodes,3);//3,4
    SystematicProperties2D_NDGAr frelaspi_sysprop("FrElas_pi", "frelaspi", disresmodes,3);//3,4
    SystematicProperties2D_NDGAr frinelpi_sysprop("FrInel_pi", "frinelpi", disresmodes,3);//3,4
    SystematicProperties2D_NDGAr frabspi_sysprop("FrAbs_pi", "frabspi", disresmodes,3);//3,4
    SystematicProperties2D_NDGAr frpiprodpi_sysprop("FrPiProd_pi", "frpiprodpi", disresmodes,3);//3,4
	*/

    SystematicProperties2D_NDGAr frcexpi_sysprop("FrCEx_pi", "frcexpi", allmodes,3);//3,4
//    SystematicProperties2D_NDGAr frelaspi_sysprop("FrElas_pi", "frelaspi", allmodes,3);//3,4
    SystematicProperties2D_NDGAr frinelpi_sysprop("FrInel_pi", "frinelpi", allmodes,3);//3,4
    SystematicProperties2D_NDGAr frabspi_sysprop("FrAbs_pi", "frabspi", allmodes,3);//3,4
    SystematicProperties2D_NDGAr frpiprodpi_sysprop("FrPiProd_pi", "frpiprodpi", allmodes,3);//3,4

    //SystematicProperties2D_NDGAr mfpn_sysprop("MFP_N", "mfpn", mostmodes,3);//1,3,4,10 **
    SystematicProperties2D_NDGAr frcexn_sysprop("FrCEx_N", "frcexn", mostmodes,3);//1,3,4,10
//    SystematicProperties2D_NDGAr frelasn_sysprop("FrElas_N", "frelasn", mostmodes,3);//1,3,4,10
    SystematicProperties2D_NDGAr frineln_sysprop("FrInel_N", "frineln", mostmodes,3);//1,3,4,10
    SystematicProperties2D_NDGAr frabsn_sysprop("FrAbs_N", "frabsn", mostmodes,3);//1,3,4,10
    SystematicProperties2D_NDGAr frpiprodn_sysprop("FrPiProd_N", "frpiprodn", mostmodes,3);//1,3,4,10
//    SystematicProperties2D_NDGAr pauli_sysprop("CCQEPauliSupViaKF", "pauli", ccqemodes,0);//1
    //SystematicProperties2D_NDGAr gauss_sysprop("Mnv2p2hGaussEnhancement", "gauss", allmodes,3);//1,3,4,5,10 *
    //SystematicProperties2D_NDGAr mkspp_sysprop("MKSPP_ReWeight", "mkspp", allmodes,3);//1,3,4,5,10**
    SystematicProperties2D_NDGAr e2anu_sysprop("E2p2h_A_nu", "e2anu", mecmodes,2);//1,3,4,5,7,10 
    SystematicProperties2D_NDGAr e2bnu_sysprop("E2p2h_B_nu", "e2bnu", mecmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr e2anubar_sysprop("E2p2h_A_nubar", "e2anubar", mecmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr e2bnubar_sysprop("E2p2h_B_nubar", "e2bnubar", mecmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nuncc2_sysprop("NR_nu_n_CC_2Pi", "nuncc2", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nuncc3_sysprop("NR_nu_n_CC_3Pi", "nuncc3", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nupcc2_sysprop("NR_nu_p_CC_2Pi", "nupcc2", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nupcc3_sysprop("NR_nu_p_CC_3Pi", "nupcc3", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nunpcc1_sysprop("NR_nu_np_CC_1Pi", "nunpcc1", allmodes,3);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nunnc1_sysprop("NR_nu_n_NC_1Pi", "nunnc1", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nunnc2_sysprop("NR_nu_n_NC_2Pi", "nunnc2", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nunnc3_sysprop("NR_nu_n_NC_3Pi", "nunnc3", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nupnc1_sysprop("NR_nu_p_NC_1Pi", "nupnc1", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nupnc2_sysprop("NR_nu_p_NC_2Pi", "nupnc2", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nupnc3_sysprop("NR_nu_p_NC_3Pi", "nupnc3", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarncc1_sysprop("NR_nubar_n_CC_1Pi", "nubarncc1", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarncc2_sysprop("NR_nubar_n_CC_2Pi", "nubarncc2", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarncc3_sysprop("NR_nubar_n_CC_3Pi", "nubarncc3", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarpcc1_sysprop("NR_nubar_p_CC_1Pi", "nubarpcc1", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarpcc2_sysprop("NR_nubar_p_CC_2Pi", "nubarpcc2", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarpcc3_sysprop("NR_nubar_p_CC_3Pi", "nubarpcc3", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarnnc1_sysprop("NR_nubar_n_NC_1Pi", "nubarnnc1", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarnnc2_sysprop("NR_nubar_n_NC_2Pi", "nubarnnc2", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarnnc3_sysprop("NR_nubar_n_NC_3Pi", "nubarnnc3", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarpnc1_sysprop("NR_nubar_p_NC_1Pi", "nubarpnc1", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarpnc2_sysprop("NR_nubar_p_NC_2Pi", "nubarpnc2", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nubarpnc3_sysprop("NR_nubar_p_NC_3Pi", "nubarpnc3", allmodes,2);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr berpaa_sysprop("BeRPA_A", "berpaa", ccqemodes,3);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr berpab_sysprop("BeRPA_B", "berpab", ccqemodes,3);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr berpad_sysprop("BeRPA_D", "berpad", ccqemodes,3);//1,3,4,5,7,10
    //SystematicProperties2D_NDGAr berpae_sysprop("BeRPA_E", "berpae", allmodes,2);//1,3,4,5,7,10 **
    //SystematicProperties2D_NDGAr lepmom_sysprop("EbFSLepMomShift", "lepmom", allmodes,3);//1,3,4,5,7,10 **
    SystematicProperties2D_NDGAr c12nu_sysprop("C12ToAr40_2p2hScaling_nu", "c12nu", mecmodes,0);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr c12nubar_sysprop("C12ToAr40_2p2hScaling_nubar", "c12nubar", mecmodes,0);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nuexsec_sysprop("nuenuebar_xsec_ratio", "nuexsec", allmodes,1);//1,3,4,5,7,10
    SystematicProperties2D_NDGAr nuemuxsec_sysprop("nuenumu_xsec_ratio", "nuemuxsec", allmodes,1);//1,3,4,5,7,10
    //SystematicProperties2D_NDGAr q2sup_sysprop("SPPLowQ2Suppression", "q2sup", allmodes,3);//1,3,4,5,10 **
    //SystematicProperties2D_NDGAr fsismear_sysprop("FSILikeEAvailSmearing", "fsismear", allmodes,3);//1,3,4,5,7,10 **


    // Put these systematics in a list
    std::vector<SystematicProperties2D_NDGAr> systProps;   	
//    systProps.push_back(maqe_sysprop);
//    systProps.push_back(vecffqeshape_sysprop);
    //systProps.push_back(mancel_sysprop);
    //systProps.push_back(etancel_sysprop);
//    systProps.push_back(mares_sysprop);
//    systProps.push_back(mvres_sysprop);
//    systProps.push_back(mancres_sysprop);
//    systProps.push_back(mvncres_sysprop);
    //systProps.push_back(rdecbr1gamma_sysprop);
    //systProps.push_back(rdecbr1eta_sysprop);
    systProps.push_back(thetadelta_sysprop);
    systProps.push_back(ahtby_sysprop);
    systProps.push_back(bhtby_sysprop);
    systProps.push_back(cv1uby_sysprop);
    systProps.push_back(cv2uby_sysprop);
    //systProps.push_back(formzone_sysprop);
    //systProps.push_back(mfppi_sysprop);
    systProps.push_back(frcexpi_sysprop);
//    systProps.push_back(frelaspi_sysprop);
    systProps.push_back(frinelpi_sysprop);
    systProps.push_back(frabspi_sysprop);
    systProps.push_back(frpiprodpi_sysprop);
    //systProps.push_back(mfpn_sysprop);
    systProps.push_back(frcexn_sysprop);
//    systProps.push_back(frelasn_sysprop);
    systProps.push_back(frineln_sysprop);
    systProps.push_back(frabsn_sysprop);
    systProps.push_back(frpiprodn_sysprop);
//    systProps.push_back(pauli_sysprop);
    //systProps.push_back(gauss_sysprop);
    //systProps.push_back(mkspp_sysprop);
    systProps.push_back(e2anu_sysprop);
    systProps.push_back(e2bnu_sysprop);
    systProps.push_back(e2anubar_sysprop);
    systProps.push_back(e2bnubar_sysprop);
    systProps.push_back(nuncc2_sysprop);
    systProps.push_back(nuncc3_sysprop);
    systProps.push_back(nupcc2_sysprop);
    systProps.push_back(nupcc3_sysprop);
    systProps.push_back(nunpcc1_sysprop);
    systProps.push_back(nunnc1_sysprop);
    systProps.push_back(nunnc2_sysprop);
    systProps.push_back(nunnc3_sysprop);
    systProps.push_back(nupnc1_sysprop);
    systProps.push_back(nupnc2_sysprop);
    systProps.push_back(nupnc3_sysprop);
    systProps.push_back(nubarncc1_sysprop);
    systProps.push_back(nubarncc2_sysprop);
    systProps.push_back(nubarncc3_sysprop);
    systProps.push_back(nubarpcc1_sysprop);
    systProps.push_back(nubarpcc2_sysprop);
    systProps.push_back(nubarpcc3_sysprop);
    systProps.push_back(nubarnnc1_sysprop);
    systProps.push_back(nubarnnc2_sysprop);
    systProps.push_back(nubarnnc3_sysprop);
    systProps.push_back(nubarpnc1_sysprop);
    systProps.push_back(nubarpnc2_sysprop);
    systProps.push_back(nubarpnc3_sysprop);
    systProps.push_back(berpaa_sysprop);
    systProps.push_back(berpab_sysprop);
    systProps.push_back(berpad_sysprop);
    //systProps.push_back(berpae_sysprop);
    //systProps.push_back(lepmom_sysprop);
    systProps.push_back(c12nu_sysprop);
    systProps.push_back(c12nubar_sysprop);
    systProps.push_back(nuexsec_sysprop);
    systProps.push_back(nuemuxsec_sysprop);
    //systProps.push_back(q2sup_sysprop);
    //systProps.push_back(fsismear_sysprop);	

    if (inputSyst) {
      std::vector<SystematicProperties2D_NDGAr> systPropsInput;
	  int numSysts = systProps.size();
      for (int i = 0; i < numSysts; i++) {
        if (systProps[i].shortName == thisSyst) {systPropsInput.push_back(systProps[i]);}
      }
    // Call the XSecVary class constructor
	systProps = systPropsInput;
	}

    XsecVary2D_NDGAr xs(wtfile, mtuple, systProps);

    // ---------------------------------------------------------------------------------- //
    //                                                                                    //
    //                End block for specifying systematics treatments.                    //
    //                                                                                    //
    // ---------------------------------------------------------------------------------- //

 
    // For binned splines: set binning then make variations.
    if (byEvent_splines == false)
      {
	// Load the flux covariance matrix and extract energy binning
	// Determine which true energy binning based on the first event
	// in the file
	xs.GetEntry(0);

	const int ebins = recoe_temps_etru[0]->GetXaxis()->GetNbins();
    const Double_t *ebinrange = recoe_temps_etru[0]->GetXaxis()->GetXbins()->GetArray();
	std::cout << "Number of true energy bins: " << ebins << "." << std::endl;
	std::cout << "True energy bins: {";
	for (int ibin=0; ibin<ebins; ibin++)
	  std::cout << ebinrange[ibin] << ", ";
	std::cout << ebinrange[ebins] << "}" << std::endl;
	
	const int rebins = recoe_temps_erec[0]->GetXaxis()->GetNbins();
	const Double_t *rebinrange = recoe_temps_erec[0]->GetXaxis()->GetXbins()->GetArray();
    std::cout << "Number of reco energy bins: " << rebins << "." << std::endl;
	std::cout << "Reco energy bins: {";
	for (int ibin=0; ibin<rebins; ibin++)
	  std::cout << rebinrange[ibin] << ", ";
	std::cout << rebinrange[rebins] << "}" << std::endl;
 
	const int rybins = recoe_temps_yrec[0]->GetXaxis()->GetNbins();
    const Double_t *rybinrange = recoe_temps_yrec[0]->GetXaxis()->GetXbins()->GetArray();
    std::cout << "Number of reco y bins: " << rybins << "." << std::endl;
	std::cout << "Reco y bins: {";
	for (int ibin=0; ibin<rybins; ibin++)
	  std::cout << rybinrange[ibin] << ", ";
	std::cout << rybinrange[rybins] << "}" << std::endl;
	
	xs.SetBinning(ebinrange,ebins,rebinrange,rebins,rybinrange,rybins);
	
	std::cout <<  "Binning implemented successfully" << std::endl; 

	xs.MakeVariations();

	xs.WriteGraphs(outfile);
        
      }
/*
    else // Event-by-event splines 
      {
	xs.InitByEventTree(13, outfile);
	xs.MakeVariations_byEvent();

	xs.WriteByEventTree();
      }
*/

    return 0;
}

// Print the cmd line syntax
void Usage(){
    cout << "Cmd line syntax should be:" << endl;
    cout << "./bin/make_xsec_response.exe [-t 'byEv' (for event by event splines)] -w wtfile -m mtuple -o outfile [-nue (for nue selection, as opposed to numu)]" << endl;
}

// Messy way to process cmd line arguments.
void ParseArgs(int argc, char **argv){
    //  std::cout << "parse args" << std::endl;
    int nargs = 1; 
    if(argc<(nargs*2+1)){ Usage(); exit(1); }
    for(int i = 1; i < argc; i+=2){
        if(std::string(argv[i]) == "-w") wtfile = argv[i+1];
        else if(std::string(argv[i]) == "-m") mtuple = argv[i+1];
        else if(std::string(argv[i]) == "-o") outfile = argv[i+1];
        else if(std::string(argv[i]) == "-s") {thisSyst = argv[i+1]; inputSyst = true;}
	else if(std::string(argv[i]) == "-t") 
	  {
	    if (std::string(argv[i+1])=="byEv") byEvent_splines=true;
	    else byEvent_splines=false;
	  }
	else if(std::string(argv[i]) == "-f")
	  {
	    if(std::string(argv[i+1])=="fine")
	      {
	      }
	  }
        else 
	  {  
            cout << "Invalid argument:" << argv[i] << " "<< argv[i+1] << endl;
            Usage();
            exit(1);
	  }
    } 
}


void ReadTemplates(TFile *ftemp, std::vector<TH1F*> &recoe_temps){
    bool found_temp = true;
    int temp_iter=0;
    while(found_temp){
        char temp_name[50];
        sprintf(temp_name,"hrecoe_%d",temp_iter); 
        TH1F *blah = NULL;
        blah = (TH1F*)(ftemp->Get(temp_name));
        if(blah!=NULL){
            recoe_temps.push_back(blah);
            temp_iter++;
	    std::cout << temp_name << std::endl
		      << blah << std::endl
		      << blah->GetName() << std::endl
		      << blah->GetXaxis()->GetNbins()
		      << std::endl;
        } else found_temp = false;
    }
}
