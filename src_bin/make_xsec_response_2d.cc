// Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_bin/make_xsec_response_sk_2019_2d.cc
//
// Usage:
//  ./bin/make_xsec_response_2d -w wtfile.root -m mtuple.root -o outfile.root -selec {nue, numu, cc1pi} -f {true, false}
// where -f {true, false} specifies whether or not fitqun is used to produce the mtuple
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
#include "TH3F.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TAxis.h"

#include "XsecVary2D.h"

using std::cout;
using std::cerr;
using std::endl;

bool nue = false;
bool NHC = true;
bool byEvent_splines = false;
bool isfitqun = false;
std::string fCovFile = "inputs/flux_covariance_13av1_prelim.root";
char * wtfile = NULL;
char * mtuple = NULL;
char * outfile = NULL;
char * selection = NULL;

char * treename = NULL;
char * horn = NULL;
char * var1name = NULL;
char * var2name = NULL;


void Usage();
void ParseArgs(int argc, char **argv);
void ReadTemplates(TFile *ftemp, std::vector<TH1F*> &recoe_temps);


int main(int argc, char *argv[])
{
    //Get the arguments
    ParseArgs(argc, argv);

	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	std::cout << "Starting spline generation using " << wtfile << " and selection " << selection << std::endl;
    // Check name of file to see if it's numode or anti-numode, and whether nue or numu selection is required //!!PD NEEDS ADAPTING TO DUNE
    std::string RHC_mode("m250");
    std::string alt_RHC_mode("rhc");
    std::string mc_version(mtuple);
    std::string nue_selec("nue");
    std::string numu_selec("numu");
    std::string cc1pi_selec("cc1pi");
    std::string str_selec(selection);
    std::size_t mc_found_RHC = mc_version.find(RHC_mode);
    std::size_t alt_mc_found_RHC = mc_version.find(alt_RHC_mode);
    if (mc_found_RHC != std::string::npos || alt_mc_found_RHC != std::string::npos)
      NHC = false;
    std::size_t selec_found_nue = str_selec.find(nue_selec);
    std::size_t selec_found_numu = str_selec.find(numu_selec);
    std::size_t selec_found_cc1pi = str_selec.find(cc1pi_selec);
    if (selec_found_nue != std::string::npos)
      nue = true;
    else if (selec_found_numu != std::string::npos)
      nue = false;
	else if (selec_found_cc1pi != std::string::npos)
      nue = true;
    else{
      std::cerr << "ERROR 2D: Incompatible selection given: " << selection << ". Something has gone wrong..." << std::endl;
      std::cerr << "ERROR 2D: You gave me a weight file " << wtfile << std::endl;
      std::cerr << "ERROR 2D: You gave me a mtuple file " << mtuple << std::endl;
      std::cerr << "ERROR 2D: You gave me a selection of " << selection << std::endl;
	  throw;
	}

    std::cout << "Detected from filename: Numode = " << NHC << std::endl;
    if (nue) std::cout << "Using nue selection" << std::endl;
    else std::cout << "Using numu selection" << std::endl;



    //set parameters and binning for numu and nue
    std::string  var1numu="erec";
    const int var1numubins=58;
    double var1numubinrange[var1numubins+1];
    for(int i=0;i<31; i++){
      var1numubinrange[i] = 1.5/double(30)*i;
    }
    for(int i=31;i<46; i++){
      var1numubinrange[i] = 1.5+(0.1*double(i-30));
    }
    var1numubinrange[46] = 3.25;
    var1numubinrange[47] = 3.5;
    var1numubinrange[48] = 3.75;
    var1numubinrange[49] = 4.0;
    var1numubinrange[50] = 4.5;
    var1numubinrange[51] = 5.0;
    var1numubinrange[52] = 5.5;
    var1numubinrange[53] = 6.0;
    var1numubinrange[54] = 7.0;
    var1numubinrange[55] = 8.0;
    var1numubinrange[56] = 9.0;
    var1numubinrange[57] = 10.0;
    var1numubinrange[58] = 30.0;
    std::string  var2numu="theta";
    //    const int var2numubins=1;
    //double var2numubinrange[var1numubins+1]={0.0, 180.0};
    const int var2numubins=15;
    double var2numubinrange[var1numubins+1]={0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 180.0};

    std::string var1nue="erec";
    const int var1nuebins=25;
    double var1nuebinrange[var1nuebins+1];
    for(int iBin=0;iBin<var1nuebins+1;iBin++){
      var1nuebinrange[iBin]=1.25/double(var1nuebins)*iBin;
    }
    std::string var2nue="theta";
    const int var2nuebins=15;
    double var2nuebinrange[var1nuebins+1]={0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 180.0};


    std::cout << "Creating splines with ";
    if (nue) std::cout << "nue binning" << std::endl;
    else std::cout << "numu binning" << std::endl;


    //Load the xsec variations, generate splines and graphs for xsec variations

    //For the nue sample ('0' is code for nue, '1' is code for numu). Actually this doesn't do anything in XsecVary at the moment, but could be used to put in the selections
    int nutype;
    if (nue) nutype = 0;
    else nutype = 1;


    // ---------------------------------------------------------------------------------- //
    //                                                                                    //
    //               Start block for specifying systematics treatments.                   //
    //                                                                                    //
    // ---------------------------------------------------------------------------------- //

    // Vectors containing interaction mode codes that will eventually be associated with a
    // particular systematic. These are the modes this systematic affects.  
    std::vector<int> ccqemode; ccqemode.push_back(0);
    std::vector<int> mecmode; mecmode.push_back(8);
    std::vector<int> pimodes; pimodes.push_back(1); pimodes.push_back(4); pimodes.push_back(5);
    std::vector<int> ccothmode; ccothmode.push_back(3);
	std::vector<int> ccmpimode; ccmpimode.push_back(10);
	std::vector<int> ccdismode; ccdismode.push_back(11);
    std::vector<int> allmodes; 
    for(int i = 0; i <= 11; i++) 
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
	//
    //    SystematicProperties maqe_sf_sysprop("NXSec_MaCCQE_SF","maqesf",ccqemode,14);
	//    SystematicProperties2D pfo_sysprop("NIWG2014a_pF_O16","pfo",ccqemode,6,true);
	//    SystematicProperties2D dismpishp_sysprop("NIWG2012a_dismpishp","dismpishp",ccothmode,6);
	//    SystematicProperties2D berpaA_sysprop("NIWG_Effective_rpaCCQE_A","erpa_A",ccqemode,6);
	//    SystematicProperties2D berpaB_sysprop("NIWG_Effective_rpaCCQE_B","erpa_B",ccqemode,6);
	//    SystematicProperties2D berpaC_sysprop("NIWG_Effective_rpaCCQE_C","erpa_D",ccqemode,6);
	//    SystematicProperties2D berpaD_sysprop("NIWG_Effective_rpaCCQE_D","erpa_E",ccqemode,6);

    SystematicProperties2D maqe_sysprop("NXSec_MaCCQE","maqe",ccqemode,3);
    SystematicProperties2D mec_sysprop("NIWGMEC_PDDWeight_O16","mecpdd",mecmode,3);
    SystematicProperties2D ca5_sysprop("NXSec_CA5RES","ca5",pimodes,3);
    SystematicProperties2D manff_sysprop("NXSec_MaRES","manff",pimodes,3);
	SystematicProperties2D bgscl_sysprop("NXSec_BgSclRES","bgscl",pimodes,3);

	// only 7 knots!
	SystematicProperties2D FSI_NCasc_FrAbs_pi_sysprop("NCasc_FrAbs_pi","FSIPIABS",allmodes,3);
	SystematicProperties2D FSI_NCasc_FrInelLow_pi_sysprop("NCasc_FrInelLow_pi","FSIINELLO",allmodes,3);
	SystematicProperties2D FSI_NCasc_FrCExLow_pi_sysprop("NCasc_FrCExLow_pi","FSICEXLO",allmodes,3);
	//SystematicProperties2D FSI_NCasc_FrCExHigh_pi_sysprop("NCasc_FrCExHigh_pi","FSICEXHI",allmodes,3);//ETA- NOT USUALLY IN THE OA
	SystematicProperties2D FSI_NCasc_FrInelHigh_pi_sysprop("NCasc_FrInelHigh_pi","FSIINELHI",allmodes,3);
	SystematicProperties2D FSI_NCasc_FrPiProd_pi_sysprop("NCasc_FrPiProd_pi","FSIPIPROD",allmodes,1);

	//ETA CCDIS and CCMPi parameters for 2020OA
	SystematicProperties2D BY_DIS_sysprop("NIWG_DIS_BY", "DISBY", ccdismode,3);
	SystematicProperties2D BY_MPi_sysprop("NIWG_MultiPi_BY", "MPiBY", ccmpimode,3);
	SystematicProperties2D Xsec_AGKY_sysprop("NIWG_MultiPi_Xsec_AGKY", "MPiAGKYXsec", ccmpimode,3);
	//ETA adding I-1/2 anti-nu parameter
	SystematicProperties2D bgsclbar_sysprop("NXSec_BgSclLMCPiBarRES", "bgsclbar", pimodes, 3);
	//ETA adding 2p2h which wasn't in this branch...
	SystematicProperties2D Edep2p2h_lowEnu_sysprop("NIWG_2p2hEdep_lowEnu","2p2hedeplowenu", mecmode, 1);
	SystematicProperties2D Edep2p2h_highEnu_sysprop("NIWG_2p2hEdep_highEnu", "2p2hedephienu", mecmode, 1);
	SystematicProperties2D Edep2p2h_lowEnubar_sysprop("NIWG_2p2hEdep_lowEnubar", "2p2hedeplowenubar", mecmode, 1);
	SystematicProperties2D Edep2p2h_highEnubar_sysprop("NIWG_2p2hEdep_highEnubar", "2p2hedephienubar", mecmode, 1);
	//ETA adding low Q2 params for VALOR                                                                                                                                                                      
	//Changing the naming to be conistent with ND
	//i.e. C++ indexing on the names!
	SystematicProperties2D LowQ2Suppression_1_sysprop("NIWGQETwk_LowQ2Suppression1", "LowQ2Suppression0", ccqemode, 3);
	SystematicProperties2D LowQ2Suppression_2_sysprop("NIWGQETwk_LowQ2Suppression2", "LowQ2Suppression1", ccqemode, 3);
	SystematicProperties2D LowQ2Suppression_3_sysprop("NIWGQETwk_LowQ2Suppression3", "LowQ2Suppression2", ccqemode, 3);
	SystematicProperties2D LowQ2Suppression_4_sysprop("NIWGQETwk_LowQ2Suppression4", "LowQ2Suppression3", ccqemode, 3);
	SystematicProperties2D LowQ2Suppression_5_sysprop("NIWGQETwk_LowQ2Suppression5", "LowQ2Suppression4", ccqemode, 3);
	SystematicProperties2D LowQ2Suppression_6_sysprop("NIWGQETwk_LowQ2Suppression6", "LowQ2Suppression5", ccqemode, 3);
	SystematicProperties2D LowQ2Suppression_7_sysprop("NIWGQETwk_LowQ2Suppression7", "LowQ2Suppression6", ccqemode, 3);
	SystematicProperties2D LowQ2Suppression_8_sysprop("NIWGQETwk_LowQ2Suppression8", "LowQ2Suppression7", ccqemode, 3);




    // Put these systematics in a list
    std::vector<SystematicProperties2D> systProps;
    systProps.push_back(maqe_sysprop);
    //systProps.push_back(maqe_sf_sysprop);
    systProps.push_back(mec_sysprop);
	//    systProps.push_back(pfo_sysprop);
    systProps.push_back(ca5_sysprop);
    systProps.push_back(manff_sysprop);
	systProps.push_back(bgscl_sysprop);
	//    systProps.push_back(dismpishp_sysprop);
	//    systProps.push_back(berpaA_sysprop);
	//    systProps.push_back(berpaB_sysprop);
	//    systProps.push_back(berpaC_sysprop);
	//    systProps.push_back(berpaD_sysprop);
    systProps.push_back(FSI_NCasc_FrAbs_pi_sysprop);
    systProps.push_back(FSI_NCasc_FrCExLow_pi_sysprop);
    //systProps.push_back(FSI_NCasc_FrCExHigh_pi_sysprop); //ETA- NOT USUALLY IN THE OA
    systProps.push_back(FSI_NCasc_FrInelLow_pi_sysprop);
    systProps.push_back(FSI_NCasc_FrPiProd_pi_sysprop);
    systProps.push_back(FSI_NCasc_FrInelHigh_pi_sysprop);
	//ETA CCDIS and CCMPi parameters for 2020OA
    systProps.push_back(BY_DIS_sysprop);
    systProps.push_back(BY_MPi_sysprop);
    systProps.push_back(Xsec_AGKY_sysprop);
	//ETA I-1/2 anti-nu parameter
	systProps.push_back(bgsclbar_sysprop);
	//ETA 2p2h energy dependent parameters
	systProps.push_back(Edep2p2h_lowEnu_sysprop);
	systProps.push_back(Edep2p2h_highEnu_sysprop);
	systProps.push_back(Edep2p2h_lowEnubar_sysprop);
	systProps.push_back(Edep2p2h_highEnubar_sysprop);
	//Q2 parmaeters
	systProps.push_back(LowQ2Suppression_1_sysprop);
	systProps.push_back(LowQ2Suppression_2_sysprop);
	systProps.push_back(LowQ2Suppression_3_sysprop);
	systProps.push_back(LowQ2Suppression_4_sysprop);
	systProps.push_back(LowQ2Suppression_5_sysprop);       
	systProps.push_back(LowQ2Suppression_6_sysprop);       
	systProps.push_back(LowQ2Suppression_7_sysprop);       
	systProps.push_back(LowQ2Suppression_8_sysprop);       

	std::cout << "Size of systProps is " << systProps.size() << std::endl;

    // Call the XSecVary class constructor
    XsecVary2D* xs;
    if (nue){
      xs=new XsecVary2D(wtfile, mtuple, systProps, nutype,"init",var2nue,false,var1nue,isfitqun);
    }
    else{
      xs=new XsecVary2D(wtfile, mtuple, systProps, nutype,"init",var2numu,false,var1numu,isfitqun);
    }
 
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
	xs->GetEntry(0);
	std::cerr<<"Add Kevin's bin hist for 2D as well as 1D"<<std::endl;
	std::cerr<<__FILE__<<":"<<__LINE__<<std::endl;
	throw;
	const int ebins = 2;//true_e_bins->GetNbins();//!!CHANGE THIS LINE
	Double_t * bins{0};//!!
	const Double_t *ebinrange = bins;//true_e_bins->GetXbins()->GetArray();//!!CHANGE THIS LINE
	std::cout << "Number of true energy bins: " << ebins << "." << std::endl;
	std::cout << "True energy bins: {";
	for (int ibin=0; ibin<ebins; ibin++){
	  std::cout << ebinrange[ibin] << ", ";
	}
	std::cout << ebinrange[ebins] << "}" << std::endl;
	
	if(!nue){
	  std::cout << "Number of var1numu bins: " << var1numubins << "." << std::endl;
	  std::cout << "Reco energy bins: {";
	  for (int ibin=0; ibin<var1numubins; ibin++)
	    std::cout << var1numubinrange[ibin] << ", ";
	  std::cout << var1numubinrange[var1numubins] << "}" << std::endl;
	  
	  std::cout << "Number of var2numu bins: " << var2numubins << "." << std::endl;
	  std::cout << "Reco energy bins: {";
	  for (int ibin=0; ibin<var2numubins; ibin++)
	    std::cout << var2numubinrange[ibin] << ", ";
	  std::cout << var2numubinrange[var2numubins] << "}" << std::endl;
	  
	  xs->SetBinning(ebinrange,ebins,var1numubinrange,var1numubins,var2numubinrange,var2numubins);
	
	}
	else{
	  std::cout << "Number of var1nue bins: " << var1nuebins << "." << std::endl;
	  std::cout << "Reco energy bins: {";
	  for (int ibin=0; ibin<var1nuebins; ibin++)
	    std::cout << var1nuebinrange[ibin] << ", ";
	  std::cout << var1nuebinrange[var1nuebins] << "}" << std::endl;
	  
	  std::cout << "Number of var2nue bins: " << var2nuebins << "." << std::endl;
	  std::cout << "Reco energy bins: {";
	  for (int ibin=0; ibin<var2nuebins; ibin++)
	    std::cout << var2nuebinrange[ibin] << ", ";
	  std::cout << var2nuebinrange[var2nuebins] << "}" << std::endl;

	  xs->SetBinning(ebinrange,ebins,var1nuebinrange,var1nuebins,var2nuebinrange,var2nuebins);
	}
	  
 
	
	
	std::cout <<  "Binning implemented successfully" << std::endl; 

	xs->MakeVariations();

	xs->WriteGraphs(outfile);
      }


    return 0;
}

// Print the cmd line syntax
void Usage(){
    cout << "Cmd line syntax should be:" << endl;
    cout << "./bin/make_xsec_response_2d.exe [-t 'byEv' (for event by event splines)] -w wtfile -m mtuple -o outfile [-nue (for nue selection, as opposed to numu)]" << endl;
}

// Messy way to process cmd line arguments.
void ParseArgs(int argc, char **argv){
    //  std::cout << "parse args" << std::endl;
    int nargs = 1; 
    if(argc<(nargs*2+1)){ Usage(); exit(1); }
    for(int i = 1; i < argc; i+=2){
        if(std::string(argv[i]) == "-selec") selection = argv[i+1];
        else if(std::string(argv[i]) == "-w") wtfile = argv[i+1];
        else if(std::string(argv[i]) == "-m") mtuple = argv[i+1];
        else if(std::string(argv[i]) == "-o") outfile = argv[i+1];
	else if(std::string(argv[i]) == "-t") 
	  {
	    if (std::string(argv[i+1])=="byEv") byEvent_splines=true;
	    else byEvent_splines=false;
	  }
	else if(std::string(argv[i]) == "-f") 
	  {
	    if (std::string(argv[i+1])=="true") isfitqun=true;
	    else isfitqun=false;
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
