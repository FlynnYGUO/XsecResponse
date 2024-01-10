// :Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_lib/XsecVary2019.cc

#define XsecVary_cxx
#include "XsecVary.h"
#include "SRGlobal.h"
#include <TH3.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>

//DUNE FD FV cut
inline bool IsInFDFV(double pos_x_cm, double pos_y_cm, double pos_z_cm) 
{
  return true;
  return (abs(pos_x_cm) < 310 && abs(pos_y_cm) < 550 && pos_z_cm > 50 && pos_z_cm < 1244);
}

inline Int_t get_modee(Int_t mode, Int_t isCC){
  Int_t modee = -999;
  if((abs(mode)==-1) && (isCC == 1)) modee = -999; //Unknown
  if((abs(mode)==0) && (isCC == 1)) modee = 1; //CCQE
  if((abs(mode)==2) && (isCC == 1)) modee = 3; //CC DIS
  if((abs(mode)==1) && (isCC == 1)) modee = 4; //CC RES
  if((abs(mode)==3) && (isCC == 1)) modee = 5; //CC COH
  if((abs(mode)==11) && (isCC == 1)) modee = 6; //CC Diffractive
  if((abs(mode)==5) && (isCC == 1)) modee = 7; //CC Electron scattering
  if((abs(mode)==9) && (isCC == 1)) modee = 9; //CC AMnuGamma
  if((abs(mode)==10) && (isCC == 1)) modee = 2; //CC MEC
  if((abs(mode)==4) && (isCC == 1)) modee = 11; //CC COH Elastic
  if((abs(mode)==7) && (isCC == 1)) modee = 12; // CC IBD
  if((abs(mode)==8) && (isCC == 1)) modee = 13; //CC Glashow RES
  if((abs(mode)==6) && (isCC == 1)) modee = 14;//CC IMD Annihalation 
  
  if((abs(mode)==0) && (isCC == 0)) modee = 15; //NCQE
  if((abs(mode)==2) && (isCC == 0)) modee = 16; //NC DIS
  if((abs(mode)==1) && (isCC == 0)) modee = 17; //NC RES
  if((abs(mode)==3) && (isCC == 0)) modee = 18; //NC COH
  if((abs(mode)==11) && (isCC == 0)) modee = 19; //NC Diffractive
  if((abs(mode)==5) && (isCC == 0)) modee = 20; //NC Electron scattering
  if((abs(mode)==9) && (isCC == 0)) modee = 21; //NC AMnuGamma
  if((abs(mode)==10) && (isCC == 0)) modee = 22; //NC MEC
  if((abs(mode)==4) && (isCC == 0)) modee = 23; //NC COH Elastic
  if((abs(mode)==7) && (isCC == 0)) modee = 24; // NC IBD
  if((abs(mode)==8) && (isCC == 0)) modee = 25; //NC Glashow RES
  if((abs(mode)==6) && (isCC == 0)) modee = 26;//NC IMD Annihalation

  return modee;
}

//Note that inside T2KReWeight the BeRPA parameters are called ABCDU, but they are then saved as ABDEU

// This function does most of the work towards creating the xsec splines
void XsecVary::MakeVariations()
{

  //ETA - moving erec definition out of header, now get this from reco information
  //removes any ambiguity as to how this was calculated when making minituples
  double erec = -999; 
  
  if (fChain == 0 ) return;

  // Check whether making variations for 1Rnue or 1Rnumu events
  switch(SampleType)
  {
    case 0: // nue
      std::cout << "Making variations for Nue splines" << std::endl;
      std::cout << "With the current code, this doesn't actually mean anything!" << std::endl;
      break;
    case 1: // numu
      std::cout << "Making variations for Numu splines" << std::endl;
      std::cout << "With the current code, this doesn't actually mean anything!" << std::endl;
      break;
    default: // Something's wrong, SampleType not set
      std::cerr << "Error: SampleType not set correctly" << std::endl
       	 << "Please use SampleType = 0 for nue, or SampleType = 1 for numu" << std::endl;
      std::cout << "With the current code, this doesn't actually mean anything!" << std::endl;
  }

  // define a "map" from interaction mode number (which is the element number within the array)
  // to a string describing the interaction type

  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break; 
    fChain->GetEntry(jentry); 
    if (cvnnue > 0.85) {erec = erec_nue;}
    else {erec = erec_numu;}
    // std::cout << "erec: " << erec << std::endl;
        
	/*if(abs(erec - Erec_check) > 1E-6){
	  std::cout << "~~~~~~~~~~~~~~~~~~~~ " << std::endl;
	  std::cout << "Erec in minituple is " << Erec_check << std::endl;
	  std::cout << "Erec is " << erec << std::endl;
	  std::cout << "Mode is " << mode << std::endl;
	  std::cout << "iclass is " << iclass << std::endl;
	  std::cout << "LepMom is " << LepMom << std::endl;
	  std::cout << "dirnu is: " << dirnu[0][0] << ", " << dirnu[0][1] << ", " << dirnu[0][2] << std::endl;
	  std::cout << "fq1rdir is: " << fq1rdir[0][1][0] << ", " << fq1rdir[0][1][1] << ", " << fq1rdir[0][1][2] << std::endl;
	  std::cout << "LepCosTh is " << LepCosTh << std::endl;	
	}	
	*/

    //Osc and fluc weight: Set to 1

    // wgtflx = 1;
    wgtosc = 1;
    

    // Determine interaction mode
    Int_t modee = get_modee(mode, isCC);
    if(modee == -999) std::cout << "unknown mode" << abs(mode) << std::endl;

	if(TMath::IsNaN(erec) || erec < 0){
	  std::cout << "erec" << erec << std::endl;
	} 

    int sysMode = 0;
    berpacv = 1;

    // loop over all systematics:
    for(int i = 0; (unsigned)i <  systematicProperties.size(); i++)
    {
      // loop over all interactions modes this systematic affects
      for(unsigned k = 0; k < systematicProperties[i].intModes.size(); k++)
      {
        // is this event one of these modes?
        if(systematicProperties[i].intModes[k] == modee)
        {
          // Does this event pass the FD Fiducial Volume cut?
          std::vector<float> weightArr = weightVec->at(systematicProperties[i].syst_idx);
          if(IsInFDFV(vtx_x, vtx_y, vtx_z))
          {
            // loop over all knots
            for(int j = 0 ; j < systematicProperties[i].GetKnotArray()->GetSize(); j++)
            {

              //Sanity checks
              if(TMath::IsNaN(weightArr[j])){
                std::cout << "Spooky weight: " << weightArr[j] << std::endl;
                continue;
              }

              float corr_weight = weightArr[j];

              if(weightArr[j] < 0){
                corr_weight = 0;
              }

              
              if(systematicProperties[i].reweight && systematicProperties[i].nominalPosition != j && abs(corr_weight - 1) > 1e-5){
                corr_weight /= wgtxsec;
              }

              // fill the histogram
              dev_tmp[sysMode+k][j]->Fill(pnu,erec,-0.5,wgtflx*berpacv*corr_weight);
              // std::cout << "Filling histo -> " << pnu << " ; " << erec << " ; " << wgtflx << " x " << berpacv << " x " << corr_weight << std::endl;
            
              // Check which modes have reponses for each systematic
              std::string loopname;
              char lloopname[100];
              sprintf(lloopname,"%s",systematicProperties[i].shortName.c_str());
              loopname = lloopname;
              // if((weightArr[j] != 1) & (systematicProperties[i].shortName.c_str() ==systematicProperties[12].shortName.c_str()  || loopname == "empty"|| loopname == "empty")) {std::cout << "Event: " << jentry << " || for Systematic: " << systematicProperties[i].shortName.c_str() <<  "  ||  Mode: "  << modee  << " || CC/NC: " << isCC << " || The knot value is: " << systematicProperties[i].GetWeightArray()->At(j) << std::endl;}
            

              //std::cout << "Osc Weight = " << wgtosc << std::endl;
	      

	    // double www = wgtflx*weightArr[j];
		//std::cout << "www is " << www << std::endl;
		//std::cout << "wgtflux " << wgtflx << ", xsec_weight " << systematicProperties[i].GetWeightArray()->At(j) << ", pnu[0]: " << pnu[0] << ", Erec: " << erec << ", mode: " << mode << std::endl;
	  //   if(TMath::IsNaN(www) || www < 0){
		//  std::cout << "Sys: "<<  systematicProperties[i].shortName.c_str() <<  "  weight " << www << ", wgtflux " << wgtflx << ", wgtosc " << wgtosc << ", xsec_weight " << weightArr[j] << ", pnu[0]: " << pnu[0] << ", Erec: " << erec << ", mode: " << mode << std::endl;
	    }
          }
        }
      }
    //}
 
      sysMode += systematicProperties[i].intModes.size();
    }
  } // event loop

  int sysMode = 0;
  
  // loop over all systematics:
  for(unsigned i = 0; i < systematicProperties.size(); i++)
  {
    // loop over all interactions modes this systematic affects
    for(unsigned k = 0; k < systematicProperties[i].intModes.size(); k++)
    {
      int iter = 0;
      // if this systematic requires extra knots to be added at the edges, do it
      if(systematicProperties[i].addKnots)
      {
        for(int y=1; y<3; y++)
        {
          for(int z = 0; z < dev_tmp[sysMode+k][y]->GetNbinsX(); z++)
          {
            for(int x = 0; x < dev_tmp[sysMode+k][y]->GetNbinsY(); x++)
            {
              int nom = systematicProperties[i].nominalPosition;
              //int size = systematicProperties[i].GetKnotArray()->GetSize();
              graphs[sysMode+k][z][x]->SetPoint(iter,(systematicProperties[i].GetKnotArray()->At(0)+2*y*systematicProperties[i].GetKnotArray()->At(nom-1)),
                                                 (dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,1) > 0. ?
                                                 dev_tmp[sysMode+k][0]->GetBinContent(z+1,x+1,1)/
                                                 dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,1) : 1.) );


	      if(TMath::IsNaN((dev_tmp[sysMode+k][nom+1]->GetBinContent(z+1,x+1,1))) || (dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,1)) < 0){
		std::cout << "ERROR x: is nan at nom: " << nom << std::endl;
	      }
	      
	      if(TMath::IsNaN(dev_tmp[sysMode+k][0]->GetBinContent(z+1,x+1,1)) || dev_tmp[sysMode+k][0]->GetBinContent(z+1,x+1,1) < 0){
		std::cout << "ERROR x: is nan at 0: " << std::endl;
	      }


            }
          }
          iter++;
        }
      }
      
      // fill the graphs with at the supplied (from the weight file) knot locations
      for(int a = 0; a < systematicProperties[i].GetKnotArray()->GetSize(); a++)
      {
        for(int l = 0; l < dev_tmp[sysMode+k][a]->GetNbinsX(); l++)
        {
          for(int j = 0; j < dev_tmp[sysMode+k][a]->GetNbinsY(); j++ )
          {
            
            //std::cout << "Bin Number: " << dev_tmp[sysMode+k][a]->GetNbinsY() << std::endl;
            int nom = systematicProperties[i].nominalPosition;
            Float_t fill_value = 1.;

            // if(l == dev_tmp[sysMode+k][a]->GetNbinsX() - 1 && j == dev_tmp[sysMode+k][a]->GetNbinsY() - 1 && a != nom){
            if(a != nom){
              fill_value = 1.00002; //Ensuring to have no syst with only flat splines
            }

            if(dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1,1) > 0.){
              fill_value = dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1,1)/dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1,1);
            }                                   
                  
            graphs[sysMode+k][l][j]->SetPoint(iter,systematicProperties[i].GetKnotArray()->At(a), fill_value);
            if(abs(fill_value) > 100){
              std::cout << "Filling graph at bin : " << sysMode + k << " ; " << l << " ; " << j << " ; " << a << " with value " << fill_value << std::endl;
            }
            if(i==30 && k == 0 && dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1,1)!=0 && l == 16 && j == 16) {std::cout <<  " Sys Name: " << systematicProperties[30].shortName.c_str() << "Value of Nom at Knot " << a << " for bins: " << l+1  << " , "  << j+1  << " is: " << dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1,1) << "  Value of Den: " << dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1) << std::endl;}
            if(TMath::IsNaN((dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1,1))) || dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1,1) < 0){
	      std::cout << "ERROR: is nan at nom: " << nom << std::endl;
	    }

	    if(TMath::IsNaN(dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1,1)) || dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1,1) < 0){
	      std::cout << "ERROR: is nan at a: " << a << std::endl;
	    }

          }
        }
        if (a == systematicProperties[i].nominalPosition ) {std::cout << "For Sys: " << systematicProperties[i].shortName.c_str() << "  for event mode: " << k <<  "  the nom bin with the most stuff = " << dev_tmp[sysMode+k][a]->GetBinContent(dev_tmp[sysMode+k][a]->GetMaximumBin())/dev_tmp[sysMode+k][a]->GetBinContent(dev_tmp[sysMode+k][a]->GetMaximumBin()) << std::endl;}
        iter++;
      }

      // if this systematic requires extra knots to be added at the edges, do it
      if(systematicProperties[i].addKnots)
      {
        for(int y=1; y<3; y++)
        {
          for(int z = 0; z < dev_tmp[sysMode+k][y]->GetNbinsX(); z++)
          {
            for(int x = 0; x < dev_tmp[sysMode+k][y]->GetNbinsY(); x++)
            {
              int nom = systematicProperties[i].nominalPosition;
              int size = systematicProperties[i].GetKnotArray()->GetSize();
              graphs[sysMode+k][z][x]->SetPoint(iter,systematicProperties[i].GetKnotArray()->At(size-1)+y*systematicProperties[i].GetKnotArray()->At(nom+1),
                                                 (dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,1) > 0. ?
                                                 dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1,1)/
                                                 dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,1) : 1.) );



	      if(TMath::IsNaN((dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,1))) || dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,1) < 0){
		std::cout << "ERROR x: is nan at nom: " << std::endl;
	      }
	      
	      if(TMath::IsNaN(dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1,1)) || dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1,1)< 0){
		std::cout << "ERROR x: is nan at 0: " << std::endl;
	      }



            }
          }
          iter++;
        }
      }
    }
    sysMode += systematicProperties[i].intModes.size();
  }

  // ---
  std::cout << "Filling splines..." << std::endl;
  //std::cout << "Size check" << dev_tmp[35][0]->GetNbinsX() << std::endl;
  // ---

  for(unsigned a = 0; a < dev_tmp.size(); a++)
  {
    for(int b = 0; b < dev_tmp[a][0]->GetNbinsX(); b++)
    {
      std::vector<TSpline3*> tmp_xbin;
      for(int c = 0; c < dev_tmp[a][0]->GetNbinsY(); c++)
      {
        char sname[1000];
        sprintf(sname,"%s_sp_%d_%d",dev_tmp_names[a].c_str(),b,c);
	      TSpline3 *tmp_ybin = NULL;
        tmp_ybin = new TSpline3(sname,graphs[a][b][c]);
        tmp_xbin.push_back(tmp_ybin);
      }

      if(splines.size() == 0) splines.resize(dev_tmp.size());
      splines[a].push_back(tmp_xbin);
    }
  }

  std::cout << "Finished making splines! Exit." << std::endl;

}

// Etrue-Erec binning (used for splines)
void XsecVary::SetBinning(const Double_t *ebins, Int_t nebins, const Double_t *rebins, Int_t nrebins)
{
  double angBins[3] = {-1, 0, 1};
  for(unsigned i=0; i<dev_tmp.size(); i++)
  {
    char hname[30];
    for(unsigned j=0; j<dev_tmp[i].size(); j++)
    {
      sprintf(hname,"dev_tmp_%d_%d",j,i);
      dev_tmp[i][j] = new TH3F(hname,hname,nebins,ebins,nrebins,rebins, 2, angBins);
    } 
  }

  for(unsigned a = 0; a < dev_tmp.size(); a++)
  {
    for(int i = 0; i<nebins; i++)
    {
      std::vector<TGraph*> tmp_ebin; 
      for(int j = 0; j<nrebins; j++)
      {
        TGraph* tmp_pbin = new TGraph();
        tmp_ebin.push_back(tmp_pbin);
      }
      if(graphs.size() == 0) graphs.resize(dev_tmp.size());
      graphs[a].push_back(tmp_ebin);
    }
  }

}

// ------------------------------------------------------------------------------------ //

void XsecVary::WriteGraphs(std::string outputname){

  TFile fout(outputname.c_str(),"RECREATE");
  fout.cd();

  for(unsigned i=0; i<dev_tmp.size(); i++)
  {
    char hname[30];
    for(unsigned j=0; j<dev_tmp[i].size(); j++)
    {
      sprintf(hname,"dev_tmp_%d_%d",j,i);
      dev_tmp[i][j]->Write(hname);
    } 
  }

  int rwbin = 0;

  //std::string modeToName[] = {"ccqe","ccqe","cccoh","ccmisc","ncpiz","ncpipm","nccoh","ncoth","mec", "nc1gamma", "ccmpi", "ccdis"};//ETA adding ccmpi and ccdis for the 2020OA
	std::string modeToName[] = {"ccqe", "ccmec", "ccdis", "ccres", "cccoh", "ccdiff", "ccnueel", "unknown", "ccamnugamma", "unknown", "cccohel", "ccibd", "ccglasres", "ccimdannihilation", "ncqe", "ncdis", "ncres", "nccoh", "ncdiff", "ncnueel", "ncamnugamma", "ncmec", "nccohel", "ncibd", "ncglasres", "ncimdannihilation"};

  for(unsigned a = 0; a < systematicProperties.size(); a++)
  {
    for(unsigned b = 0; b < systematicProperties[a].intModes.size(); b++)
    {
      for(unsigned i = 0; i < graphs[rwbin+b].size(); i++)
      {
        for(unsigned j = 0; j < (unsigned)graphs[rwbin+b][i].size(); j++)
        {
          char grname[100];
          //std::cout << "Index is " << systematicProperties[a].intModes[b]-1 << std::endl;
          sprintf(grname,"dev_%s_%s_gr_%d_%d",systematicProperties[a].shortName.c_str(),modeToName[systematicProperties[a].intModes[b]-1].c_str(),i,j);
		  //ETA - some problem with writing large files so trying removing writing TGraphs (think it might be a problem with the size of the file)
          // graphs[rwbin+b][i][j]->Write(grname);
        }
      }
    }
    rwbin += systematicProperties[a].intModes.size();
  }

  rwbin = 0;

  std::cout << "Writing splines." << std::endl;
  int nb = 0;
  for(unsigned a = 0; a < systematicProperties.size(); a++)
  {
    std::cout << "Systematic " << systematicProperties[a].shortName.c_str() << std::endl;
    for(unsigned b = 0; b < systematicProperties[a].intModes.size(); b++)
    {
      // char spname2[100];
      // sprintf(spname2,"dev_%s_%s_%d",systematicProperties[a].shortName.c_str(),modeToName[systematicProperties[a].intModes[b] -1].c_str(), systematicProperties[a].intModes[b] -1);
      // std::cout << spname2 << std::endl;
      // continue;
      for(unsigned i = 0; i < splines[rwbin+b].size(); i++)
      {
        for(unsigned j = 0; j < splines[rwbin+b][i].size(); j++)
        {
          for(uint k=0; k < 2; k++){
            char spname[100];
            nb++;
            sprintf(spname,"dev_%s_%s_sp_%d_%d_%d",systematicProperties[a].shortName.c_str(),modeToName[systematicProperties[a].intModes[b] -1].c_str(),i,j,k);
            // if(splines[rwbin+b][i][j]->Eval(2) != 1)
            //   std::cout << "Spline value at 2: " << splines[rwbin+b][i][j]->Eval(2) << std::endl;
            splines[rwbin+b][i][j]->Write(spname);
          //   if(a == 0 && b == 0){
          //     Double_t x, y, s, t, u;
          //     for(uint i = 0; i < 7u; i++){
          //       splines[rwbin+b][i][j]->GetCoeff(i, x, y, s, t, u);
          //       std::cout << systematicProperties[a].shortName.c_str() << " : " << modeToName[systematicProperties[a].intModes[b] -1] << " -> " <<  x << " " << y << " " << s << " " << t << " " << u << std::endl;
          //     }
          //   }
          }
          
		  //if(systematicProperties[a].intModes[b] == 10){
			//std::cout << spname << std::endl;
		  //}
        }
      }
    }
    rwbin += systematicProperties[a].intModes.size();
  }

  fout.Close();

  std::cout << "File closed. Wrote " << nb << " splines" << std::endl;

}

XsecVary::XsecVary(std::string sum_file, std::vector<SystematicProperties> systProps, int sampletype)
: XsecVary(sum_file, sum_file, systProps, sampletype) {}

XsecVary::XsecVary(std::string weight_file, std::string kinematics_file, std::vector<SystematicProperties> systProps, int sampletype)
{ 
  std::cout << "XsecVary constructor with sampletype = " << sampletype << std::endl;
  // Set SampleType (should be 0 for nue or 1 for numu)
  SampleType = sampletype;

  systematicProperties = systProps;

  weightFile = weight_file;

  fChain = new TChain("cafTree");//h1/mtuple/minituple
  if(fChain)   fChain->Add(kinematics_file.c_str());
  if(fChain->GetEntries()==0)
  {
    std::cout << "Looking for tree 'h1' instead." << std::endl;
    fChain = new TChain("h1");
    if(fChain)   fChain->Add(kinematics_file.c_str());
  }
  if(fChain->GetEntries()==0)
  {
    std::cout << "Looking for tree 'minituple' instead." << std::endl;
    fChain = new TChain("minituple");
    if(fChain)   fChain->Add(kinematics_file.c_str());
  }

  weightVec = nullptr;
  Init();
  // abort();
  AddSystematics(systematicProperties);

}

XsecVary::~XsecVary()
{
   //if (!fChain_w_NXSec_MaCCQE || !fChain) return;
  //  if (!fChain_weights.at(0) || !fChain) return;
  //  if(fChain_weights.at(0)) delete fChain_weights.at(0)->GetCurrentFile();
   if(fChain) delete fChain->GetCurrentFile();

   std::cout << "Deleting splines " << std::endl;
   for(unsigned i = 0; i < splines.size(); i++)
   {
     for(unsigned j = 0; j < splines[i].size(); j++)
     {
       for(unsigned k = 0; k < splines[i][j].size(); k++)
       {
         if(splines[i][j][k]) splines[i][j][k]->Delete();
       }
     }
   }

   std::cout << "Deleting graphs " << std::endl;
   for(unsigned i = 0; i < graphs.size(); i++)
   {
     for(unsigned j = 0; j < graphs[i].size(); j++)
     {
       for(unsigned k = 0; k < graphs[i][j].size(); k++)
       {
         if(graphs[i][j][k]) graphs[i][j][k]->Delete();
       }
     }
   }

   std::cout << "Deleting histograms " << std::endl;
   for(unsigned i = 0; i < dev_tmp.size(); i++)
   {
     for(unsigned j = 0; j < dev_tmp[i].size(); j++)
     {
       if(dev_tmp[i][j]) dev_tmp[i][j]->Delete();
     }
   }

}

SystematicProperties::SystematicProperties(std::string name, std::string sname, std::vector<int> modes, int nomLoc)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = false;  // default
  reweight          = false;  // default
}

SystematicProperties::SystematicProperties(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = add;
  reweight          = false;
  std::cout << "USING THE ADDKNOT CONSTRUCTOR" << std::endl;
}

SystematicProperties::SystematicProperties(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add, bool _reweight)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = add;
  reweight          = _reweight;
  // reweight          = false;
  std::cout << "USING THE REWEIGHT CONSTRUCTOR" << std::endl;
}

void SystematicProperties::SetBranches(TBranch* wbranch, TBranch* kbranch)
{
  weightBranch = wbranch;
  knotBranch   = kbranch;
}

void SystematicProperties::SetArrays(TArrayF* warr, TArrayF* karr)
{
  weightArray = warr;
  knotArray   = karr;
}

const TArrayF * SystematicProperties::GetWeightArray(){return weightArray;}
const TArrayF * SystematicProperties::GetKnotArray(){return knotArray;}

Int_t XsecVary::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t XsecVary::GetEntries()
{
  // Return number of entries in tree
  if (!fChain) return 0;
  return fChain->GetEntries();
}

void XsecVary::AddSystematics(std::vector<SystematicProperties> &systProps)
{

  //Get Systematics from file
  TChain  *tempChain = new TChain("globalTree");
  tempChain->Add(weightFile.c_str());
  // tempChain->Add("caf_1.root");
  tempChain->Print();
  caf::SRGlobal *global = nullptr;
  std::vector<std::vector<double>> *shifts = nullptr;
  tempChain->SetBranchAddress("global", &global);
  tempChain->SetBranchAddress("shifts", &shifts);

  tempChain->GetEntry(0);


  // Fill vectors before setting branch addresses, because the vectors will reallocate memory
  // while being filled. We need to know the final address to read in the weight file
  int numSysts = systProps.size();
  for(int i = 0; i < numSysts; i++)
  {

    std::string sysName      = systProps[i].systName;
    std::string sysTreeName  = sysName + "_tree";

   // std::cout << sysTreeName << std::endl;

    TArrayF *tempWeightsArr = 0;
    TArrayF *tempKnotsArr = 0;
    TBranch *tempWeightsBranch = 0;
    TBranch *tempKnotsBranch   = 0;
    // TChain  *tempChain = new TChain(sysTreeName.c_str());
    // tempChain->Add(weightFile.c_str());
    // fChain->AddFriend(sysTreeName.c_str(),weightFile.c_str());
    sysNames.push_back(sysName);
    weightArrays.push_back(tempWeightsArr);
    knotArrays.push_back(tempKnotsArr);
    weightBranches.push_back(tempWeightsBranch);
    knotBranches.push_back(tempKnotsBranch);
    // fChain_weights.push_back(tempChain);

    // Set the arrays to read in the weight file.
    systProps[i].SetArrays(tempWeightsArr, tempKnotsArr);
   
  }


  int numSplines = 0;

  std::vector<caf::SRSystParamHeader> &names = global->wgts.params;

  // Now set the branch addresses.
  for(int j = 0; j < numSysts; j++)
  {
    std::string weightBrName = sysNames.at(j) + "_weights";
    std::string knotBrName   = sysNames.at(j) + "_knots";

    systProps[j].weightArray  = 0;
    systProps[j].knotArray    = 0;
    // systProps[j].weightBranch = 0;
    // systProps[j].knotBranch   = 0;

    //Getting the matching if between IFILE systs and code systs
    int syst_idx = -1;

    for(uint k=0; k<names.size(); k++){
      if(names[k].name == sysNames.at(j)){
        syst_idx = k;
        break;
      }
    }

    if(syst_idx == -1){
      std::cout << "ERROR: Could not find syst " << sysNames.at(j) << " in the input file" << std::endl;
      abort();
    }

    systProps[j].syst_idx = syst_idx;
    std::cout << sysNames.at(j) << " -> " << syst_idx << std::endl;

    //Filling the knots
    std::vector<float> knots = std::vector<float>(shifts->at(syst_idx).begin(), shifts->at(syst_idx).end());
    systProps[j].knotArray = new TArrayF(knots.size(), knots.data());
    std::cout << "Knots for syst " << sysNames.at(j) << " -> ";
    for(const float& knot : knots){
      std::cout << knot << " ";
    }
    std::cout << std::endl;


    // fChain->SetBranchAddress(weightBrName.c_str(),&(systProps[j].weightArray),&(systProps[j].weightBranch));
    // fChain->SetBranchAddress(knotBrName.c_str(),&(systProps[j].knotArray),&(systProps[j].knotBranch));
  
    numSplines += systProps[j].intModes.size();

  }


  // Set up the histogram vector:
  dev_tmp.resize(numSplines);
  unsigned int itty = 0;
  //std::string modeToName[] = {"ccqe","cc1pi","cccoh","ccmisc","ncpiz","ncpipm","nccoh","ncoth","mec", "nc1gamma", "ccmpi", "ccdis"};//ETA adding ccmpi and ccdis for 2020OA, ccoth now ccmisc also was nc1gamma missing??
  std::string modeToName[] = {"ccqe", "ccmec", "ccdis", "ccres", "cccoh", "ccdiff", "ccnueel", "unknown", "ccamnugamma", "unknown", "cccohel", "ccibd", "ccglasres", "ccimdannihilation", "ncqe", "ncmec", "ncdis", "ncres", "nccoh", "ncdiff", "ncnueel", "ncamnugamma", "nccohel", "ncibd", "ncglasres", "ncimdannihilation"};
  for(int k = 0; k < numSysts; k++)
  {
    for(unsigned l = 0; l < systProps[k].intModes.size(); l++)
    {
      //std::cout << systProps[k].weightArray << std::endl;
      dev_tmp[itty].resize(systProps[k].GetKnotArray()->GetSize());
      itty++;
      char tmpname[1000];
      sprintf(tmpname,"dev_%s_%s",systProps[k].shortName.c_str(),modeToName[systProps[k].intModes[l]-1].c_str());
      dev_tmp_names.push_back(std::string(tmpname));
    }
  }
}

Long64_t XsecVary::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent_w) { // check this part, currently matches k
     //if (chain->GetTreeNumber() != fCurrent) {
     fCurrent = chain->GetTreeNumber(); 
     Notify();
   }
   return centry;
}


void XsecVary::Init()
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Note that for the for the weight chain, this is done in AddSystematics()

  // Set branch addresses and branch pointers
  fCurrent = -1;
  // fChain->SetMakeClass(1);


  // DUNE FD CAF Variables
  fChain->SetBranchAddress("cvnnumu", &cvnnumu);
  fChain->SetBranchAddress("cvnnue", &cvnnue);
  fChain->SetBranchAddress("Ev_reco_numu", &erec_numu);
  fChain->SetBranchAddress("Ev_reco_nue", &erec_nue);
  fChain->SetBranchAddress("nuPDG", &ipnu, &b_ipnu);
  fChain->SetBranchAddress("Ev", &pnu, &b_pnu);
  fChain->SetBranchAddress("mode", &mode, &b_mode);
  // fChain->SetBranchAddress("wgtosc", &wgtosc);
  fChain->SetBranchAddress("weight", &wgtflx);
  fChain->SetBranchAddress("genie_weight", &wgtxsec);
  fChain->SetBranchAddress("isCC", &isCC);
  fChain->SetBranchAddress("vtx_x", &vtx_x);
  fChain->SetBranchAddress("vtx_y", &vtx_y);
  fChain->SetBranchAddress("vtx_z", &vtx_z);
  fChain->SetBranchAddress("xsSyst_wgt", &weightVec);
  // fChain->SetBranchAddress("cvwgt", &weightVec2);
  // fChain->SetBranchAddress("BeRPA_A_cvwgt", &berpacv);
  Notify();
}

Bool_t XsecVary::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void XsecVary::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}