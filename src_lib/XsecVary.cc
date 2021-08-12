// Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_lib/XsecVary2019.cc

#define XsecVary_cxx
#include "XsecVary.h"
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>

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

    //Get LepMom and LepCosTh
    //erec = CalculateErec(iclass, LepMom, Ibound, LepCosTh);
    if (cvnnue > 0.85) {erec = erec_nue;}
    else {erec = erec_numu;}

        
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

    wgtflx = 1;
    wgtosc = 1;
    

    // Determine interaction mode
    Int_t modee = -999;
    if(abs(mode)==-1) modee = -999; //Unknown
    if(abs(mode)==0) modee = 1; //CCQE
    if(abs(mode)==2) modee = 3; //DIS
    if(abs(mode)==1) modee = 4; //RES
    if(abs(mode)==3) modee = 5; //COH
    if(abs(mode)==11) modee = 6; //Diffractive
    if(abs(mode)==5) modee = 7; //Nu(
    if(abs(mode)==9) modee = 9; //AMnuGamma
    if(abs(mode)==10) modee = 10; //MEC
    if(abs(mode)==4) modee = 11; //COH Elastic
    if(abs(mode)==7) modee = 12; // IBD
    if(abs(mode)==8)modee = 13; //Glashow RES
    if(abs(mode)==6)modee = 14;//IMD Annihalation 
   
    if(modee == -999) std::cout << "unknown mode" << abs(mode) << std::endl;

	if(TMath::IsNaN(erec) || erec < 0){
	  std::cout << "erec" << erec << std::endl;
	} 

    int sysMode = 0;

    // loop over all systematics:
    for(int i = 0; i < 1; i++)
    {
      // loop over all interactions modes this systematic affects
      for(unsigned k = 0; k < systematicProperties[i].intModes.size(); k++)
      {
        // is this event one of these modes?
        if(systematicProperties[i].intModes[k] == modee)
        {
          // loop over all knots
          for(int j = 0 ; j < systematicProperties[i].GetWeightArray()->GetSize(); j++)
          {
            // fill the histogram
            dev_tmp[sysMode+k][j]->Fill(pnu[0],erec,wgtflx*wgtosc*systematicProperties[i].GetWeightArray()->At(j));
            //std::cout << "BinContent: " << dev_tmp[sysMode+k][j]->FindBin(pnu[0], erec) << std::endl;
            if(dev_tmp[sysMode+k][j]->GetBinContent(dev_tmp[sysMode+k][j]->FindBin(pnu[0], erec)) != 1) {std::cout << "Mode: "  << modee  << std::endl;}
            //dev_tmp[sysMode+k][j]->Fill(pnu[0],erec,wgtflx*systematicProperties[i].GetWeightArray()->At(j));
	    if(TMath::IsNaN(pnu[0])){
	      std::cout << "pnu[0]" << pnu[0] << std::endl;
	    } 

	    double www = wgtflx*wgtosc*systematicProperties[i].GetWeightArray()->At(j);
		//std::cout << "erec is " << erec << std::endl;
		//std::cout << "Pnu is " << pnu[0] << std::endl;
		//std::cout << "costheta is " << LepCosTh << std::endl;
	    //double www = wgtflx*systematicProperties[i].GetWeightArray()->At(j);
		//std::cout << "www is " << www << std::endl;
		//std::cout << "wgtflux " << wgtflx << ", xsec_weight " << systematicProperties[i].GetWeightArray()->At(j) << ", pnu[0]: " << pnu[0] << ", Erec: " << erec << ", mode: " << mode << std::endl;
	    if(TMath::IsNaN(www) || www < 0){
		  std::cout << "weight " << www << ", wgtflux " << wgtflx << ", wgtosc " << wgtosc << ", xsec_weight " << systematicProperties[i].GetWeightArray()->At(j) << ", pnu[0]: " << pnu[0] << ", Erec: " << erec << ", mode: " << mode << std::endl;
	    }
          }
        }
      }
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
                                                 (dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1) > 0. ?
                                                 dev_tmp[sysMode+k][0]->GetBinContent(z+1,x+1)/
                                                 dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1) : 1.) );


	      if(TMath::IsNaN((dev_tmp[sysMode+k][nom+1]->GetBinContent(z+1,x+1))) || (dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1)) < 0){
		std::cout << "ERROR x: is nan at nom: " << nom << std::endl;
	      }
	      
	      if(TMath::IsNaN(dev_tmp[sysMode+k][0]->GetBinContent(z+1,x+1)) || dev_tmp[sysMode+k][0]->GetBinContent(z+1,x+1) < 0){
		std::cout << "ERROR x: is nan at 0: " << std::endl;
	      }


            }
          }
          iter++;
        }
      }
      int counter = 0;
      // fill the graphs with at the supplied (from the weight file) knot locations
      for(int a = 0; a < systematicProperties[i].GetKnotArray()->GetSize(); a++)
      {
        for(int l = 0; l < dev_tmp[sysMode+k][a]->GetNbinsX(); l++)
        {
          for(int j = 0; j < dev_tmp[sysMode+k][a]->GetNbinsY(); j++ )
          {
            
            //std::cout << "Bin Number: " << dev_tmp[sysMode+k][a]->GetNbinsY() << std::endl;
            int nom = systematicProperties[i].nominalPosition;
            graphs[sysMode+k][l][j]->SetPoint(iter,systematicProperties[i].GetKnotArray()->At(a),
                                               (dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1) > 0. ?
                                                dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1)/
                                                dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1)  :  1.));
            
            if (dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1) != 0) {std::cout << "Value:"  << dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1) << "|| Bin number : " << l+1 << " , "  << j+1 << std::endl;
            counter++;}

 
            if(TMath::IsNaN((dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1))) || dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1) < 0){
	      std::cout << "ERROR: is nan at nom: " << nom << std::endl;
	    }

	    if(TMath::IsNaN(dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1)) || dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1) < 0){
	      std::cout << "ERROR: is nan at a: " << a << std::endl;
	    }

          }
        }
        if (a == 0) {std::cout << "Total bins with stuff " << counter << std::endl;}
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
                                                 (dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1) > 0. ?
                                                 dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1)/
                                                 dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1) : 1.) );



	      if(TMath::IsNaN((dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1))) || dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1) < 0){
		std::cout << "ERROR x: is nan at nom: " << std::endl;
	      }
	      
	      if(TMath::IsNaN(dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1)) || dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1)< 0){
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

  for(unsigned i=0; i<dev_tmp.size(); i++)
  {
    char hname[30];
    for(unsigned j=0; j<dev_tmp[i].size(); j++)
    {
      sprintf(hname,"dev_tmp_%d_%d",j,i);
      dev_tmp[i][j] = new TH2F(hname,hname,nebins,ebins,nrebins,rebins);
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

  std::string modeToName[] = {"ccqe","cc1pi","cccoh","ccmisc","ncpiz","ncpipm","nccoh","ncoth","mec", "nc1gamma", "ccmpi", "ccdis"};//ETA adding ccmpi and ccdis for the 2020OA
  for(unsigned a = 0; a < systematicProperties.size(); a++)
  {
    for(unsigned b = 0; b < systematicProperties[a].intModes.size(); b++)
    {
      for(unsigned i = 0; i < graphs[rwbin+b].size(); i++)
      {
        for(unsigned j = 0; j < (unsigned)graphs[rwbin+b][i].size(); j++)
        {
          char grname[50];
          sprintf(grname,"dev_%s_%s_gr_%d_%d",systematicProperties[a].shortName.c_str(),modeToName[systematicProperties[a].intModes[b]].c_str(),i,j);
          graphs[rwbin+b][i][j]->Write(grname);
        }
      }
    }
    rwbin += systematicProperties[a].intModes.size();
  }

  rwbin = 0;

  std::cout << "Writing splines." << std::endl;
  for(unsigned a = 0; a < systematicProperties.size(); a++)
  {
    for(unsigned b = 0; b < systematicProperties[a].intModes.size(); b++)
    {
      for(unsigned i = 0; i < splines[rwbin+b].size(); i++)
      {
        for(unsigned j = 0; j < splines[rwbin+b][i].size(); j++)
        {
          char spname[50];
          sprintf(spname,"dev_%s_%s_sp_%d_%d",systematicProperties[a].shortName.c_str(),modeToName[systematicProperties[a].intModes[b]].c_str(),i,j);
          splines[rwbin+b][i][j]->Write(spname);
        }
      }
    }
    rwbin += systematicProperties[a].intModes.size();
  }

  fout.Close();

  std::cout << "File closed." << std::endl;

}
