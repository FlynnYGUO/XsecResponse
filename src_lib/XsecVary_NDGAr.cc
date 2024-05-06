// :Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_lib/XsecVary2019.cc

#define XsecVary_NDGAr_cxx
#include "XsecVary_NDGAr.h"
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>

//Note that inside T2KReWeight the BeRPA parameters are called ABCDU, but they are then saved as ABDEU

// This function does most of the work towards creating the xsec splines
void XsecVary_NDGAr::MakeVariations()
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

    cvnnumu = (double)(sr->common.ixn.dlp[0].nuhyp.cvn.numu);
    cvnnue = (double)(sr->common.ixn.dlp[0].nuhyp.cvn.nue);
    isCC = (int)(sr->mc.nu[0].iscc);
//    vtx_x;
//    vtx_y;
//    vtx_z;
    berpacv = 1;
/*    if (cvnnue > 0.85) {erec = erec_nue;}
    else {erec = erec_numu;}
*/
    erec = (double)(sr->mc.nu[0].E);
//    erec = (double)(sr->common.ixn.dlp[0].Enu.calo);   
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
    //wgtosc = 1;
    

    // Determine interaction mode
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
    
    if(modee == -999) std::cout << "unknown mode" << abs(mode) << std::endl;

	if(TMath::IsNaN(erec) || erec < 0){
	  std::cout << "erec" << erec << std::endl;
	} 

    int sysMode = 0;

    // loop over all systematics:
    for(int i = 0; (unsigned)i <  systematicProperties.size(); i++)
    {
      // loop over all interactions modes this systematic affects
      for(unsigned k = 0; k < systematicProperties[i].intModes.size(); k++)
      {
        // is this event one of these modes?
        //if(systematicProperties[i].intModes[k] == modee)
        //{
          // Does this event pass the FD Fiducial Volume cut?
         // if(IsInFDFV(vtx_x, vtx_y, vtx_z))
         // {
            // loop over all knots
            for(int j = 0 ; j < systematicProperties[i].GetWeightArray()->GetSize(); j++)
            {
            // fill the histogram
	      dev_tmp[sysMode+k][j]->Fill(pnu[0],erec,wgtflx*berpacv*systematicProperties[i].GetWeightArray()->At(j));
            
              // Check which modes have reponses for each systematic
              std::string loopname;
              char lloopname[100];
              sprintf(lloopname,"%s",systematicProperties[i].shortName.c_str());
              loopname = lloopname;
              if((systematicProperties[i].GetWeightArray()->At(j) != 1) & (systematicProperties[i].shortName.c_str() ==systematicProperties[12].shortName.c_str()  || loopname == "empty"|| loopname == "empty")) {std::cout << "Event: " << jentry << " || for Systematic: " << systematicProperties[i].shortName.c_str() <<  "  ||  Mode: "  << modee  << " || CC/NC: " << isCC << " || The knot value is: " << systematicProperties[i].GetWeightArray()->At(j) << std::endl;}
            

              //std::cout << "Osc Weight = " << wgtosc << std::endl;
	      if(TMath::IsNaN(systematicProperties[i].GetWeightArray()->At(j))){
	        std::cout << "Spooky weight: " << systematicProperties[i].GetWeightArray()->At(j) << std::endl;
	      } 

	    double www = wgtflx*systematicProperties[i].GetWeightArray()->At(j);
		//std::cout << "www is " << www << std::endl;
		//std::cout << "wgtflux " << wgtflx << ", xsec_weight " << systematicProperties[i].GetWeightArray()->At(j) << ", pnu[0]: " << pnu[0] << ", Erec: " << erec << ", mode: " << mode << std::endl;
	    if(TMath::IsNaN(www) || www < 0){
		 std::cout << "Sys: "<<  systematicProperties[i].shortName.c_str() <<  "  weight " << www << ", wgtflux " << wgtflx << ", wgtosc " << wgtosc << ", xsec_weight " << systematicProperties[i].GetWeightArray()->At(j) << ", pnu[0]: " << pnu[0] << ", Erec: " << erec << ", mode: " << mode << std::endl;
	    }
          }
       // }
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
                                                dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1)  
                                                :  1.));
            if(i==30 && k == 0 && dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1)!=0 && l == 16 && j == 16) {std::cout <<  " Sys Name: " << systematicProperties[30].shortName.c_str() << "Value of Nom at Knot " << a << " for bins: " << l+1  << " , "  << j+1  << " is: " << dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1) << "  Value of Den: " << dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1) << std::endl;}
            if(TMath::IsNaN((dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1))) || dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1) < 0){
	      std::cout << "ERROR: is nan at nom: " << nom << std::endl;
	    }

	    if(TMath::IsNaN(dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1)) || dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1) < 0){
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
//	TSpline3 *tmp_ybin = NULL;
        std::cout<<"sname: "<<sname<<" graphs[a][b][c]: "<<graphs[a][b][c]<<std::endl;
        TSpline3* tmp_ybin = new TSpline3(sname,graphs[a][b][c]);
        tmp_xbin.push_back(tmp_ybin);
//        delete tmp_ybin;
      }

      if(splines.size() == 0) splines.resize(dev_tmp.size());
      splines[a].push_back(tmp_xbin);
    }
  }

  std::cout << "Finished making splines! Exit." << std::endl;

}

// Etrue-Erec binning (used for splines)
void XsecVary_NDGAr::SetBinning(const Double_t *ebins, Int_t nebins, const Double_t *rebins, Int_t nrebins)
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

void XsecVary_NDGAr::WriteGraphs(std::string outputname){

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
          //graphs[rwbin+b][i][j]->Write(grname);
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
          char spname[100];
          sprintf(spname,"dev_%s_%s_sp_%d_%d",systematicProperties[a].shortName.c_str(),modeToName[systematicProperties[a].intModes[b] -1].c_str(),i,j);
          splines[rwbin+b][i][j]->Write(spname);
		  //if(systematicProperties[a].intModes[b] == 10){
			//std::cout << spname << std::endl;
		  //}
        }
      }
    }
    rwbin += systematicProperties[a].intModes.size();
  }

  fout.Close();

  std::cout << "File closed." << std::endl;

}
