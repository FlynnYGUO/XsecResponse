// Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_lib/XsecVary20192D.cc

#define XsecVary2D_cxx
#include "XsecVary2D.h"
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>

// This function does most of the work towards creating the xsec splines
void XsecVary2D::MakeVariations()
{

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


    if(cvnnue>0.85){erec=erec_nue;}
    else {erec=erec_numu;}
    if((iclass == 11 || iclass ==13) && (abs(Erec-erec) > 1E-5)){std::cout << "ERROR::Calculated Erec of " << erec << " but in minituples it's " << Erec << std::endl;}

    double firstpar=0;
    if(erec_1d&&!p_1d)
    {
      firstpar=erec/1000.; // MeV
    }
    else
    {
      std::cerr << "parameter is not o!!" << std::endl;
      exit(1);
    }

    double secondpar = 0;
    if(theta_2d && !Q2_2d && !N_2d)
    {
      std::cerr<<"Second parameter not implemented yet"<<std::endl;
      std::cerr<<__FILE__<<":"<<__LINE__<<std::endl;
      throw;
    }
    else
    {
      std::cerr << "2nd parameter is neither Q2 or Theta!!" << std::endl;
      exit(1);
    }

    //Fill up the histograms with weights
    
    // Determine interaction mode
    Int_t modee = -1;

    if(abs(mode)==1) modee = 0; //CCQE
    if(abs(mode)>=11&&abs(mode)<=13) modee = 1; //CC1pi
    if(abs(mode)==16) modee = 2; //CCCOH
    if( (abs(mode)>=17&&abs(mode)<30 && abs(mode) != 21 && abs(mode) != 26) || abs(mode) == 15) modee = 3; //CCMisc ETA: making changes for 2020OA where we have splines that specifically effect CCDIS and CCMpi
    if(abs(mode)==31||abs(mode)==32) modee = 4; //NCpi0
    if(abs(mode)==33||abs(mode)==34) modee = 5; //NCpi+/-
    if(abs(mode)==36) modee = 6; //NCCOH
    if((abs(mode)>=37&&abs(mode)<=52) || abs(mode)==35) modee = 7; //NCOTH includes NCEL
    if(abs(mode)==2) modee = 8; //MEC
    if(abs(mode)==38 || abs(mode)==39) modee = 9; //NC1gamma
	if(abs(mode) == 21){modee = 10;} //CCMpi
	if(abs(mode) == 26){modee = 11;} //CCDIS
   
    if(modee == -1) std::cout << "unknown mode" << std::endl;

    int sysMode = 0;

    // loop over all systematics:
    for(unsigned i = 0; i < systematicProperties.size(); i++)
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
            //dev_tmp[sysMode+k][j]->Fill(pnu[0],firstpar,secondpar,wgtflx*wgtosc*systematicProperties[i].GetWeightArray()->At(j));
            dev_tmp[sysMode+k][j]->Fill(pnu[0],firstpar,secondpar,wgtflx*wgtosc*systematicProperties[i].GetWeightArray()->At(j));
	    if(TMath::IsNaN(wgtflx*wgtosc*systematicProperties[i].GetWeightArray()->At(j))){
	      std::cout << "weight" << wgtflx << ", " << wgtosc << ", " << systematicProperties[i].GetWeightArray()->At(j) << ", pnu[0]: " << pnu[0] << ", firstpar: " << firstpar << ", secondpar: " << secondpar << std::endl;
            }

	    double www = wgtflx*wgtosc*systematicProperties[i].GetWeightArray()->At(j);
		//double www = wgtflx*systematicProperties[i].GetWeightArray()->At(j);
	    if(TMath::IsNaN(www) || www < 0){
	      std::cout << "weight " << www << ", wgtflux " << wgtflx << ", wgtosc " << wgtosc << ", xsec_weight " << systematicProperties[i].GetWeightArray()->At(j) << ", pnu[0]: " << pnu[0] << ", firstpar: " << firstpar << ", secondpar: " << secondpar << ", mode: " << mode << std::endl;
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
              for(int w = 0; w < dev_tmp[sysMode+k][y]->GetNbinsZ(); w++)
              {
 
                int nom = systematicProperties[i].nominalPosition;
                //int size = systematicProperties[i].GetKnotArray()->GetSize();
                graphs[sysMode+k][z][x][w]->SetPoint(iter,(systematicProperties[i].GetKnotArray()->At(0)+2*y*systematicProperties[i].GetKnotArray()->At(nom-1)),
                                                   (dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,w+1) > 0. ?
                                                   dev_tmp[sysMode+k][0]->GetBinContent(z+1,x+1,w+1)/
                                                   dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,w+1) : 1.) );


		if(TMath::IsNaN((dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1, w+1))) || dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1, w+1) < 0){
		std::cout << "ERROR 2D x: is nan at nom: " << std::endl;
	      }
	      
		if(TMath::IsNaN(dev_tmp[sysMode+k][0]->GetBinContent(z+1,x+1,w+1)) || dev_tmp[sysMode+k][0]->GetBinContent(z+1,x+1,w+1)< 0){
		std::cout << "ERROR 2D x : is nan at 0: " << std::endl;
	      }

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
            for(int w = 0; w < dev_tmp[sysMode+k][a]->GetNbinsZ(); w++)
            {
              int nom = systematicProperties[i].nominalPosition;
              graphs[sysMode+k][l][j][w]->SetPoint(iter,systematicProperties[i].GetKnotArray()->At(a),
                                                 (dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1,w+1) > 0. ?
                                                  dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1,w+1)/
                                                  dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1,w+1)  :  1.));

	      if(TMath::IsNaN((dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1, w+1))) || dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1, w+1) < 0){
		std::cout << "ERROR 2D: is nan at nom: " << std::endl;
	      }
	      
	      if(TMath::IsNaN(dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1,w+1)) || dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1,w+1)< 0){
			std::cout << "ERROR 2D: is nan at knot a: " << a <<  std::endl;
			std::cout << "Systematic is " << i << std::endl; 
			std::cout << "Are you returning a NaN? " << TMath::IsNaN(dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1,w+1)) << std::endl;
			std::cout << "Do you have a negative bin content? " << dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1,w+1) << std::endl;
			std::cout << "currently working on sysMode " << sysMode << std::endl;

	      }


            }
          }
        }
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
              for(int w = 0; w < dev_tmp[sysMode+k][y]->GetNbinsZ(); w++)
              {
                int nom = systematicProperties[i].nominalPosition;
                int size = systematicProperties[i].GetKnotArray()->GetSize();
                graphs[sysMode+k][z][x][w]->SetPoint(iter,systematicProperties[i].GetKnotArray()->At(size-1)+y*systematicProperties[i].GetKnotArray()->At(nom+1),
                                                   (dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,w+1) > 0. ?
                                                   dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1,w+1)/
                                                   dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1,w+1) : 1.) );

		if(TMath::IsNaN((dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1, w+1))) || dev_tmp[sysMode+k][nom]->GetBinContent(z+1,x+1, w+1) < 0){
		  std::cout << "ERROR 2D x: is nan at nom: " << std::endl;
		}
		
		if(TMath::IsNaN(dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1,w+1)) || dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1,w+1) < 0){
		  std::cout << "ERROR 2D x: is nan at a: " << size-1 << std::endl;
		  std::cout << "Are you returning a NaN? " << TMath::IsNaN(dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1,w+1)) << std::endl;
		  std::cout << "Do you have a negative bin content? " << dev_tmp[sysMode+k][size-1]->GetBinContent(z+1,x+1,w+1) << 0 << std::endl;
		  std::cout << "currently working on sysMode " << sysMode << std::endl;
		}

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
      std::vector< std::vector< TSpline3* > > tmp_xbin;
      for(int c = 0; c < dev_tmp[a][0]->GetNbinsY(); c++)
      {
        std::vector<TSpline3*> tmp_ybin;
        for(int d = 0; d < dev_tmp[a][0]->GetNbinsZ(); d++)
        {
          char sname[1000];
          sprintf(sname,"%s_sp_%d_%d_%d",dev_tmp_names[a].c_str(),b,c,d);
          TSpline3 *tmp_zbin = new TSpline3(sname,graphs[a][b][c][d]);
          tmp_ybin.push_back(tmp_zbin);
        }
        tmp_xbin.push_back(tmp_ybin);
      }
      if(splines.size() == 0) splines.resize(dev_tmp.size());
      splines[a].push_back(tmp_xbin);
    }
  }

  std::cout << "Finished making splines!" << std::endl;   
 
}

// Etrue-Erec binning (used for splines)
void XsecVary2D::SetBinning(const Double_t *ebins, Int_t nebins, const Double_t *rebins, Int_t nrebins,const Double_t* qsqdbins, Int_t nqsqdbins)
{

  for(unsigned i=0; i<dev_tmp.size(); i++)
  {
    char hname[1000];
    for(unsigned j=0; j<dev_tmp[i].size(); j++)
    {
      sprintf(hname,"dev_tmp_%d_%d",j,i);
      dev_tmp[i][j] = new TH3F(hname,hname,nebins,ebins,nrebins,rebins,nqsqdbins,qsqdbins);
    } 
  }

  for(unsigned a = 0; a < dev_tmp.size(); a++)
  {
    for(int i = 0; i<nebins; i++)
    {
      std::vector< std::vector< TGraph* > > tmp_ebin;
      for(int j = 0; j<nrebins; j++)
      {
        std::vector< TGraph* > tmp_pbin;
        for(int k = 0; k < nqsqdbins; k++)
        {
	  TGraph* tmp_qsqdbin = new TGraph();
	  tmp_pbin.push_back(tmp_qsqdbin);
        }
        tmp_ebin.push_back(tmp_pbin);
      } 
      if(graphs.size() == 0) graphs.resize(dev_tmp.size());
      graphs[a].push_back(tmp_ebin);
    } 
  }

}

// ------------------------------------------------------------------------------------ //

void XsecVary2D::WriteGraphs(char *outputname){

  TFile fout(outputname,"RECREATE");
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

  std::string modeToName[] = {"ccqe","cc1pi","cccoh","ccmisc","ncpiz","ncpipm","nccoh","ncoth","mec", "nc1gamma", "ccmpi", "ccdis"};//ETA adding ccmpi and ccdis for the 2020OA
  int rwbin = 0;
  
  for(unsigned a = 0; a < systematicProperties.size(); a++)
  {
    for(unsigned b = 0; b < systematicProperties[a].intModes.size(); b++)
    {
      for(unsigned i = 0; i < graphs[rwbin+b].size(); i++)
      {
        for(unsigned j = 0; j < (unsigned)graphs[rwbin+b][i].size(); j++)
        {
          for(unsigned k = 0; k < (unsigned)graphs[rwbin+b][i][j].size(); k++)
          {
            char grname[50];
            sprintf(grname,"dev_%s_%s_gr_%d_%d_%d",systematicProperties[a].shortName.c_str(),modeToName[systematicProperties[a].intModes[b]].c_str(),i,j,k);
            graphs[rwbin+b][i][j][k]->Write(grname);
          }
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
        for(unsigned j = 0; j < (unsigned)splines[rwbin+b][i].size(); j++)
        {
          for(unsigned k = 0; k < (unsigned)splines[rwbin+b][i][j].size(); k++)
          {
            char spname[50];
            sprintf(spname,"dev_%s_%s_sp_%d_%d_%d",systematicProperties[a].shortName.c_str(),modeToName[systematicProperties[a].intModes[b]].c_str(),i,j,k);
            splines[rwbin+b][i][j][k]->Write(spname);
          }
        }
      }
    }
    rwbin += systematicProperties[a].intModes.size();
  }

  fout.Close();

  std::cout << "File closed." << std::endl;

}
  
double XsecVary2D::CalcQ2(double Enu,double Elep,double Plep,double costheta,double Mlep)
{
  // Enu      --> neutrino reco energy
  // Elep     --> lepton reco energy
  // Plep     --> lepton reco momentum
  // costheta --> angle between neutrino and lepton

  // Q2 in the laboratory reference system (Giunti section 5.3.1)
  double Q2 = 2*Enu*(Elep-Plep*costheta) - Mlep*Mlep;

  return Q2;
}

double XsecVary2D::CalcXSecTerm(double Enu,double Elep,double Plep,double costheta,double Mlep,double Mn)
{    
  // Term of CCQE dsigma/dQ^2 which sign is different for nu/nubar (arXiv:1305.7513v1, eq. 57-61)
  // Enu      --> neutrino reco energy
  // Elep     --> lepton reco energy
  // Plep     --> lepton reco momentum
  // costheta --> angle between neutrino and lepton
  // Mn       --> nucleon mass

  double MA = 1.2257; // Axial Mass (BANFF post-fit)
  double gA = 1.2694; // FA(Q2=0)

  double Q2 = this->CalcQ2(Enu,Elep,Plep,costheta,Mlep);

  double F1 = 0.5;
  double F2 = 0.5;
  double FA = gA / TMath::Power( 1 + Q2/(MA*MA), 2. );
  double B = Q2/(Mn*Mn)*FA*(F1+F2);

  double XSecTerm = ((4*Mn*Enu - Q2 - Mlep*Mlep) / (Mn*Mn)) * B;

  return XSecTerm;
}
