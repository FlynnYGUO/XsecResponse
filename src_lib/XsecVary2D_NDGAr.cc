// Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_lib/XsecVary20192D.cc

#define XsecVary2D_NDGAr_cxx
#include "XsecVary2D_NDGAr.h"
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>

// This function does most of the work towards creating the xsec splines
void XsecVary2D_NDGAr::MakeVariations()
{

  if (fChain == 0 ) return;

  // define a "map" from interaction mode number (which is the element number within the array)
  // to a string describing the interaction type
  
  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry); 
    double firstpar=0;
    firstpar = (double)(sr->mc.nu[0].E);
//    firstpar=erec_caf; // MeV
	
    double secondpar = 0;
    secondpar = (double)((sr->mc.nu[0].E - sr->mc.nu[0].prim[0].p.E)/(sr->mc.nu[0].E));
	//Yrec
//	secondpar = (erec_caf - elep_reco)/erec_caf;
	//std::cout << "Etrue = " << pnu[1] << std::endl;
	//std::cout << "Erec = " << firstpar << std::endl;
	//std::cout << "Yrec = " << secondpar << std::endl;
	//std::cout << "Elep = " << elep_reco << std::endl;

    //Fill up the histograms with weights
    isCC = (int)(sr->mc.nu[0].iscc);
    berpacv = 1;
    mode = sr->mc.nu[0].mode;

    bool simb_mode = false;

	wgtflx = 1;
	wgtosc = 1;
    // Determine interaction mode
    Int_t modee = -999;
	if(simb_mode) {

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

      if(modee == -999) std::cout << "unknown mode: " << abs(mode) << " IsCC: " << isCC << std::endl;
	}

	else {

      if((abs(mode)==-1) && (isCC == 1)) modee = -999; //Unknown
      if((abs(mode)==1) && (isCC == 1)) modee = 1; //CCQE
      if((abs(mode)==3) && (isCC == 1)) modee = 3; //CC DIS
      if((abs(mode)==4) && (isCC == 1)) modee = 4; //CC RES
      if((abs(mode)==5) && (isCC == 1)) modee = 5; //CC COH
      if((abs(mode)==6) && (isCC == 1)) modee = 6; //CC Diffractive
      if((abs(mode)==7) && (isCC == 1)) modee = 7; //CC Electron scattering
      if((abs(mode)==9) && (isCC == 1)) modee = 9; //CC AMnuGamma
      if((abs(mode)==10) && (isCC == 1)) modee = 2; //CC MEC
      if((abs(mode)==11) && (isCC == 1)) modee = 11; //CC COH Elastic
      if((abs(mode)==12) && (isCC == 1)) modee = 12; // CC IBD
      if((abs(mode)==13) && (isCC == 1)) modee = 13; //CC Glashow RES
      if((abs(mode)==14) && (isCC == 1)) modee = 14;//CC IMD Annihalation 
      if((abs(mode)==8) && (isCC == 1)) modee = 14;//CC IMD Annihalation 
   
      if((abs(mode)==1) && (isCC == 0)) modee = 15; //NCQE
      if((abs(mode)==3) && (isCC == 0)) modee = 16; //NC DIS
      if((abs(mode)==4) && (isCC == 0)) modee = 17; //NC RES
      if((abs(mode)==5) && (isCC == 0)) modee = 18; //NC COH
      if((abs(mode)==6) && (isCC == 0)) modee = 19; //NC Diffractive
      if((abs(mode)==7) && (isCC == 0)) modee = 20; //NC Electron scattering
      if((abs(mode)==9) && (isCC == 0)) modee = 21; //NC AMnuGamma
      if((abs(mode)==10) && (isCC == 0)) modee = 22; //NC MEC
      if((abs(mode)==11) && (isCC == 0)) modee = 23; //NC COH Elastic
      if((abs(mode)==12) && (isCC == 0)) modee = 24; // NC IBD
      if((abs(mode)==13) && (isCC == 0)) modee = 25; //NC Glashow RES
      if((abs(mode)==14) && (isCC == 0)) modee = 26;//NC IMD Annihalation 
      if((abs(mode)==8) && (isCC == 0)) modee = 26;//NC IMD Annihalation 
    
      if(modee == -999) std::cout << "unknown mode: " << abs(mode) << std::endl;
	}
    
    int sysMode = 0;

    // loop over all systematics:
    for(unsigned i = 0; i < systematicProperties.size(); i++)
    {
      // loop over all interactions modes this systematic affects
      for(unsigned k = 0; k < systematicProperties[i].intModes.size(); k++)
      {
        // is this event one of these modes and does the event pass the selection cut
//        if(systematicProperties[i].intModes[k] == modee && IsInNDFV(vtx_x, vtx_y, vtx_z) && IsCCInclusive(reco_numu, muon_contained, muon_tracker, Ehad_veto))
//        {
          // loop over all knots
          for(int j = 0 ; j < systematicProperties[i].GetWeightArray()->GetSize(); j++)
          {
            // fill the histogram
            //dev_tmp[sysMode+k][j]->Fill(pnu[1],firstpar,secondpar,wgtflx*wgtosc*systematicProperties[i].GetWeightArray()->At(j));
            dev_tmp[sysMode+k][j]->Fill(pnu[1],firstpar,secondpar,wgtflx*wgtosc*berpacv*systematicProperties[i].GetWeightArray()->At(j));

			//std::cout << "Event = " << jentry << " || Etrue = " << pnu[1] << "|| Erec = " << firstpar << " || Yrec = " << secondpar << " || WEIGHT = " << wgtflx*wgtosc*berpacv*systematicProperties[i].GetWeightArray()->At(j) << std::endl;

	    if(TMath::IsNaN(wgtflx*systematicProperties[i].GetWeightArray()->At(j))){
	      std::cout << "weight" << wgtflx << ", " << wgtosc << ", " << systematicProperties[i].GetWeightArray()->At(j) << ", pnu[1]: " << pnu[1] << ", firstpar: " << firstpar << ", secondpar: " << secondpar << std::endl;
            }

	    double www = wgtflx*wgtosc*systematicProperties[i].GetWeightArray()->At(j);
		//double www = wgtflx*systematicProperties[i].GetWeightArray()->At(j);
	    if(TMath::IsNaN(www) || www < 0){
	      std::cout << "weight " << www << ", wgtflux " << wgtflx << ", wgtosc " << wgtosc << ", xsec_weight " << systematicProperties[i].GetWeightArray()->At(j) << ", pnu[1]: " << pnu[1] << ", firstpar: " << firstpar << ", secondpar: " << secondpar << ", mode: " << mode << std::endl;
	    }


          }
//        }
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
			  std::cout << "Current mode = " << sysMode+k << " || " << " nom bin = " << nom << " || xbin = " << l+1  << " || ybin = " << j+1 << " || zbin = " << w+1 << std::endl;
			  std::cout << "Whats in this? : " << dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1,w+1) << std::endl;
              graphs[sysMode+k][l][j][w]->SetPoint(iter,systematicProperties[i].GetKnotArray()->At(a),
                                                 (dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1,w+1) > 0. ?
                                                  dev_tmp[sysMode+k][a]->GetBinContent(l+1,j+1,w+1)/
                                                  dev_tmp[sysMode+k][nom]->GetBinContent(l+1,j+1,w+1)  :  1.));
              
              std::cout<<"graphs getN(): "<<graphs[sysMode+k][l][j][w]->GetN()<<std::endl;
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
          std::cout<<"a = "<<a<<", b = "<<b<<", c ="<<c<<", d = "<<d<<" max d = "<<dev_tmp[a][0]->GetNbinsZ()<<std::endl;
          sprintf(sname,"%s_sp_%d_%d_%d",dev_tmp_names[a].c_str(),b,c,d);
          std::cout<<"N: "<<graphs[a][b][c][d]->GetN()<<" X mean: "<<graphs[a][b][c][d]->GetMean(1)<<" Y mean: "<<graphs[a][b][c][d]->GetMean(2)<<std::endl;         
          std::cout<<"sname: "<<sname<<std::endl;
          tmp_zbin = new TSpline3(sname,graphs[a][b][c][d]);
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
void XsecVary2D_NDGAr::SetBinning(const Double_t *ebins, Int_t nebins, const Double_t *rebins, Int_t nrebins,const Double_t* qsqdbins, Int_t nqsqdbins)
{
  std::cout<<"SetBinning"<<std::endl;
  for(unsigned i=0; i<dev_tmp.size(); i++)
  {
    char hname[1000];
    for(unsigned j=0; j<dev_tmp[i].size(); j++)
    { 
      sprintf(hname,"dev_tmp_%d_%d",j,i);
      std::cout<<"dev_tmp_"<<j<<"_"<<i<<" nebins: "<<nebins<<" ebins: "<<*ebins<<" nrebins: "<<nrebins<<" rebins: "<<*rebins<<" nqsqdbins: "<<nqsqdbins<< "qsqdbins: "<<*qsqdbins<<std::endl;
      dev_tmp[i][j] = new TH3F(hname,hname,nebins,ebins,nrebins,rebins,nqsqdbins,qsqdbins);
      const double *netrue = dev_tmp[i][j]->GetXaxis()->GetXbins()->GetArray();
      std::cout<<"netrue bins: "<<netrue[0]<<" "<<netrue[4]<<std::endl;
      const double *nerec = dev_tmp[i][j]->GetYaxis()->GetXbins()->GetArray();
      std::cout<<"netrue bins: "<<nerec[0]<<" "<<nerec[4]<<std::endl;
      const double *nyrec = dev_tmp[i][j]->GetZaxis()->GetXbins()->GetArray();
      std::cout<<"netrue bins: "<<nyrec[0]<<" "<<nyrec[4]<<std::endl;
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

void XsecVary2D_NDGAr::WriteGraphs(std::string outputname){

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

	std::string modeToName[] = {"ccqe", "ccmec", "ccdis", "ccres", "cccoh", "ccdiff", "ccnueel", "unknown", "ccamnugamma", "unknown", "cccohel", "ccibd", "ccglasres", "ccimdannihilation", "ncqe", "ncdis", "ncres", "nccoh", "ncdiff", "ncnueel", "ncamnugamma", "ncmec", "nccohel", "ncibd", "ncglasres", "ncimdannihilation"};
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
            sprintf(grname,"dev_%s_%s_gr_%d_%d_%d",systematicProperties[a].shortName.c_str(),modeToName[systematicProperties[a].intModes[b]-1].c_str(),i,j,k);
            graphs[rwbin+b][i][j][k]->Write(grname);
            graphs[rwbin+b][i][j][k]->Delete();
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
            char spname[1000];
            sprintf(spname,"dev_%s_%s_sp_%d_%d_%d",systematicProperties[a].shortName.c_str(),modeToName[systematicProperties[a].intModes[b]-1].c_str(),i,j,k);
            //std::cout<<"a: "<<a<<" b: "<<b<<" i: "<<i<<" j: "<<j<<" k: "<<k<<std::endl;
            //std::cout<<"spname: "<<spname<<std::endl;
            //std::string problem = "dev_thetadelta_ccres_sp_0_5_4";
            //if(spname == problem){std::cout<< "here"<<std::endl; splines[rwbin+b][i+1][j][k]->Write(spname); continue;}
            if(splines[rwbin+b][i][j][k]){ 
            splines[rwbin+b][i][j][k]->Write(spname);
            splines[rwbin+b][i][j][k]->Delete();
            }
            else {std::cout<<"spname: "<<spname<<" NO SPLINE"<<std::endl;}
          }
        }
      }
    }
    rwbin += systematicProperties[a].intModes.size();
  }

  fout.Close();

  std::cout << "File closed." << std::endl;

}
  
double XsecVary2D_NDGAr::CalcQ2(double Enu,double Elep,double Plep,double costheta,double Mlep)
{
  // Enu      --> neutrino reco energy
  // Elep     --> lepton reco energy
  // Plep     --> lepton reco momentum
  // costheta --> angle between neutrino and lepton

  // Q2 in the laboratory reference system (Giunti section 5.3.1)
  double Q2 = 2*Enu*(Elep-Plep*costheta) - Mlep*Mlep;

  return Q2;
}

double XsecVary2D_NDGAr::CalcXSecTerm(double Enu,double Elep,double Plep,double costheta,double Mlep,double Mn)
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
