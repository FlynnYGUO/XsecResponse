// Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_lib/XsecVary20192D.cc

#ifndef XsecVary2D_NDGAr_h
#define XsecVary2D_NDGAr_h
#include <stdlib.h>
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TArrayF.h>
#include <TH3F.h>
#include <TH2F.h>
#include <vector>
#include <TGraph.h>
#include <TSpline.h>
#include <TMath.h>
#include "duneanaobj/StandardRecord/StandardRecord.h"


class SystematicProperties2D_NDGAr
{
  public:

    // name of this systematic
    std::string systName;

    // short hand name for this systematic that will be featured in spline titles
    std::string shortName;

    // which interaction modes does this systematic affect?
    std::vector<int> intModes;

    // pointers to the arrays we will read in from the weights file
    TArrayF* weightArray;
    TArrayF* knotArray;

    // which element # in the knot array correspondes to this systematic's nominal value?
    int nominalPosition;

    // do we need to include extra knots outside of the range given by the weight file
    // in order to keep splines flat(ish) outside this range? e.g. pfo in 2017/18 OA
    bool addKnots;

    TBranch* weightBranch;
    TBranch* knotBranch;

    // constructor
    SystematicProperties2D_NDGAr(std::string name, std::string sname, std::vector<int> modes, int nomLoc);
    SystematicProperties2D_NDGAr(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add);

    // getters
    const TArrayF * GetWeightArray();
    const TArrayF * GetKnotArray();
    std::vector<int> GetInteractionModes();

    //setters
    void SetArrays(TArrayF* warr, TArrayF* karr);
    void SetBranches(TBranch* wbranch, TBranch* kbranch);
 
};

class XsecVary2D_NDGAr {
  public :

    // Vector of systematics for which we want splines made, containing info regarding how
    // we want them made
    std::vector<SystematicProperties2D_NDGAr> systematicProperties;
    std::vector< std::string > sysNames;

  
    // Pointer to the analyzed TTree or TChain
    TChain          *fChain;
    Int_t           fCurrent;

    // File to save output tree 
    TFile           *fout;
  
    // Input file containing weights, properly formatted
    std::string               weightFile;
    std::vector< TChain* >    fChain_weights;
    Int_t                     fCurrent_w;
    std::vector< TArrayF* >   weightArrays;
    std::vector< TArrayF* >   knotArrays;
    std::vector< TBranch* >   weightBranches;
    std::vector< TBranch* >   knotBranches;
  
    bool isfitqun;
    bool erec_1d;
    bool p_1d;
    bool Q2_2d;
    bool N_2d;
    bool theta_2d;
    bool is_rhc;
   
    Float_t         Enu;
    Float_t         Pmu;
    Float_t         CosThetamu; 
  
    // Declaration of leaf types: fChain
    Double_t        wgtflx;
    Double_t        wgtosc;
    //Float_t        wgtflx;
    //Float_t        wgtosc;

    Float_t         dir[10][3];
    Float_t         fq1rdir[10][10][3];
    Float_t         amome[10];
    Float_t         amomm[10];
    Float_t         fqmome;
    Float_t         fqmomm;
	//ETA adding in fiqun multi-ring variables needed for numuCC1pi sample(s)
	Float_t         fqmrdir[200][6][3];
	Float_t         fqmrmom[200][6];
	Float_t         fqmreloss[200][6];
    Int_t           mode;
    Int_t           ipnu[50];
    Float_t         pnu[50];
    Double_t        erec;
    Double_t        erec_numu;
    Double_t        erec_nue;
   
    Float_t         Etrue;
  
    // List of branches: fChain
    TBranch        *b_iclass;   //!
	TBranch        *b_Ibound;
    TBranch        *b_erec;   //!
    TBranch        *b_mode;   //!
    TBranch        *b_ipnu;   //!
    TBranch        *b_pnu;   //!

  // CAF Standard record
    caf::StandardRecord * sr = new caf::StandardRecord();
 
 //DUNE VARIABLES ND
    Double_t        erec_caf;
    Double_t        elep_reco;
    Int_t           reco_numu;
    Int_t           reco_nue;
    Int_t           muon_contained;
    Int_t           muon_tracker;
	Double_t        Ehad_veto;
    Int_t           isCC;
    Double_t        vtx_x;
    Double_t        vtx_y;
    Double_t        vtx_z;
	Double_t        berpacv;


    // Vectors of splines, graphs, and histograms
    std::vector< std::vector< TH3F* > > dev_tmp;
    std::vector< std::string >          dev_tmp_names;
    std::vector< std::vector< std::vector< std::vector< TSpline3* > > > > splines;
    std::vector< std::vector< std::vector< std::vector< TGraph* > > > >   graphs;
    
    TSpline3* tmp_zbin;  
    // Option for making nue or numu splines (0=nue, 1=numu, any other value indicates an error) [Note: this is currently not used any more - 22/01/2015]
    // Default value (set in constructor if no other value given) is 3
    int SampleType;
   
    XsecVary2D_NDGAr(std::string weight_file, std::string kinematics_file, std::vector<SystematicProperties2D_NDGAr> systProps, int sampletype=3); 
    virtual ~XsecVary2D_NDGAr();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Int_t    GetEntries();
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init();
    virtual void     AddSystematics(std::vector<SystematicProperties2D_NDGAr> &systProps);
    virtual void     MakeVariations();
    virtual void     WriteGraphs(std::string outputname);
    virtual void     SetBinning(const Double_t *ebins, Int_t nebins, 
       		       const Double_t *rebins, Int_t nrebins, const Double_t *qsqdbins, Int_t nqsqdbins);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    virtual double   CalcQ2(double Enu, double Elep, double Plep, double costheta, double Mlep);
    virtual double   CalcXSecTerm(double Enu, double Elep, double Plep, double costheta, double Mlep, double Mn);
    
};

#endif

#ifdef XsecVary2D_NDGAr_cxx
XsecVary2D_NDGAr::XsecVary2D_NDGAr(std::string weight_file, std::string kinematics_file, std::vector<SystematicProperties2D_NDGAr> systProps , int sampletype)
{ 
  std::cout << "XsecVary2D constructor with sampletype = " << sampletype << std::endl;
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

  Init();

  AddSystematics(systematicProperties);

}


XsecVary2D_NDGAr::~XsecVary2D_NDGAr()
{
  if (!fChain_weights.at(0) || !fChain) return;
  if(fChain_weights.at(0)) delete fChain_weights.at(0)->GetCurrentFile();
  if(fChain) delete fChain->GetCurrentFile();
/*
  std::cout << "Deleting splines " << std::endl;
  for(unsigned i = 0; i < splines.size(); i++)
  {
    for(unsigned j = 0; j < splines[i].size(); j++)
    {
      for(unsigned k = 0; k < splines[i][j].size(); k++)
      {
        for(unsigned l = 0; l < splines[i][j][k].size(); l++)
        {
          if(splines[i][j][k][l]) splines[i][j][k][l]->Delete();
        }
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
        for(unsigned l = 0; l < graphs[i][j][k].size(); l++)
        {
          if(graphs[i][j][k][l]) graphs[i][j][k][l]->Delete();
        }
      }
    }
  }
*/
  std::cout << "Deleting histograms " << std::endl;
  for(unsigned i = 0; i < dev_tmp.size(); i++)
  {
    for(unsigned j = 0; j < dev_tmp[i].size(); j++)
    {
      if(dev_tmp[i][j]) dev_tmp[i][j]->Delete();
    }
  }

}

SystematicProperties2D_NDGAr::SystematicProperties2D_NDGAr(std::string name, std::string sname, std::vector<int> modes, int nomLoc)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = false;  // default
}

SystematicProperties2D_NDGAr::SystematicProperties2D_NDGAr(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = add;
}

void SystematicProperties2D_NDGAr::SetBranches(TBranch* wbranch, TBranch* kbranch)
{
  weightBranch = wbranch;
  knotBranch   = kbranch;
}

void SystematicProperties2D_NDGAr::SetArrays(TArrayF* warr, TArrayF* karr)
{
  weightArray = warr;
  knotArray   = karr;
}

const TArrayF * SystematicProperties2D_NDGAr::GetWeightArray(){return weightArray;}
const TArrayF * SystematicProperties2D_NDGAr::GetKnotArray(){return knotArray;}

Int_t XsecVary2D_NDGAr::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Int_t XsecVary2D_NDGAr::GetEntries()
{
  // Return number of entries in tree
  if (!fChain) return 0;
  return fChain->GetEntries();
}

void XsecVary2D_NDGAr::AddSystematics(std::vector<SystematicProperties2D_NDGAr> &systProps)
{
  // Fill vectors before setting branch addresses, because the vectors will reallocate memory
  // while being filled. We need to know the final address to read in the weight file
  int numSysts = systProps.size();
  for(int i = 0; i < numSysts; i++)
  {
    std::string sysName      = systProps[i].systName;
    std::string sysTreeName  = sysName + "_tree";

    std::cout << sysTreeName << std::endl;

    TArrayF *tempWeightsArr = 0;
    TArrayF *tempKnotsArr = 0;
    TBranch *tempWeightsBranch = 0;
    TBranch *tempKnotsBranch   = 0;
    TChain  *tempChain = new TChain(sysTreeName.c_str());
    tempChain->Add(weightFile.c_str());
    fChain->AddFriend(sysTreeName.c_str(),weightFile.c_str());
  
    sysNames.push_back(sysName);
    weightArrays.push_back(tempWeightsArr);
    knotArrays.push_back(tempKnotsArr);
    weightBranches.push_back(tempWeightsBranch);
    knotBranches.push_back(tempKnotsBranch);
    fChain_weights.push_back(tempChain);

    // Set the arrays to read in the weight file.
    systProps[i].SetArrays(tempWeightsArr, tempKnotsArr);
   
  }

  int numSplines = 0;

  // Now set the branch addresses.
  for(int j = 0; j < numSysts; j++)
  {
    std::string weightBrName = sysNames.at(j) + "_weights";
    std::string knotBrName   = sysNames.at(j) + "_knots";

    systProps[j].weightArray  = 0;
    systProps[j].knotArray    = 0;
    systProps[j].weightBranch = 0;
    systProps[j].knotBranch   = 0;
  
    fChain->SetBranchAddress(weightBrName.c_str(),&(systProps[j].weightArray),&(systProps[j].weightBranch));
    fChain->SetBranchAddress(knotBrName.c_str(),&(systProps[j].knotArray),&(systProps[j].knotBranch));

    numSplines += systProps[j].intModes.size();
  }

  // Set up the histogram vector:
  dev_tmp.resize(numSplines);
  fChain->GetEntry(0);
  unsigned int itty = 0;
  std::string modeToName[] = {"ccqe", "ccmec", "ccdis", "ccres", "cccoh", "ccdiff", "ccnueel", "unknown", "ccamnugamma", "unknown", "cccohel", "ccibd", "ccglasres", "ccimdannihilation", "ncqe", "ncmec", "ncdis", "ncres", "nccoh", "ncdiff", "ncnueel", "ncamnugamma", "nccohel", "ncibd", "ncglasres", "ncimdannihilation"};
  for(int k = 0; k < numSysts; k++)
  {
    for(unsigned l = 0; l < systProps[k].intModes.size(); l++)
    {
      dev_tmp[itty].resize(systProps[k].GetKnotArray()->GetSize());
      itty++;
      char tmpname[1000];
      sprintf(tmpname,"dev_%s_%s",systProps[k].shortName.c_str(),modeToName[systProps[k].intModes[l]-1].c_str());
      dev_tmp_names.push_back(std::string(tmpname));
    }
  }

}

Long64_t XsecVary2D_NDGAr::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent_w)
  { // check this part, currently matches k
    fCurrent = chain->GetTreeNumber(); 
    Notify();
  }
  return centry;
}


void XsecVary2D_NDGAr::Init()
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
  fChain->SetMakeClass(0);

  fChain->SetBranchAddress("rec", &sr);

/*  fChain->SetBranchAddress("mode", &mode, &b_mode);
  fChain->SetBranchAddress("reco_numu", &reco_numu);
  fChain->SetBranchAddress("reco_nue", &reco_nue);
  fChain->SetBranchAddress("Ev_reco", &erec_caf);
  fChain->SetBranchAddress("Elep_reco", &elep_reco);
  fChain->SetBranchAddress("nuPDG", ipnu, &b_ipnu);
  fChain->SetBranchAddress("Ev", pnu, &b_pnu);
  fChain->SetBranchAddress("mode", &mode, &b_mode);
  fChain->SetBranchAddress("isCC", &isCC);
  fChain->SetBranchAddress("vtx_x", &vtx_x);
  fChain->SetBranchAddress("vtx_y", &vtx_y);
  fChain->SetBranchAddress("vtx_z", &vtx_z);
  fChain->SetBranchAddress("BeRPA_A_cvwgt", &berpacv);
  fChain->SetBranchAddress("muon_contained", &muon_contained);
  fChain->SetBranchAddress("muon_tracker", &muon_tracker);
  fChain->SetBranchAddress("Ehad_veto", &Ehad_veto);
*/ 
  Notify();
}

Bool_t XsecVary2D_NDGAr::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

   return kTRUE;
}

void XsecVary2D_NDGAr::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

//DUNE ND CC inclusive Selection Cut
inline bool IsCCInclusive(int reco_numu, int muon_contained, int muon_tracker, double Ehad_veto) 
{
  return (reco_numu && (muon_contained || muon_tracker) && Ehad_veto < 30);
}



//DUNE ND FV cut
inline bool IsInNDFV(double pos_x_cm, double pos_y_cm, double pos_z_cm) 
{
  bool inDeadRegion = false;
  for (int i = -3; i <= 3; ++i) {
    // 0.5cm cathode in the middle of each module, plus 0.5cm buffer
    double cathode_center = i * 102.1;
    if (pos_x_cm > cathode_center - 0.75 && pos_x_cm < cathode_center + 0.75) {
      inDeadRegion = true;
	}

    // 1.6cm dead region between modules (0.5cm module wall and 0.3cm pixel
    // plane, x2) don't worry about outer boundary because events are only
    // generated in active Ar + insides
    double module_boundary = i * 102.1 + 51.05;
    if (i <= 2 && pos_x_cm > module_boundary - 1.3 && pos_x_cm < module_boundary + 1.3) {
       inDeadRegion = true; 
	}
  }
  for (int i = 1; i <= 4; ++i) {
    // module boundaries in z are 1.8cm (0.4cm ArCLight plane + 0.5cm module
    // wall, x2) module is 102.1cm wide, but only 101.8cm long due to cathode
    // (0.5cm) being absent in length but ArCLight is 0.1cm thicker than pixel
    // plane so it's 0.3cm difference positions are off-set by 0.6 because I
    // defined 0 to be the upstream edge based on the active volume by
    // inspecting a plot, and aparently missed by 3 mm, but whatever add 8mm =
    // 2 pad buffer due to worse position resolution in spatial dimension z
    // compared to timing direction x so total FV gap will be 1.8 + 2*0.8
    // = 3.4cm
    double module_boundary = i * 101.8 - 0.6;
    if (pos_z_cm > module_boundary - 1.7 && pos_z_cm < module_boundary + 1.7) {
	  inDeadRegion = true;
    }
  }
    
  return (abs(pos_x_cm) < 200 && abs(pos_y_cm) < 100 && pos_z_cm > 50 &&
		              pos_z_cm < 350 && !inDeadRegion);
}

#endif // #ifdef XsecVary_NDGAr_cxx
