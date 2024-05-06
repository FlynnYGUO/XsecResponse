// Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_lib/XsecVary2019.cc

#ifndef XsecVary_NDGAr_h
#define XsecVary_NDGAr_h
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
#include <string>
#include "duneanaobj/StandardRecord/StandardRecord.h"

class SystematicProperties_NDGAr
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
    SystematicProperties_NDGAr(std::string name, std::string sname, std::vector<int> modes, int nomLoc);
    SystematicProperties_NDGAr(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add);

    // getters
    const TArrayF * GetWeightArray();
    const TArrayF * GetKnotArray();
    std::vector<int> GetInteractionModes();

    //setters
    void SetArrays(TArrayF* warr, TArrayF* karr);
    void SetBranches(TBranch* wbranch, TBranch* kbranch);
 
};


class XsecVary_NDGAr {
  public :

    // Vector of systematics for which we want splines made, containing info regarding how
    // we want them made
    std::vector<SystematicProperties_NDGAr> systematicProperties;
    std::vector< std::string > sysNames;

	Double_t CalculateErec(int iclass, double LepMom, int Ibound, double LepCosTh);

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
  
    Float_t         Enu;
    Float_t         Pmu;
    Float_t         CosThetamu;
  
    // Declaration of leaf types: fChain
    //Double_t        erec;
    Double_t        wgtflx;
    Double_t        wgtosc;
    //Float_t         wgtosc;
    Float_t         dir[10][3];
    Float_t         fq1rdir[10][7][3];
    Float_t         amome[10];
	//ETA adding in fiqun multi-ring variables needed for numuCC1pi sample(s)
	Float_t         fqmrdir[200][6][3];
	Float_t         fqmrmom[200][6];
	Float_t         fqmreloss[200][6];
	Float_t         fqmomm; //Single ring reco muon mom
	Float_t         fqmome; //Single ring reco electron mom
    Int_t           numnu;
    Int_t           mode;
    Int_t           ipnu[50];
    Double_t         pnu[50];
    
    // List of branches: fChain
    TBranch        *b_iclass;   //!
	TBranch        *b_Ibound;
    TBranch        *b_erec;   //!
    TBranch        *b_mode;   //!
    TBranch        *b_ipnu;   //!
    TBranch        *b_pnu;   //!
    
  // CAF Standard record
    caf::StandardRecord * sr = new caf::StandardRecord();

  //DUNE VARIABLES FD
    Double_t        erec_numu;
    Double_t        erec_nue;
    Double_t        cvnnumu;
    Double_t        cvnnue;
    Int_t           isCC;
    Double_t        vtx_x;
    Double_t        vtx_y;
    Double_t        vtx_z;
    Double_t        berpacv;

    // Declaration of leaf types: fTree_byEv
    Float_t         Etrue;
    Float_t         Erec_check;
  
  
    // Vectors of splines, graphs, and histograms
    std::vector< std::vector< TH2F* > > dev_tmp;
    std::vector< std::string >          dev_tmp_names;
    std::vector< std::vector< std::vector< TSpline3*  >  >  > splines;
    std::vector< std::vector< std::vector< TGraph*  >  >  >   graphs;
  
    // Option for making nue or numu splines (0=nue, 1=numu, any other value indicates an error) [Note: this is currently not used any more - 22/01/2015]
    // Default value (set in constructor if no other value given) is 3
    int SampleType;
  
    XsecVary_NDGAr(std::string weight_file, std::string kinematics_file, std::vector<SystematicProperties_NDGAr> systProps, int sampletype=3); 
    virtual ~XsecVary_NDGAr();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Int_t    GetEntries();
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init();
    virtual void     AddSystematics(std::vector<SystematicProperties_NDGAr> &systProps);
    virtual void     MakeVariations();
    virtual void     WriteGraphs(std::string outputname);
    virtual void     SetBinning(const Double_t *ebins, Int_t nebins, 
         		       const Double_t *rebins, Int_t nrebins);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
  
};

#endif

#ifdef XsecVary_NDGAr_cxx
XsecVary_NDGAr::XsecVary_NDGAr(std::string weight_file, std::string kinematics_file, std::vector<SystematicProperties_NDGAr> systProps, int sampletype)
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

  Init();

  AddSystematics(systematicProperties);

}

XsecVary_NDGAr::~XsecVary_NDGAr()
{
   //if (!fChain_w_NXSec_MaCCQE || !fChain) return;
   if (!fChain_weights.at(0) || !fChain) return;
   if(fChain_weights.at(0)) delete fChain_weights.at(0)->GetCurrentFile();
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
   delete sr;
}

SystematicProperties_NDGAr::SystematicProperties_NDGAr(std::string name, std::string sname, std::vector<int> modes, int nomLoc)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = false;  // default
}

SystematicProperties_NDGAr::SystematicProperties_NDGAr(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = add;
}

void SystematicProperties_NDGAr::SetBranches(TBranch* wbranch, TBranch* kbranch)
{
  weightBranch = wbranch;
  knotBranch   = kbranch;
}

void SystematicProperties_NDGAr::SetArrays(TArrayF* warr, TArrayF* karr)
{
  weightArray = warr;
  knotArray   = karr;
}

const TArrayF * SystematicProperties_NDGAr::GetWeightArray(){return weightArray;}
const TArrayF * SystematicProperties_NDGAr::GetKnotArray(){return knotArray;}

Int_t XsecVary_NDGAr::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t XsecVary_NDGAr::GetEntries()
{
  // Return number of entries in tree
  if (!fChain) return 0;
  return fChain->GetEntries();
}

void XsecVary_NDGAr::AddSystematics(std::vector<SystematicProperties_NDGAr> &systProps)
{
  // Fill vectors before setting branch addresses, because the vectors will reallocate memory
  // while being filled. We need to know the final address to read in the weight file
  int numSysts = systProps.size();
  for(int i = 0; i < numSysts; i++)
  {

    std::string sysName      = systProps[i].systName;
    std::string sysTreeName  = sysName + "_tree";

    std::cout << "sysTreeName = "<<sysTreeName << std::endl;

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

  std::cout<<"numSplines: "<< numSplines<<std::endl;
  // Set up the histogram vector:
  dev_tmp.resize(numSplines);
  std::cout<<"get Entry before"<<std::endl;
  fChain->GetEntry(0);
  std::cout<<"get Entry after"<<std::endl;
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

Long64_t XsecVary_NDGAr::LoadTree(Long64_t entry)
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


void XsecVary_NDGAr::Init()
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
  std::cout<<"XsecVary Init()"<<std::endl;
  fCurrent = -1;
  fChain->SetMakeClass(0);

  fChain->SetBranchAddress("rec", &sr);
/*  
  // DUNE FD CAF Variables
  fChain->SetBranchAddress("cvnnumu", &cvnnumu);
  fChain->SetBranchAddress("cvnnue", &cvnnue);
  fChain->SetBranchAddress("Ev_reco_numu", &erec_numu);
  fChain->SetBranchAddress("Ev_reco_nue", &erec_nue);
  fChain->SetBranchAddress("nuPDG", ipnu, &b_ipnu);
  fChain->SetBranchAddress("Ev", pnu, &b_pnu);
  fChain->SetBranchAddress("mode", &mode, &b_mode);
  fChain->SetBranchAddress("wgtosc", &wgtosc);
  fChain->SetBranchAddress("isCC", &isCC);
  fChain->SetBranchAddress("vtx_x", &vtx_x);
  fChain->SetBranchAddress("vtx_y", &vtx_y);
  fChain->SetBranchAddress("vtx_z", &vtx_z);
  fChain->SetBranchAddress("BeRPA_A_cvwgt", &berpacv);
*/
  Notify();
}

Bool_t XsecVary_NDGAr::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void XsecVary_NDGAr::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

//DUNE FD FV cut
inline bool IsInFDFV(double pos_x_cm, double pos_y_cm, double pos_z_cm) 
{
  return (abs(pos_x_cm) < 310 && abs(pos_y_cm) < 550 && pos_z_cm > 50 && pos_z_cm < 1244);
}


#endif // #ifdef XsecVary_NDGAr_cx
