// Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_lib/XsecVary2019.cc

#ifndef XsecVary_h
#define XsecVary_h
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

class SystematicProperties
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
    SystematicProperties(std::string name, std::string sname, std::vector<int> modes, int nomLoc);
    SystematicProperties(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add);

    // getters
    const TArrayF * GetWeightArray();
    const TArrayF * GetKnotArray();
    std::vector<int> GetInteractionModes();

    //setters
    void SetArrays(TArrayF* warr, TArrayF* karr);
    void SetBranches(TBranch* wbranch, TBranch* kbranch);
 
};


class XsecVary {
  public :

    // Vector of systematics for which we want splines made, containing info regarding how
    // we want them made
    std::vector<SystematicProperties> systematicProperties;
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
    Int_t           iclass;
	Int_t           Ibound;
    //Double_t        erec;
    Double_t        pizmass;
    Double_t        pizmom;
    Double_t        pizcosb;
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
	Int_t           ipp; //index of 2R pi+pi+ fit
    Int_t           npar;
    Float_t         wallv;
    Float_t         ipv[50];
    Float_t         posv[3];
    Float_t         dirv[50][3];
    Float_t         pmomv[50];
    Int_t           npar2;
    Int_t           ipv2[50];
    Int_t           numnu;
    Int_t           mode;
    Int_t           ipnu[50];
    Double_t         pnu[50];
    Float_t         dirnu[50][3];
    Int_t           Npvc;
    Int_t           Ipvc[100];
    Int_t           Ichvc[100];
    Int_t           Iorgvc[100];
    Int_t           Iflvc[100];
    Float_t         Abspvc[100];
    Float_t         Pvc[100][3];
    Float_t         Crsx;
    Float_t         Crsy;
    Float_t         Crsz;
    Float_t         Crsphi;
    Int_t           Nvert;
    Float_t         Posvert[300][3];
    Int_t           Iflgvert[300];
    Int_t           Nvcvert;
    Float_t         Dirvert[900][3];
    Float_t         Abspvert[900];
    Float_t         Abstpvert[900];
    Int_t           Ipvert[900];
    Int_t           Iverti[900];
    Int_t           Ivertf[900];
    Float_t         Fsiprob;
    Float_t         posc[2][6];
    Int_t           EventType;
    Double_t        FlxWeight;
    Double_t        OscWeight;
    Float_t         pi0mass[2];
    Double_t        erec_numu;
    Double_t        erec_nue;
    Double_t        cvnnumu;
    Double_t        cvnnue;
    // Declaration of leaf types: fTree_byEv
    Float_t         Etrue;
    Float_t         Erec_check;
  
    // List of branches: fChain
    TBranch        *b_iclass;   //!
	TBranch        *b_Ibound;
	TBranch        *b_fqmomm;
	TBranch        *b_fqmome;
	TBranch        *b_fqmrdir;
	TBranch        *b_fq1rdir;
	TBranch        *b_fqmrmom;
	TBranch        *b_fqmreloss;
	//ETA - need ipp i.e. the index of the 2Rpipi fit for the numuCC1pi samples
	TBranch        *b_ipp;
    TBranch        *b_erec;   //!
    TBranch        *b_pizmass;   //!
    TBranch        *b_pizmom;   //!
    TBranch        *b_pizcosb;   //!
    TBranch        *b_wgtflx;   //!
    TBranch        *b_wgtosc;   //!
    TBranch        *b_dir;   //!
    TBranch        *b_amome;   //!
    TBranch        *b_npar;   //!
    TBranch        *b_wallv;   //!
    TBranch        *b_ipv;   //!
    TBranch        *b_posv;   //!
    TBranch        *b_dirv;   //!
    TBranch        *b_pmomv;   //!
    TBranch        *b_npar2;   //!
    TBranch        *b_ipv2;   //!
    TBranch        *b_numnu;   //!
    TBranch        *b_mode;   //!
    TBranch        *b_ipnu;   //!
    TBranch        *b_pnu;   //!
    TBranch        *b_dirnu;   //!
    TBranch        *b_Npvc;   //!
    TBranch        *b_Ipvc;   //!
    TBranch        *b_Ichvc;   //!
    TBranch        *b_Iorgvc;   //!
    TBranch        *b_Iflvc;   //!
    TBranch        *b_Abspvc;   //!
    TBranch        *b_Pvc;   //!
    TBranch        *b_Crsx;   //!
    TBranch        *b_Crsy;   //!
    TBranch        *b_Crsz;   //!
    TBranch        *b_Crsphi;   //!
    TBranch        *b_Nvert;   //!
    TBranch        *b_Posvert;   //!
    TBranch        *b_Iflgvert;   //!
    TBranch        *b_Nvcvert;   //!
    TBranch        *b_Dirvert;   //!
    TBranch        *b_Abspvert;   //!
    TBranch        *b_Abstpvert;   //!
    TBranch        *b_Ipvert;   //!
    TBranch        *b_Iverti;   //!
    TBranch        *b_Ivertf;   //!
    TBranch        *b_Fsiprob;   //!
    TBranch        *b_posc;   //!
    TBranch        *b_EventType;
    TBranch        *b_FlxWeight;
    TBranch        *b_OscWeight;
    TBranch        *b_pi0mass;
  
    // Vectors of splines, graphs, and histograms
    std::vector< std::vector< TH2F* > > dev_tmp;
    std::vector< std::string >          dev_tmp_names;
    std::vector< std::vector< std::vector< TSpline3*  >  >  > splines;
    std::vector< std::vector< std::vector< TGraph*  >  >  >   graphs;
  
    // Option for making nue or numu splines (0=nue, 1=numu, any other value indicates an error) [Note: this is currently not used any more - 22/01/2015]
    // Default value (set in constructor if no other value given) is 3
    int SampleType;
  
    XsecVary(std::string weight_file, std::string kinematics_file, std::vector<SystematicProperties> systProps, int sampletype=3); 
    virtual ~XsecVary();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Int_t    GetEntries();
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init();
    virtual void     AddSystematics(std::vector<SystematicProperties> &systProps);
    virtual void     MakeVariations();
    virtual void     WriteGraphs(std::string outputname);
    virtual void     SetBinning(const Double_t *ebins, Int_t nebins, 
         		       const Double_t *rebins, Int_t nrebins);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
  
};

#endif

#ifdef XsecVary_cxx
XsecVary::XsecVary(std::string weight_file, std::string kinematics_file, std::vector<SystematicProperties> systProps, int sampletype)
{ 
  std::cout << "XsecVary constructor with sampletype = " << sampletype << std::endl;
  // Set SampleType (should be 0 for nue or 1 for numu)
  SampleType = sampletype;

  systematicProperties = systProps;

  weightFile = weight_file;

  fChain = new TChain("caf");//h1/mtuple/minituple
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

XsecVary::~XsecVary()
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

}

SystematicProperties::SystematicProperties(std::string name, std::string sname, std::vector<int> modes, int nomLoc)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = false;  // default
}

SystematicProperties::SystematicProperties(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = add;
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
  //std::string modeToName[] = {"ccqe","cc1pi","cccoh","ccmisc","ncpiz","ncpipm","nccoh","ncoth","mec", "nc1gamma", "ccmpi", "ccdis"};//ETA adding ccmpi and ccdis for 2020OA, ccoth now ccmisc also was nc1gamma missing??
  std::string modeToName[] = {"qe", "mec", "dis", "res", "coh", "diff", "nueel", "unknown", "amnugamma", "unknown", "cohel", "ibd", "glasres", "imdannihilation"};

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
  fChain->SetMakeClass(1);

  //fChain->SetBranchAddress("iclass", &iclass, &b_iclass);
  //fChain->SetBranchAddress("Ibound", &Ibound, &b_Ibound);
  //fChain->SetBranchAddress("fqmomm", &fqmomm, &b_fqmomm);
  //fChain->SetBranchAddress("fqmome", &fqmome, &b_fqmome);
  //fChain->SetBranchAddress("fqmrmom", &fqmrmom, &b_fqmrmom);
  //fChain->SetBranchAddress("fqmreloss", &fqmreloss, &b_fqmreloss);
  //fChain->SetBranchAddress("fqmrdir", &fqmrdir, &b_fqmrdir);
  //fChain->SetBranchAddress("fq1rdir", &fq1rdir, &b_fq1rdir);
  //fChain->SetBranchAddress("ipp", &ipp, &b_ipp);
  //fChain->SetBranchAddress("erec", &Erec_check, &b_erec);
  //fChain->SetBranchAddress("pizmass", &pizmass, &b_pizmass);
  //fChain->SetBranchAddress("pizmom", &pizmom, &b_pizmom);
  //fChain->SetBranchAddress("pizcosb", &pizcosb, &b_pizcosb);
  //fChain->SetBranchAddress("wgtflx", &wgtflx, &b_wgtflx);
  //fChain->SetBranchAddress("wgtosc", &wgtosc, &b_wgtosc);
  //fChain->SetBranchAddress("dir", dir, &b_dir);
  //fChain->SetBranchAddress("amome", amome, &b_amome);
  //fChain->SetBranchAddress("npar", &npar, &b_npar);
  //fChain->SetBranchAddress("wallv", &wallv, &b_wallv);
  //fChain->SetBranchAddress("ipv", ipv, &b_ipv);
  //fChain->SetBranchAddress("posv", posv, &b_posv);
  //fChain->SetBranchAddress("dirv", dirv, &b_dirv);
  //fChain->SetBranchAddress("pmomv", pmomv, &b_pmomv);
  //fChain->SetBranchAddress("npar2", &npar2, &b_npar2);
  //fChain->SetBranchAddress("ipv2", ipv2, &b_ipv2);
  //fChain->SetBranchAddress("numnu", &numnu, &b_numnu);
  //fChain->SetBranchAddress("mode", &mode, &b_mode);
  //fChain->SetBranchAddress("nuPDG", ipnu, &b_ipnu);
  //fChain->SetBranchAddress("Ev", pnu, &b_pnu);
  //fChain->SetBranchAddress("dirnu", dirnu, &b_dirnu);
  //fChain->SetBranchAddress("Npvc", &Npvc, &b_Npvc);
  //fChain->SetBranchAddress("Ipvc", Ipvc, &b_Ipvc);
  //fChain->SetBranchAddress("Ichvc", Ichvc, &b_Ichvc);
  //fChain->SetBranchAddress("Iorgvc", Iorgvc, &b_Iorgvc);
  //fChain->SetBranchAddress("Iflvc", Iflvc, &b_Iflvc);
  //fChain->SetBranchAddress("Abspvc", Abspvc, &b_Abspvc);
  //fChain->SetBranchAddress("Pvc", Pvc, &b_Pvc);
  //fChain->SetBranchAddress("Crsx", &Crsx, &b_Crsx);
  //fChain->SetBranchAddress("Crsy", &Crsy, &b_Crsy);
  //fChain->SetBranchAddress("Crsz", &Crsz, &b_Crsz);
  //fChain->SetBranchAddress("Crsphi", &Crsphi, &b_Crsphi);
  //fChain->SetBranchAddress("Nvert", &Nvert, &b_Nvert);
  //fChain->SetBranchAddress("Posvert", Posvert, &b_Posvert);
  //fChain->SetBranchAddress("Iflgvert", Iflgvert, &b_Iflgvert);
  //fChain->SetBranchAddress("Nvcvert", &Nvcvert, &b_Nvcvert);
  //fChain->SetBranchAddress("Dirvert", Dirvert, &b_Dirvert);
  //fChain->SetBranchAddress("Abspvert", Abspvert, &b_Abspvert);
  //fChain->SetBranchAddress("Abstpvert", Abstpvert, &b_Abstpvert);
  //fChain->SetBranchAddress("Ipvert", Ipvert, &b_Ipvert);
  //fChain->SetBranchAddress("Iverti", Iverti, &b_Iverti);
  //fChain->SetBranchAddress("Ivertf", Ivertf, &b_Ivertf);
  //fChain->SetBranchAddress("Fsiprob", &Fsiprob, &b_Fsiprob);
  //fChain->SetBranchAddress("posc", posc, &b_posc);
  /* fChain->SetBranchAddress("EventType", &EventType, &b_EventType); */
  /* fChain->SetBranchAddress("FlxWeight", &FlxWeight, &b_FlxWeight); */
  /* fChain->SetBranchAddress("OscWeight", &OscWeight, &b_OscWeight); */
  /* fChain->SetBranchAddress("pi0mass", &pi0mass, &b_pi0mass); */

  // DUNE CAF Variables
  fChain->SetBranchAddress("cvnnumu", &cvnnumu);
  fChain->SetBranchAddress("cvnnue", &cvnnue);
  fChain->SetBranchAddress("Ev_reco_numu", &erec_numu);
  fChain->SetBranchAddress("Ev_reco_nue", &erec_nue);
  fChain->SetBranchAddress("nuPDG", ipnu, &b_ipnu);
  fChain->SetBranchAddress("Ev", pnu, &b_pnu);
  fChain->SetBranchAddress("mode", &mode, &b_mode);
  fChain->SetBranchAddress("wgtosc", &wgtosc);
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

#endif // #ifdef XsecVary_cxx
