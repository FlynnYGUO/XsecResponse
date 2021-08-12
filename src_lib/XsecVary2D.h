// Patrick Dunne, p.dunne12@imperial.ac.uk
// August 11, 2021
//
// Adapted from T2K  XsecResponse/src_lib/XsecVary20192D.cc

#ifndef XsecVary2D_h
#define XsecVary2D_h
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

class SystematicProperties2D
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
    SystematicProperties2D(std::string name, std::string sname, std::vector<int> modes, int nomLoc);
    SystematicProperties2D(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add);

    // getters
    const TArrayF * GetWeightArray();
    const TArrayF * GetKnotArray();
    std::vector<int> GetInteractionModes();

    //setters
    void SetArrays(TArrayF* warr, TArrayF* karr);
    void SetBranches(TBranch* wbranch, TBranch* kbranch);
 
};

class XsecVary2D {
  public :

    // Vector of systematics for which we want splines made, containing info regarding how
    // we want them made
    std::vector<SystematicProperties2D> systematicProperties;
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
    Int_t           iclass;
	Int_t           Ibound;
    //Double_t        erec;
    Double_t        pizmass;
    Double_t        pizmom;
    Double_t        pizcosb;
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
    Float_t         pnu[50];
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
    Float_t         t2kposc[2][6];
    Int_t           EventType;
    Double_t        FlxWeight;
    Double_t        OscWeight;
    Float_t         pi0mass[2];
    Double_t        erec;
    Double_t        erec_numu;
    Double_t        erec_nue;
    Double_t        cvnnumu;
    Double_t        cvnnue;
   
    Float_t         Etrue;
    //Float_t         Erec;
    Double_t         Erec;
  
    // List of branches: fChain
    TBranch        *b_iclass;   //!
	TBranch        *b_Ibound;
    TBranch        *b_erec;   //!
    TBranch        *b_pizmass;   //!
    TBranch        *b_pizmom;   //!
    TBranch        *b_pizcosb;   //!
    TBranch        *b_wgtflx;   //!
    TBranch        *b_wgtosc;   //!
    TBranch        *b_dir;   //!
    TBranch        *b_fq1rdir;   //!
    TBranch        *b_amome;   //!
    TBranch        *b_amomm;   //!
    TBranch        *b_fqmome;   //!
    TBranch        *b_fqmomm;   //!
	TBranch        *b_fqmrdir; //!
	TBranch        *b_fqmrmom; //!
	TBranch        *b_fqmreloss; //!
	TBranch        *b_ipp; //!
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
    TBranch        *b_t2kposc;   //!
    TBranch        *b_EventType;
    TBranch        *b_FlxWeight;
    TBranch        *b_OscWeight;
    TBranch        *b_pi0mass;
     
    // Vectors of splines, graphs, and histograms
    std::vector< std::vector< TH3F* > > dev_tmp;
    std::vector< std::string >          dev_tmp_names;
    std::vector< std::vector< std::vector< std::vector< TSpline3* > > > > splines;
    std::vector< std::vector< std::vector< std::vector< TGraph* > > > >   graphs;
  
    // Option for making nue or numu splines (0=nue, 1=numu, any other value indicates an error) [Note: this is currently not used any more - 22/01/2015]
    // Default value (set in constructor if no other value given) is 3
    int SampleType;
   
    XsecVary2D(std::string weight_file, std::string kinematics_file, std::vector<SystematicProperties2D> systProps, int sampletype=3, std::string treename="init",std::string second_par="theta",bool horn=0,std::string first_par="erec",bool is_fitqun=false); 
    virtual ~XsecVary2D();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Int_t    GetEntries();
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init();
    virtual void     AddSystematics(std::vector<SystematicProperties2D> &systProps);
    virtual void     MakeVariations();
    virtual void     WriteGraphs(char *outputname);
    virtual void     SetBinning(const Double_t *ebins, Int_t nebins, 
       		       const Double_t *rebins, Int_t nrebins, const Double_t *qsqdbins, Int_t nqsqdbins);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    virtual double   CalcQ2(double Enu, double Elep, double Plep, double costheta, double Mlep);
    virtual double   CalcXSecTerm(double Enu, double Elep, double Plep, double costheta, double Mlep, double Mn);
    
};

#endif

#ifdef XsecVary2D_cxx
XsecVary2D::XsecVary2D(std::string weight_file, std::string kinematics_file, std::vector<SystematicProperties2D> systProps , int sampletype, std::string treename, std::string second_par,bool horn,std::string first_par,bool is_fitqun)
{ 
  std::cout << "XsecVary2D constructor with sampletype = " << sampletype << std::endl;
  // Set SampleType (should be 0 for nue or 1 for numu)
  SampleType = sampletype;

  systematicProperties = systProps;

  weightFile = weight_file;

  isfitqun=is_fitqun;

  fChain = new TChain("mtuple");//h1/mtuple/minituple
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
  if(fChain->GetEntries()==0)
  {
    std::cout<<"Looking for tree named "<<treename<<" instead."<<std::endl;
    fChain = new TChain(treename.c_str());
    if(fChain) fChain->Add(kinematics_file.c_str());
  }

  Init();

  AddSystematics(systematicProperties);

  theta_2d=false;
  Q2_2d=false;
  N_2d=false;
  erec_1d=false;
  p_1d=false;
  if(std::string(second_par) == "theta") theta_2d = true;
  if(std::string(second_par) == "Q2") Q2_2d = true;
  if(std::string(second_par) == "xsec") N_2d = true;
  if(std::string(first_par) == "erec") erec_1d = true;
  if(std::string(first_par) == "p") p_1d = true;
}


XsecVary2D::~XsecVary2D()
{
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

  std::cout << "Deleting histograms " << std::endl;
  for(unsigned i = 0; i < dev_tmp.size(); i++)
  {
    for(unsigned j = 0; j < dev_tmp[i].size(); j++)
    {
      if(dev_tmp[i][j]) dev_tmp[i][j]->Delete();
    }
  }

}

SystematicProperties2D::SystematicProperties2D(std::string name, std::string sname, std::vector<int> modes, int nomLoc)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = false;  // default
}

SystematicProperties2D::SystematicProperties2D(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add)
{
  systName          = name;
  shortName         = sname;
  intModes          = modes;
  nominalPosition   = nomLoc;
  addKnots          = add;
}

void SystematicProperties2D::SetBranches(TBranch* wbranch, TBranch* kbranch)
{
  weightBranch = wbranch;
  knotBranch   = kbranch;
}

void SystematicProperties2D::SetArrays(TArrayF* warr, TArrayF* karr)
{
  weightArray = warr;
  knotArray   = karr;
}

const TArrayF * SystematicProperties2D::GetWeightArray(){return weightArray;}
const TArrayF * SystematicProperties2D::GetKnotArray(){return knotArray;}

Int_t XsecVary2D::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Int_t XsecVary2D::GetEntries()
{
  // Return number of entries in tree
  if (!fChain) return 0;
  return fChain->GetEntries();
}

void XsecVary2D::AddSystematics(std::vector<SystematicProperties2D> &systProps)
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
  std::cout << "numSplines = " << numSplines;
  fChain->GetEntry(0);
  unsigned int itty = 0;
  std::string modeToName[] = {"ccqe","cc1pi","cccoh","ccmisc","ncpiz","ncpipm","nccoh","ncoth","mec", "nc1gamma", "ccmpi", "ccdis"};//ETA : adding ccmpi and ccdis for 2020OA, ccoth now ccmisc
  for(int k = 0; k < numSysts; k++)
  {
    for(unsigned l = 0; l < systProps[k].intModes.size(); l++)
    {
      dev_tmp[itty].resize(systProps[k].GetKnotArray()->GetSize());
      itty++;
      char tmpname[1000];
      sprintf(tmpname,"dev_%s_%s",systProps[k].shortName.c_str(),modeToName[systProps[k].intModes[l]].c_str());
      dev_tmp_names.push_back(std::string(tmpname));
    }
  }

}

Long64_t XsecVary2D::LoadTree(Long64_t entry)
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


void XsecVary2D::Init()
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

  fChain->SetBranchAddress("iclass", &iclass, &b_iclass);
  fChain->SetBranchAddress("Ibound", &Ibound, &b_Ibound);
  fChain->SetBranchAddress("erec", &Erec, &b_erec);
  fChain->SetBranchAddress("pizmass", &pizmass, &b_pizmass);
  fChain->SetBranchAddress("pizmom", &pizmom, &b_pizmom);
  fChain->SetBranchAddress("pizcosb", &pizcosb, &b_pizcosb);
  fChain->SetBranchAddress("wgtflx", &wgtflx, &b_wgtflx);
  fChain->SetBranchAddress("wgtosc", &wgtosc, &b_wgtosc);
  fChain->SetBranchAddress("dir", dir, &b_dir);
  if(isfitqun)
  {
    fChain->SetBranchAddress("fq1rdir", fq1rdir, &b_fq1rdir);
    fChain->SetBranchAddress("fqmome", &fqmome, &b_fqmome);
    fChain->SetBranchAddress("fqmomm", &fqmomm, &b_fqmomm);
	//ETA - adding in fitqun multi-ring variables for numuCC1pi samples(s)
	fChain->SetBranchAddress("ipp", &ipp, &b_ipp);
	fChain->SetBranchAddress("fqmrdir", &fqmrdir, &b_fqmrdir);
	fChain->SetBranchAddress("fqmrmom", &fqmrmom, &b_fqmrmom);
	fChain->SetBranchAddress("fqmreloss", &fqmreloss, &b_fqmreloss);
	fChain->SetBranchAddress("ipp", &ipp, &b_ipp);
  }
  fChain->SetBranchAddress("amome", amome, &b_amome);
  fChain->SetBranchAddress("amomm", amomm, &b_amomm);
  fChain->SetBranchAddress("npar", &npar, &b_npar);
  fChain->SetBranchAddress("wallv", &wallv, &b_wallv);
  fChain->SetBranchAddress("ipv", ipv, &b_ipv);
  fChain->SetBranchAddress("posv", posv, &b_posv);
  fChain->SetBranchAddress("dirv", dirv, &b_dirv);
  fChain->SetBranchAddress("pmomv", pmomv, &b_pmomv);
  fChain->SetBranchAddress("npar2", &npar2, &b_npar2);
  fChain->SetBranchAddress("ipv2", ipv2, &b_ipv2);
  fChain->SetBranchAddress("numnu", &numnu, &b_numnu);
  fChain->SetBranchAddress("mode", &mode, &b_mode);
  fChain->SetBranchAddress("ipnu", ipnu, &b_ipnu);
  fChain->SetBranchAddress("pnu", pnu, &b_pnu);
  fChain->SetBranchAddress("dirnu", dirnu, &b_dirnu);
  fChain->SetBranchAddress("Npvc", &Npvc, &b_Npvc);
  fChain->SetBranchAddress("Ipvc", Ipvc, &b_Ipvc);
  fChain->SetBranchAddress("Ichvc", Ichvc, &b_Ichvc);
  fChain->SetBranchAddress("Iorgvc", Iorgvc, &b_Iorgvc);
  fChain->SetBranchAddress("Iflvc", Iflvc, &b_Iflvc);
  fChain->SetBranchAddress("Abspvc", Abspvc, &b_Abspvc);
  fChain->SetBranchAddress("Pvc", Pvc, &b_Pvc);
  fChain->SetBranchAddress("Crsx", &Crsx, &b_Crsx);
  fChain->SetBranchAddress("Crsy", &Crsy, &b_Crsy);
  fChain->SetBranchAddress("Crsz", &Crsz, &b_Crsz);
  fChain->SetBranchAddress("Crsphi", &Crsphi, &b_Crsphi);
  fChain->SetBranchAddress("Nvert", &Nvert, &b_Nvert);
  fChain->SetBranchAddress("Posvert", Posvert, &b_Posvert);
  fChain->SetBranchAddress("Iflgvert", Iflgvert, &b_Iflgvert);
  fChain->SetBranchAddress("Nvcvert", &Nvcvert, &b_Nvcvert);
  fChain->SetBranchAddress("Dirvert", Dirvert, &b_Dirvert);
  fChain->SetBranchAddress("Abspvert", Abspvert, &b_Abspvert);
  fChain->SetBranchAddress("Abstpvert", Abstpvert, &b_Abstpvert);
  fChain->SetBranchAddress("Ipvert", Ipvert, &b_Ipvert);
  fChain->SetBranchAddress("Iverti", Iverti, &b_Iverti);
  fChain->SetBranchAddress("Ivertf", Ivertf, &b_Ivertf);
  fChain->SetBranchAddress("Fsiprob", &Fsiprob, &b_Fsiprob);
  fChain->SetBranchAddress("t2kposc", t2kposc, &b_t2kposc);
  /* fChain->SetBranchAddress("EventType", &EventType, &b_EventType); */
  /* fChain->SetBranchAddress("FlxWeight", &FlxWeight, &b_FlxWeight); */
  /* fChain->SetBranchAddress("OscWeight", &OscWeight, &b_OscWeight); */
  /* fChain->SetBranchAddress("pi0mass", &pi0mass, &b_pi0mass); */

  fChain->SetBranchAddress("cvnnumu", &cvnnumu);
  fChain->SetBranchAddress("cvnnue", &cvnnue);
  fChain->SetBranchAddress("Ev_reco_numu", &erec_numu);
  fChain->SetBranchAddress("Ev_reco_nue", &erec_nue);
 
  Notify();
}

Bool_t XsecVary2D::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

   return kTRUE;
}

void XsecVary2D::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

#endif // #ifdef XsecVary_cxx
