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

typedef Double_t FLOAT_T;

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

    //If true, the spline weights will be normalized by the event weight. Should be used for GENIEReweight systs
    bool reweight;

    //Index of the corresponding syst in the ifile
    int syst_idx;

    TBranch* weightBranch;
    TBranch* knotBranch;

    // constructor
    SystematicProperties(std::string name, std::string sname, std::vector<int> modes, int nomLoc);
    SystematicProperties(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add);
    SystematicProperties(std::string name, std::string sname, std::vector<int> modes, int nomLoc, bool add, bool reweight=false);

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
    Int_t                     fCurrent_w = 0;
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
    Double_t        wgtxsec;
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
    Int_t           ipnu;
    FLOAT_T         pnu;
    
    // List of branches: fChain
    TBranch        *b_iclass;   //!
	TBranch        *b_Ibound;
    TBranch        *b_erec;   //!
    TBranch        *b_mode;   //!
    TBranch        *b_ipnu;   //!
    TBranch        *b_pnu;   //!
    

  //DUNE VARIABLES FD
    FLOAT_T        erec_numu;
    FLOAT_T        erec_nue;
    FLOAT_T        cvnnumu;
    FLOAT_T        cvnnue;
    Int_t           isCC;
    FLOAT_T        vtx_x;
    FLOAT_T        vtx_y;
    FLOAT_T        vtx_z;
    FLOAT_T        berpacv;
    std::vector<std::vector<float>> *weightVec;

    // Declaration of leaf types: fTree_byEv
    FLOAT_T         Etrue;
    FLOAT_T         Erec_check;
  
  
    // Vectors of splines, graphs, and histograms
    std::vector< std::vector< TH3F* > > dev_tmp;
    std::vector< std::string >          dev_tmp_names;
    std::vector< std::vector< std::vector< TSpline3*  >  >  > splines;
    std::vector< std::vector< std::vector< TGraph*  >  >  >   graphs;
  
    // Option for making nue or numu splines (0=nue, 1=numu, any other value indicates an error) [Note: this is currently not used any more - 22/01/2015]
    // Default value (set in constructor if no other value given) is 3
    int SampleType;
  
    XsecVary(std::string weight_file, std::string kinematics_file, std::vector<SystematicProperties> systProps, int sampletype=3); 
    XsecVary(std::string sum_file, std::vector<SystematicProperties> systProps, int sampletype=3);
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

// #ifdef XsecVary_cxx
// #endif // #ifdef XsecVary_cxx
