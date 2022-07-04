#ifndef Muon_h
#define Muon_h

#include "Lepton.h"
#include "TString.h"

class Muon: public Lepton{

public:

  Muon();
  ~Muon();

  void SetTypeBit(unsigned int typebit);
  void SetIDBit(unsigned int idbit);

  enum Selector {
    CutBasedIdLoose        = 1UL<< 0,  
    CutBasedIdMedium       = 1UL<< 1,  
    CutBasedIdMediumPrompt = 1UL<< 2,  // medium with IP cuts
    CutBasedIdTight        = 1UL<< 3,  
    CutBasedIdGlobalHighPt = 1UL<< 4,  // high pt muon for Z',W' (better momentum resolution)
    CutBasedIdTrkHighPt    = 1UL<< 5,  // high pt muon for boosted Z (better efficiency)
    PFIsoVeryLoose         = 1UL<< 6,  // reliso<0.40
    PFIsoLoose             = 1UL<< 7,  // reliso<0.25
    PFIsoMedium            = 1UL<< 8,  // reliso<0.20
    PFIsoTight             = 1UL<< 9,  // reliso<0.15
    PFIsoVeryTight         = 1UL<<10,  // reliso<0.10
    TkIsoLoose             = 1UL<<11,  // reliso<0.10
    TkIsoTight             = 1UL<<12,  // reliso<0.05
    SoftCutBasedId         = 1UL<<13,  
    SoftMvaId              = 1UL<<14,  
    MvaLoose               = 1UL<<15,  
    MvaMedium              = 1UL<<16,  
    MvaTight               = 1UL<<17,
    MiniIsoLoose           = 1UL<<18,  // reliso<0.40
    MiniIsoMedium          = 1UL<<19,  // reliso<0.20
    MiniIsoTight           = 1UL<<20,  // reliso<0.10
    MiniIsoVeryTight       = 1UL<<21,  // reliso<0.05
    TriggerIdLoose         = 1UL<<22,  // robust selector for HLT
    InTimeMuon             = 1UL<<23,   
    PFIsoVeryVeryTight     = 1UL<<24,  // reliso<0.05
    MultiIsoLoose          = 1UL<<25,  // miniIso with ptRatio and ptRel 
    MultiIsoMedium         = 1UL<<26   // miniIso with ptRatio and ptRel 
  };

  enum Type {
    GlobalMuon     =  1<<1,
    TrackerMuon    =  1<<2,
    StandAloneMuon =  1<<3,
    CaloMuon =  1<<4,
    PFMuon =  1<<5,
    RPCMuon =  1<<6,
    GEMMuon =  1<<7,
    ME0Muon = 1<<8
  };

  inline bool PassSelector( unsigned int s ) const { return (j_IDBit & s)==s; }
  inline bool IsType( unsigned int t ) const { return (j_TypeBit & t); }

  inline bool isPOGTight() const {return PassSelector(CutBasedIdTight);}
  inline bool isPOGMedium() const {return PassSelector(CutBasedIdMedium);}
  inline bool isPOGLoose() const {return PassSelector(CutBasedIdLoose);}

  void SetPOGMediumHIP(bool ismedium_hip, bool ismedium_nohip);
  inline bool isPOGMedium_hip() const {return j_ismedium_hip;} //hip mitigated ID. Default for 2016preVFP data
  inline bool isPOGMedium_nohip() const {return j_ismedium_nohip;} //ordinary medium ID. Default for 2016 MC

  //==== TODO isOLDPOGHighPt returns values from bit, which is before the update
  //==== In UltraLegacy, this should be removed
  inline bool isOLDPOGHighPt() const {return PassSelector(CutBasedIdGlobalHighPt);}
  void SetisPOGHighPt(bool b);
  inline bool isPOGHighPt() const {return j_isPOGHighPt;}

  void SetIso(double ch04, double nh04, double ph04, double pu04, double trkiso);
  void CalcPFRelIso();
  inline double TrkIso() const {return j_trkiso;}
  double EA();

  void SetChi2(double chi2);
  inline double Chi2() const { return j_chi2; }

  void SetMiniAODPt(double d);
  void SetMiniAODTunePPt(double d);
  inline double MiniAODPt() const {return j_MiniAODPt;}
  inline double MiniAODTunePPt() const {return j_MiniAODTunePPt;}

  void SetMomentumScaleUpDown(double pt_up, double pt_down);
  inline double MomentumShift(int s) const {
    if(s==0) return Pt();
    else if(s>0) return j_MomentumScaleUp;
    else         return j_MomentumScaleDown;
  }

  void SetTuneP4(double pt, double pt_err, double eta, double phi, double q);
  inline Particle TuneP4() const {return j_TuneP4;}
  inline double TunePPtError() const {return j_TunePPtError;}

  //==== TODO lowptMVA not supported in CMSSW_10_2_10, only filling MVA now
  void SetMVA(double MVA);
  inline double MVA() const { return j_MVA; }
  inline double lowptMVA() const { return j_lowptMVA; }
  inline double softMVA() const { return j_softMVA; }

  //==== ID
  bool PassID(TString ID) const;
  bool Pass_POGTightWithTightIso() const;
  bool Pass_POGHighPtWithLooseTrkIso() const;
  bool Pass_TESTID() const;
  
  void SetTrackerLayers(int n);
  inline int TrackerLayers() const { return j_trackerLayers; }

  void SetFilterBits(ULong64_t filterbits){ j_filterbits=filterbits; }
  void SetPathBits(ULong64_t pathbits){ j_pathbits=pathbits; }
  bool PassFilter(TString filter) const;
  bool PassPath(TString path) const;


private:

  unsigned int j_TypeBit, j_IDBit;
  bool j_isPOGHighPt;
  bool j_ismedium_hip;
  bool j_ismedium_nohip;
  double j_chi2;
  double j_PFCH04, j_PFNH04, j_PFPH04, j_PU04, j_trkiso;
  double j_MiniAODPt, j_MiniAODTunePPt, j_MomentumScaleUp, j_MomentumScaleDown;
  Particle j_TuneP4;
  double j_TunePPtError;
  double j_MVA, j_lowptMVA, j_softMVA;
  int j_trackerLayers;

  ULong64_t j_filterbits;
  ULong64_t j_pathbits;

  ClassDef(Muon,1);
};

#endif
