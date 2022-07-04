#ifndef Jet_h
#define Jet_h

#include "Particle.h"
#include "JetTaggingParameters.h"

class Jet: public Particle {
public:

  Jet();
  ~Jet();

  void SetArea(double area);
  void SetGenFlavours(int pf, int hf);
  inline int partonFlavour() const { return j_partonFlavour; }
  inline int hadronFlavour() const { return j_hadronFlavour; }
  void SetGenHFHadronMatcher(int flavour, int origin);
  inline int GenHFHadronMatcherFlavour() const { return j_GenHFHadronMatcher_flavour; }
  inline int GenHFHadronMatcherOrigin() const { return j_GenHFHadronMatcher_origin; }
  inline bool IsGenMatched() const { return j_GenHFHadronMatcher_flavour!=-999; }
  void SetTaggerResults(std::vector<double> ds);
  void SetEnergyFractions(double cH, double nH, double nEM, double cEM, double muE);
  void SetMultiplicities(double cM, double nM);
  void SetPileupJetId(double v);

  void SetEnShift(double en_up, double en_down);
  inline double EnShift(int s) const {
    if(s==0) return 1.;
    else if(s>0) return j_En_up;
    else return j_En_down;
  }
  void SetResShift(double res_up, double res_down);
  inline double ResShift(int s) const {
    if(s==0) return 1.;
    else if(s>0) return j_Res_up;
    else return j_Res_down;
  }
  void SetBJetNNCorrection(double bJetNN_corr, double bJetNN_res);
  inline double BJetNNCorrection() const { return j_bJetNN_corr; }
  inline double BJetNNResolution() const { return j_bJetNN_res; }
  void SetCJetNNCorrection(double cJetNN_corr, double cJetNN_res);
  inline double CJetNNCorrection() const { return j_cJetNN_corr; }
  inline double CJetNNResolution() const { return j_cJetNN_res; }

  void SetTightJetID(double b);
  void SetTightLepVetoJetID(double b);
  inline bool Pass_tightJetID() const { return j_tightJetID; }
  inline bool Pass_tightLepVetoJetID() const { return j_tightLepVetoJetID; }

  bool PassID(TString ID) const;

  double GetTaggerResult(JetTagging::Tagger tg) const;

private:

  double  j_area;
  int j_partonFlavour;
  int j_hadronFlavour;
  int j_GenHFHadronMatcher_flavour;
  int j_GenHFHadronMatcher_origin;
  double  j_DeepCSV;
  double  j_DeepCSV_CvsL;
  double  j_DeepCSV_CvsB;
  double  j_DeepJet;
  double  j_DeepJet_CvsL;
  double  j_DeepJet_CvsB;
  double  j_chargedHadronEnergyFraction;
  double  j_neutralHadronEnergyFraction;
  double  j_neutralEmEnergyFraction;
  double  j_chargedEmEnergyFraction;
  double j_muonEnergyFraction;
  int j_chargedMultiplicity;
  int j_neutralMultiplicity;
  double j_PileupJetId;
  double j_En_up;
  double j_En_down;
  double j_Res_up;
  double j_Res_down;
  double j_bJetNN_corr;
  double j_bJetNN_res;
  double j_cJetNN_corr;
  double j_cJetNN_res;
  bool j_tightJetID, j_tightLepVetoJetID;

  ClassDef(Jet,1)
};

#endif
