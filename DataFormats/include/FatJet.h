#ifndef FatJet_h
#define FatJet_h

#include "Particle.h"
#include "JetTaggingParameters.h"

class FatJet: public Particle {
public:

  FatJet();
  ~FatJet();

  void SetArea(double area);
  void SetGenFlavours(double pf, double hf);
  inline int partonFlavour() const { return j_partonFlavour; }
  inline int hadronFlavour() const { return j_hadronFlavour; }
  void SetTaggerResults(std::vector<double> ds);
  void SetEnergyFractions(double cH, double nH, double nEM, double cEM, double muE);
  void SetMultiplicities(double cM, double nM);
  void SetLSF(double lsf, int lsf_PID);
  double LSF() const;
  double LSF_PID() const;

  void SetEnShift(double en_up, double en_down);
  double EnShift(int s){
    if(s==0) return 1.;
    else if(s>0) return j_En_up;
    else return j_En_down;
  }
  void SetResShift(double res_up, double res_down);
  double ResShift(int s){
    if(s==0) return 1.;
    else if(s>0) return j_Res_up;
    else return j_Res_down;
  }

  void SetTightJetID(double b);
  void SetTightLepVetoJetID(double b);
  inline bool Pass_tightJetID() const { return j_tightJetID; }
  inline bool Pass_tightLepVetoJetID() const { return j_tightLepVetoJetID; }

  bool PassID(TString ID) const;

  double GetTaggerResult(JetTagging::Tagger tg) const;

  void SetPuppiTaus(double t1, double t2, double t3, double t4);
  double PuppiTau1() const { return j_puppi_tau1; }
  double PuppiTau2() const { return j_puppi_tau2; }
  double PuppiTau3() const { return j_puppi_tau3; }
  double PuppiTau4() const { return j_puppi_tau4; }

  void SetSDMass(double m);
  double SDMass() const { return j_SDMass; }

private:

  double  j_area;
  int j_partonFlavour;
  int j_hadronFlavour;
  double  j_DeepCSV;
  double  j_DeepCSV_CvsL;
  double  j_DeepCSV_CvsB;
  double  j_particleNet_TvsQCD;
  double  j_particleNet_WvsQCD;
  double  j_particleNet_ZvsQCD;
  double  j_particleNet_HbbvsQCD;
  double  j_particleNet_HccvsQCD;
  double  j_particleNet_H4qvsQCD;
  double  j_particleNet_QCD;
  double  j_particleNetMD_Xbb;
  double  j_particleNetMD_Xcc;
  double  j_particleNetMD_Xqq;
  double  j_particleNetMD_QCD;
  double  j_chargedHadronEnergyFraction;
  double  j_neutralHadronEnergyFraction;
  double  j_neutralEmEnergyFraction;
  double  j_chargedEmEnergyFraction;
  double j_muonEnergyFraction;
  int j_chargedMultiplicity;
  int j_neutralMultiplicity;
  double j_lsf;
  int j_lsf_pid;
  double j_En_up;
  double j_En_down;
  double j_Res_up;
  double j_Res_down;
  bool j_tightJetID, j_tightLepVetoJetID;
  double j_puppi_tau1, j_puppi_tau2, j_puppi_tau3, j_puppi_tau4;
  double j_SDMass;

  ClassDef(FatJet,1);
};



#endif
