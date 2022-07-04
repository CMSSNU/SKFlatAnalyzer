#ifndef MCCorrection_h
#define MCCorrection_h

#include <map>
#include <vector>
#include <fstream>
#include <sstream>

#include "TKey.h"
#include "TFile.h"
#include "TString.h"
#include "TRegexp.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TObjString.h"
#include "TSystem.h"

#include "Lepton.h"
#include "Muon.h"
#include "Electron.h"
#include "Photon.h"
#include "FatJet.h"
#include "Jet.h"
#include "Gen.h"

#include "TDirectoryHelper.h"
#include "GenFinderForDY.h"
#include "BTagCalibrationStandalone.h"
#include "TRandom3.h"

class MCCorrection{

public:

  MCCorrection();
  ~MCCorrection();

  TDirectory *histDir;
  static vector<TString> Split(TString s,TString del);
  void ReadHistograms();

  TString MCSample;
  void SetMCSample(TString s);

  int DataYear;
  TString DataEra;
  void SetEra(TString era){
    DataEra=era;
    DataYear=TString(era(0,4)).Atoi();
  }
  TString GetEra() const { return DataEra; }
  int GetYear() const { return DataYear; }

  bool IsDATA;
  void SetIsDATA(bool b);

  int run,lumi,event;
  void SetEventInfo(int r, int l, int e);

  bool IsFastSim;
  void SetIsFastSim(bool b);

  bool IgnoreNoHist;

  double MuonReco_SF(TString key, double eta, double p, int sys=0);
  double MuonID_SF(TString ID, double eta, double pt, int sys=0);
  double MuonISO_SF(TString ID, double eta, double pt, int sys=0);
  double MuonTrigger_Eff(TString ID, TString trig, int DataOrMC, double eta, double pt, int sys=0);
  double MuonTrigger_SF(TString ID, TString trig, const std::vector<Muon>& muons, int sys=0);
  double MuonTrigger_SF(TString ID, TString trig, const std::vector<Muon *>& muons, int sys=0);

  std::map< TString, TH2F* > map_hist_Muon;

  double ElectronReco_SF(double sceta, double pt, int sys=0);
  double ElectronID_SF(TString ID, double sceta, double pt, int sys=0);
  double ElectronTrigger_Eff(TString ID, TString trig, int DataOrMC, double eta, double pt, int sys=0);
  double ElectronTrigger_SF(TString ID, TString trig, const std::vector<Electron>& electrons, int sys=0);
  double ElectronTrigger_SF(TString ID, TString trig, const std::vector<Electron *>& electrons, int sys=0);
  std::map< TString, TH2F* > map_hist_Electron;
  std::map< TString, TGraphAsymmErrors* > map_graph_Electron;

  std::map< TString, TH2F* > map_hist_prefire;
  double GetPrefireWeight(const std::vector<Photon>& photons, const std::vector<Jet>& jets, int sys);

  std::map< TString, TH1D* > map_hist_pileup;
  double GetPileUpWeightBySampleName(int N_pileup, int syst);
  double GetPileUpWeight(int N_pileup, int syst);

  double GetTopPtReweight(const std::vector<Gen>& gens);

  TH2D *hist_DYPtReweight_2D;
  GenFinderForDY *genFinderDY;
  double GetOfficialDYReweight(const vector<Gen>& gens, int sys=0);

  //==== b tagging

  double GetJetTaggingCutValue(JetTagging::Tagger tagger, JetTagging::WP wp);

  std::vector<JetTagging::Parameters> jetTaggingPars;
  void SetJetTaggingParameters(std::vector<JetTagging::Parameters> v);

  std::map< std::string, BTagCalibrationReader* > map_BTagCalibrationReader;
  void SetupJetTagging();

  TH2D *hist_JetTagEff_B;
  TH2D *hist_JetTagEff_C;
  TH2D *hist_JetTagEff_Light;

  std::map< TString, TH2F* > map_hist_mcjet;
  void SetupMCJetTagEff();
  double GetMCJetTagEff(JetTagging::Tagger tagger, JetTagging::WP wp, int JetFlavor, double JetPt, double JetEta, int sys=0);
  double GetJetTaggingSF(JetTagging::Parameters jtp, int JetFlavor, double JetPt, double JetEta, double Jetdiscr, string Syst="central");

  //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods

  //==== 1a) Event reweighting using scale factors and MC b-tagging efficiencies
  double GetBTaggingReweight_1a(const vector<Jet>& jets, JetTagging::Parameters jtp, string Syst="central");
  //==== 1d) Event reweighting using discriminator-dependent scale factors
  double GetBTaggingReweight_1d(const vector<Jet>& jets, JetTagging::Parameters jtp, string Syst="central");

  //==== 2a) Jet-by-jet updating of the b-tagging status
  bool IsBTagged_2a(JetTagging::Parameters jtp, const Jet& jet, string Syst="central");
  bool IsBTagged_2a(JetTagging::Parameters jtp, const FatJet& jet, string Syst="central");

  double TopPtReweight(const std::vector<Gen>& gens );
};

#endif
