#ifndef ttbar_analysis_h
#define ttbar_analysis_h

#include "AnalyzerCore.h"

class ttbar_analysis : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  std::vector<TString> muonTightIDs, muonLooseIDs;
  std::vector<TString> electronTightIDs, electronLooseIDs;
  std::vector<TString> jetIDs;

  std::vector<TString> muonTriggers, electronTriggers, emuTriggers, mueTriggers;

  double muonPtCut1, electronPtCut1, emuPtCut1, muePtCut1;
  double muonPtCut2, electronPtCut2, emuPtCut2, muePtCut2;
  double muonPtCut, electronPtCut, jetPtCut;

  std::vector<Muon> allMuons;
  std::vector<Electron> allElectrons;
  std::vector<Jet> allJets;

  TLorentzVector FindZCandidate(std::vector<Lepton> leptons);
  TLorentzVector FindZCandidate(std::vector<Muon> leptons);
  TLorentzVector FindZCandidate(std::vector<Electron> leptons);

  ttbar_analysis();
  ~ttbar_analysis();

};

#endif
