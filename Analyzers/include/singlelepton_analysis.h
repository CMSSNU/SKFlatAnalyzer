#ifndef singlelepton_analysis_h
#define singlelepton_analysis_h

#include "AnalyzerCore.h"

class singlelepton_analysis : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  std::vector<TString> muonTightIDs, muonLooseIDs;
  std::vector<TString> electronTightIDs, electronLooseIDs;
  std::vector<TString> jetIDs, fatjetIDs;

  std::vector<JetTagging::Parameters> bTaggings;
  JetTagging::Parameters bTaggingLooseWP, bTaggingMediumWP, bTaggingTightWP;

  std::vector<TString> muonTriggers, electronTriggers;

  double muonPtCut, electronPtCut, jetPtCut, fatjetPtCut;
  double leptonPtCut;
  void DrawObjectPlots(TString region, double weight, std::vector<Muon> muons, std::vector<Electron> electrons, std::vector<Jet> jets, std::vector<Jet> nonbjets, std::vector<Jet> bjets, std::vector<FatJet> fatjets, std::vector<FatJet> nonbbfatjets, std::vector<FatJet> bbfatjets, Particle missingEt);

  std::vector<Muon> allMuons;
  std::vector<Electron> allElectrons;
  std::vector<Jet> allJets;
  std::vector<FatJet> allFatJets;

  double GetReconstructedNeutrinoDet(Lepton lepton, Particle missingEt);
  Particle GetReconstructedNeutrino(Lepton lepton, Particle missingEt);

  TLorentzVector FindZCandidate(std::vector<Lepton> leptons);
  TLorentzVector FindZCandidate(std::vector<Muon> leptons);
  TLorentzVector FindZCandidate(std::vector<Electron> leptons);

  singlelepton_analysis();
  ~singlelepton_analysis();

};

#endif
