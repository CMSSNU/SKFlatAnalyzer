#ifndef SkimTree_HighPt1LJets_h
#define SkimTree_HighPt1LJets_h

#include "AnalyzerCore.h"

class SkimTree_HighPt1LJets : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_HighPt1LJets();
  ~SkimTree_HighPt1LJets();

  TTree *newtree;

  vector<TString> triggers;
  void WriteHist();

  double LeptonPtCut, AK4JetPtCut, AK8JetPtCut;

};



#endif

