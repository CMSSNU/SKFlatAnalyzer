#ifndef SkimTree_HighPt1L_h
#define SkimTree_HighPt1L_h

#include "AnalyzerCore.h"

class SkimTree_HighPt1L : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_HighPt1L();
  ~SkimTree_HighPt1L();

  TTree *newtree;

  vector<TString> triggers;
  void WriteHist();

  double LeptonPtCut, AK4JetPtCut, AK8JetPtCut;

};



#endif

