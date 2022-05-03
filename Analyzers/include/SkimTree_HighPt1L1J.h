#ifndef SkimTree_HighPt1L1J_h
#define SkimTree_HighPt1L1J_h

#include "AnalyzerCore.h"

class SkimTree_HighPt1L1J : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_HighPt1L1J();
  ~SkimTree_HighPt1L1J();

  TTree *newtree;

  vector<TString> triggers;
  void WriteHist();

  double LeptonPtCut, AK8JetPtCut;

};



#endif

