#ifndef signal_study_h
#define signal_study_h

#include "AnalyzerCore.h"

class signal_study : public AnalyzerCore {

public:

    void initializeAnalyzer();

    void executeEventFromParameter(AnalyzerParameter param);
    void executeEvent();

    std::vector<JetTagging::Parameters> bTaggings;
    JetTagging::Parameters bTaggingLooseWP, bTaggingMediumWP, bTaggingTightWP;

    std::vector<Muon> allMuons;
    std::vector<Electron> allElectrons;
    std::vector<Jet> allJets;
    std::vector<FatJet> allFatJets;
    std::vector<TString> muonTriggersHighPt, muonTriggersIso;

    int gen_check(double weight);

    void fill_histograms(TString selection, double weight, std::vector<Muon> muons, Particle missingEt);

    signal_study();
    ~signal_study();

};

#endif
