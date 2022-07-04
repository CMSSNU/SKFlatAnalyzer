#include "signal_study.h"

signal_study::signal_study(){

}

void signal_study::initializeAnalyzer(){

    muonTriggersHighPt.clear();
    muonTriggersIso.clear();

    if (DataYear == 2017){
        muonTriggersHighPt.push_back("HLT_Mu50_v");
        muonTriggersHighPt.push_back("HLT_OldMu100_v");
        muonTriggersHighPt.push_back("HLT_TkMu100_v");
        muonTriggersIso.push_back("HLT_IsoMu27_v");
    }

    bTaggingLooseWP = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb);
    bTaggingMediumWP = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
    bTaggingTightWP = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Tight, JetTagging::incl, JetTagging::comb);

    bTaggings.push_back(bTaggingLooseWP);
    bTaggings.push_back(bTaggingMediumWP);
    bTaggings.push_back(bTaggingTightWP);

    mcCorr->SetJetTaggingParameters(bTaggings);

}


signal_study::~signal_study(){

}

void signal_study::executeEvent(){

    allMuons = GetAllMuons();
    allElectrons = GetAllElectrons();
    allJets = GetAllJets();
    allFatJets = puppiCorr->Correct(GetAllFatJets());

    AnalyzerParameter param; param.Clear();
    param.syst_ = AnalyzerParameter::Central;

    executeEventFromParameter(param);

}

void signal_study::executeEventFromParameter(AnalyzerParameter param){

    Event event = GetEvent();

    double weight = 1.;
    if (!IsDATA){
        weight = weight * event.MCweight();
        if (MCSample.Contains("TypeIHeavyN")){
            if (gen_check(weight) < 0) return;
        }
    }

    if(!PassMETFilter()) return;

//    if (!(event.PassTrigger(muonTriggers) || event.PassTrigger(electronTriggers))) return;

    bool passHighPtHLT = false, passIsoHLT = false;
    float muonPtCut = -1.;
    std::vector<Muon> muonsHighPt, muonsHighPtLoose, muonsIso, muonsIsoLoose;
    if (event.PassTrigger(muonTriggersHighPt)){
        passHighPtHLT = true;
        muonPtCut = 55.;
        muonsHighPtLoose = SelectMuons(UseTunePMuon(allMuons), "POGHighPt", muonPtCut, 2.4);
        muonsHighPt = SelectMuons(muonsHighPtLoose, "POGHighPtWithLooseTrkIso", muonPtCut, 2.4);
    }
    if (event.PassTrigger(muonTriggersIso)){
        passIsoHLT = true;
        muonPtCut = 30.;
        muonsIsoLoose = SelectMuons(allMuons, "POGMedium", muonPtCut, 2.4);
        for (unsigned int i=0; i<muonsIsoLoose.size(); i++){
            Muon muon = muonsIsoLoose.at(i);
            if (muon.MiniRelIso() < 0.2) muonsIso.push_back(muon);
        }
    }

    std::vector<Electron> electrons = SelectElectrons(allElectrons, "passMediumID", 20., 2.4);

    Particle missingEt = event.GetMETVector();

    if (!(electrons.size() == 0)) return;

    if (passHighPtHLT){
        std::vector<Muon> muons = muonsHighPt;
        if (muons.size() == 1){
            fill_histograms("HighPt", weight, muons, missingEt);
        }
    }
    if (passIsoHLT){
        std::vector<Muon> muons = muonsIso;
        if (muons.size() == 1){
            fill_histograms("Iso", weight, muons, missingEt);
        }
    }

    return;
}

void signal_study::fill_histograms(TString selection, double weight, std::vector<Muon> muons, Particle missingEt){

     FillHist(selection + "_pt_muon", muons.at(0).Pt(), weight, 2000, 0., 2000.);
     FillHist(selection + "_eta_muon", muons.at(0).Eta(), weight, 50, -5., 5.);
     FillHist(selection + "_reliso_muon", muons.at(0).RelIso(), weight, 20, 0., 1.);
     FillHist(selection + "_miniiso_muon", muons.at(0).MiniRelIso(), weight, 20, 0., 1.);
     FillHist(selection + "_met", missingEt.Pt(), weight, 2000, 0., 2000.);
     FillHist(selection + "_masst_muonmet", (muons.at(0) + missingEt).Mt(), weight, 2000, 0., 2000.);

}

int signal_study::gen_check(double weight){

    std::vector<Gen> allGens = GetGens();

    std::vector<Gen> quarksGEN, bquarksGEN;
    std::vector<Gen> muonGEN, heavyneutrinoGEN, secondarybosonGEN;

    for (unsigned int i=0 ; i < allGens.size() ; i++ ){
        Gen particle = allGens.at(i);
        if (particle.isHardProcess() && particle.Status() != 21){
            int motherIndex = particle.MotherIndex();
            int motherPdgId = abs(allGens.at(motherIndex).PID());
            int motherStatus = allGens.at(motherIndex).Status();
            int particlePdgId = abs(particle.PID());
            if (particlePdgId == 9900012){
                heavyneutrinoGEN.push_back(particle);
            }
        }
    }

    for (unsigned int i=0 ; i < allGens.size() ; i++ ){
        Gen particle = allGens.at(i);
        if (particle.isHardProcess() && particle.Status() != 21){
            int motherIndex = particle.MotherIndex();
            int motherPdgId = abs(allGens.at(motherIndex).PID());
            int motherStatus = allGens.at(motherIndex).Status();
            int particlePdgId = abs(particle.PID());
            if (particlePdgId == 13){
                if (heavyneutrinoGEN.at(0).MotherIndex() == particle.MotherIndex()){
                    muonGEN.push_back(particle);
                }
                if (motherPdgId == 9900012){
                    muonGEN.push_back(particle);
                }
            }
        }
    }


    for (unsigned int i=0 ; i < allGens.size() ; i++ ){
        Gen particle = allGens.at(i);
        if (particle.isHardProcess() && particle.Status() != 21){
            int motherIndex = particle.MotherIndex();
            int motherPdgId = abs(allGens.at(motherIndex).PID());
            int motherStatus = allGens.at(motherIndex).Status();
            int particlePdgId = abs(particle.PID());
            if (particlePdgId == 23 || particlePdgId == 24 || particlePdgId == 25){
                if (motherPdgId == 9900012){
                    secondarybosonGEN.push_back(particle);
                }
            }
        }
    }

    for (unsigned int i=0 ; i < allGens.size() ; i++ ){
        Gen particle = allGens.at(i);
        if (particle.isHardProcess() && particle.Status() != 21){
            int motherIndex = particle.MotherIndex();
            int motherPdgId = abs(allGens.at(motherIndex).PID());
            int motherStatus = allGens.at(motherIndex).Status();
            int particlePdgId = abs(particle.PID());
            if (particlePdgId <= 4){ 
                if (motherPdgId == abs(secondarybosonGEN.at(0).PID())){
                    quarksGEN.push_back(particle);
                }
            }
            if (particlePdgId == 5){
                if (motherPdgId == abs(secondarybosonGEN.at(0).PID())){
                    bquarksGEN.push_back(particle);
                }
            }
        }
    }

    bool qqchannel = false, bbZchannel = false, bbHchannel = false;
    if (muonGEN.size() == 1){
        if (heavyneutrinoGEN.size() == 1){
            if (secondarybosonGEN.size() == 1){
                if (quarksGEN.size() == 2){
                    qqchannel = true;
                }
                if (bquarksGEN.size() == 2){
                    if (secondarybosonGEN.at(0).PID() == 25){
                        bbHchannel = true;
                    }
                    if (secondarybosonGEN.at(0).PID() == 23){
                        bbZchannel = true;
                    }
                }
            }
        }
    }

    if (qqchannel + bbHchannel + bbZchannel >= 2){
        FillHist("GEN_error_both_turned_on_weighted", 0., weight, 1, 0., 1.);
        FillHist("GEN_error_both_turned_on_unweighted", 0., 1, 1, 0., 1.);
        for (unsigned int i=0 ; i < allGens.size() ; i++ ){
            Gen particle = allGens.at(i);
        }
        return (-1);
    }
    if (qqchannel + bbHchannel + bbZchannel == 0){
        FillHist("GEN_no_signal_weighted", 0., weight, 1, 0., 1.);
        FillHist("GEN_no_signal_unweighted", 0., 1, 1, 0., 1.);
        for (unsigned int i=0 ; i < allGens.size() ; i++ ){
            Gen particle = allGens.at(i);
        }
        return (-1);
    }

    double channel_index = -1.;
    if (qqchannel) channel_index = 0.;
    if (bbZchannel) channel_index = 1.;
    if (bbHchannel) channel_index = 2.;

    FillHist("GEN_singal_weighted", channel_index, weight, 3, 0., 3.);
    FillHist("GEN_singal_unweighted", channel_index, 1, 3, 0., 3.);

    if (qqchannel){
        FillHist("GEN_mindeltar_leptonq", std::min(muonGEN.at(0).DeltaR(quarksGEN.at(0)), muonGEN.at(0).DeltaR(quarksGEN.at(1))), weight, 50, 0., 5.);
    }

    return channel_index;

}
