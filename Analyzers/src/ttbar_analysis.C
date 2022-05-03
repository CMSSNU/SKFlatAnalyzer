#include "ttbar_analysis.h"

ttbar_analysis::ttbar_analysis(){

}

void ttbar_analysis::initializeAnalyzer(){

    muonTightIDs.clear();
    muonLooseIDs.clear();
    electronTightIDs.clear();
    electronLooseIDs.clear();
    jetIDs.clear();

    muonTightIDs.push_back("POGTight");
    muonLooseIDs.push_back("POGMedium");

    electronTightIDs.push_back("passMVAID_iso_WP80");
    electronLooseIDs.push_back("passMVAID_iso_WP90");

    jetIDs.push_back("tight");

    muonTriggers.clear();
    electronTriggers.clear();
    emuTriggers.clear();
    mueTriggers.clear();

    muonPtCut = 10.; electronPtCut = 10.;

    if (DataYear == 2016){
//        muonTriggersBtoG.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"); // 27267.591112919 
//        muonTriggersBtoG.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"); // 27267.591112919
//        muonTriggersH.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"); // 8650.628380028
//        muonTriggersH.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"); // 8650.628380028
        muonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"); // 35918.219492947
        muonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"); // 35918.219492947 
        muonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"); // 8650.628380028
        muonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"); // 8650.628380028
        emuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"); 
        emuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
        mueTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
        mueTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
        muonPtCut1 = 20.; muonPtCut2 = 10.;
        electronPtCut1 = 25.; electronPtCut2 = 15.;
        emuPtCut1 = 25.; emuPtCut2 = 10.;
        muePtCut1 = 25.; muePtCut2 = 15.;
    }
    else if (DataYear == 2017){
        muonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
        electronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
        emuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
        mueTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
        muonPtCut1 = 20.; muonPtCut2 = 10.;
        electronPtCut1 = 25.; electronPtCut2 = 15.;
        emuPtCut1 = 25.; emuPtCut2 = 10.;
        muePtCut1 = 25.; muePtCut2 = 15.;
    }
    else if(DataYear==2018){
        muonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
        electronTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
        emuTriggers.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
        mueTriggers.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
        muonPtCut1 = 20.; muonPtCut2 = 10.;
        electronPtCut1 = 25.; electronPtCut2 = 15.;
        emuPtCut1 = 25.; emuPtCut2 = 10.;
        muePtCut1 = 25.; muePtCut2 = 15.;
    }

    jetPtCut = 20.;

}


ttbar_analysis::~ttbar_analysis(){

}

void ttbar_analysis::executeEvent(){

    allMuons = GetAllMuons();
    allElectrons = GetAllElectrons();
    allJets = GetAllJets();

    AnalyzerParameter param; param.Clear();

    param.syst_ = AnalyzerParameter::Central;

    param.Muon_Tight_ID = muonTightIDs.at(0);
    param.Muon_Loose_ID = muonLooseIDs.at(0);
    param.Electron_Tight_ID = electronTightIDs.at(0);
    param.Electron_Loose_ID = electronLooseIDs.at(0);
    param.Jet_ID = jetIDs.at(0);

    executeEventFromParameter(param);

}

void ttbar_analysis::executeEventFromParameter(AnalyzerParameter param){

    if(!PassMETFilter()) return;

    Event event = GetEvent();

    bool isDoubleMuon = false;
    bool isDoubleEG = false;
    bool isMuonEG = false;

    std::vector<Muon> muonsLoose = SelectMuons(allMuons, param.Muon_Loose_ID, muonPtCut, 2.4);
    std::vector<Muon> muons = SelectMuons(muonsLoose, param.Muon_Tight_ID, muonPtCut, 2.4);
    std::sort(muons.begin(), muons.end(), PtComparing);
    std::vector<Electron> electronsLoose = SelectElectrons(allElectrons, param.Electron_Loose_ID, electronPtCut, 2.5);
    std::vector<Electron> electrons = SelectElectrons(electronsLoose, param.Electron_Tight_ID, electronPtCut, 2.5);
    std::sort(electrons.begin(), electrons.end(), PtComparing);
    std::vector<Jet> jetsLoose = SelectJets(allJets, param.Jet_ID, jetPtCut, 2.5);
    std::vector<Jet> jets = JetsVetoLeptonInside(jetsLoose, electronsLoose, muonsLoose, 0.4);

    if (!(event.PassTrigger(muonTriggers) 
          || event.PassTrigger(electronTriggers)
          || event.PassTrigger(emuTriggers)
          || event.PassTrigger(mueTriggers))) return;

    std::vector<Lepton> leptons; leptons.clear();
    for (unsigned int i=0; i < muons.size(); i++) leptons.push_back(muons.at(i));
    for (unsigned int i=0; i < electrons.size(); i++) leptons.push_back(electrons.at(i));
    std::sort(leptons.begin(), leptons.end(), PtComparing);

    if (event.PassTrigger(muonTriggers)){
        if (muons.size() >= 2){
            if (muons.at(0).Pt() > muonPtCut1 && muons.at(1).Pt() > muonPtCut2) isDoubleMuon = true;
        }
    }
    if (event.PassTrigger(electronTriggers)){
        if (electrons.size() >= 2){
            if (electrons.at(0).Pt() > electronPtCut1 && electrons.at(1).Pt() > electronPtCut2) isDoubleEG = true;
        }
    }
    if (event.PassTrigger(emuTriggers)){
        if (muons.size() >= 1 && electrons.size() >= 1){
            if (electrons.at(0).Pt() > emuPtCut1 && muons.at(0).Pt() > emuPtCut2) isMuonEG = true;
        }
    }
    if (event.PassTrigger(mueTriggers)){
        if (muons.size() >= 1 && electrons.size() >= 1){
            if (muons.at(0).Pt() > muePtCut1 && electrons.at(0).Pt() > muePtCut2) isMuonEG = true;
        }
    }

    if (!(isDoubleMuon || isDoubleEG || isMuonEG)) return;

    if (IsDATA){
        if (DataStream.Contains("DoubleMuon")){
            if (!isDoubleMuon) return;
        }
        if (DataStream.Contains("DoubleEG")){
            if (isDoubleMuon) return; // priority goes to DoubleMuon dataset
            if (!isDoubleEG) return;
        }
        if (DataStream.Contains("MuonEG")){
            if (isDoubleMuon) return; // priority goes to DoubleMuon dataset
            if (isDoubleEG) return; // second priority goes to DoubleEG dataset
            if (!isMuonEG) return;
        }
    }

    double weight = 1.;
    if (!IsDATA){
        weight = weight * event.MCweight();
        weight = weight * GetPrefireWeight(0);
//        weight = weight * GetPileUpWeight(nPileUp,0); 
        weight = weight * weight_norm_1invpb * event.GetTriggerLumi("Full");
        if (HasFlag("reweightTopPt")){
            weight = weight * 0.;//work on this later
        }
    }

    TString leptonConfig;
    if (leptons.size() == 2){
        if (muons.size() == 2 && electrons.size() == 0 && isDoubleMuon) leptonConfig = "2l_2m0e";
        if (muons.size() == 1 && electrons.size() == 1 && isMuonEG) leptonConfig = "2l_1m1e";
        if (muons.size() == 0 && electrons.size() == 2 && isDoubleEG) leptonConfig = "2l_0m2e";
    }
    else if (leptons.size() == 3){
        if (muons.size() == 3 && electrons.size() == 0 && isDoubleMuon) leptonConfig = "3l_3m0e";
        if (muons.size() == 2 && electrons.size() == 1 && (isDoubleMuon || isMuonEG)) leptonConfig = "3l_2m1e";
        if (muons.size() == 1 && electrons.size() == 2 && (isMuonEG || isDoubleEG)) leptonConfig = "3l_1m2e";
        if (muons.size() == 0 && electrons.size() == 3 && isDoubleEG) leptonConfig = "3l_0m3e";
    }
    else if (leptons.size() == 4){
        if ((muons.size() == 4 && electrons.size() == 0 && isDoubleMuon)
             || (muons.size() == 3 && electrons.size() == 1 && (isDoubleMuon || isMuonEG))
             || (muons.size() == 2 && electrons.size() == 2 && (isDoubleMuon || isMuonEG || isDoubleEG))
             || (muons.size() == 1 && electrons.size() == 3 && (isMuonEG || isDoubleEG))
             || (muons.size() == 0 && electrons.size() == 4 && isDoubleEG)) leptonConfig = "4l_xmxe";
    }
    else leptonConfig = "NULL";

    leptonConfig = leptonConfig + "_";

    if (leptonConfig == "2l_2m0e_" || leptonConfig == "2l_1m1e_" || leptonConfig == "2l_0m2e_"){

        FillHist(leptonConfig + "nevents_weighted", 0., weight, 1, 0., 1.);
        FillHist(leptonConfig + "nevents_unweighted", 0., 1., 1, 0., 1.);

        FillHist(leptonConfig + "ptl1", leptons.at(0).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl2", leptons.at(1).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl", leptons.at(0).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl", leptons.at(1).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "etal1", leptons.at(0).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal2", leptons.at(1).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal", leptons.at(0).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal", leptons.at(1).Eta(), weight, 40, -5., 5.);

        if (leptonConfig != "2l_1m1e_"){
            TLorentzVector zCandidate = FindZCandidate(leptons);
            FillHist(leptonConfig + "masszcand", zCandidate.M(), weight, 400, 0., 400.);
        }
    }
    else if (leptonConfig == "3l_3m0e_" || leptonConfig == "3l_2m1e_" || leptonConfig == "3l_1m2e_" || leptonConfig == "3l_0m3e_"){

        FillHist(leptonConfig + "nevents_weighted", 0., weight, 1, 0., 1.);
        FillHist(leptonConfig + "nevents_unweighted", 0., 1., 1, 0., 1.);

        FillHist(leptonConfig + "ptl1", leptons.at(0).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl2", leptons.at(1).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl3", leptons.at(2).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl", leptons.at(0).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl", leptons.at(1).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl", leptons.at(2).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "etal1", leptons.at(0).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal2", leptons.at(1).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal3", leptons.at(2).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal", leptons.at(0).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal", leptons.at(1).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal", leptons.at(2).Eta(), weight, 40, -5., 5.);

        if (leptonConfig == "3l_3m0e_" || leptonConfig == "3l_2m1e_"){
            TLorentzVector zCandidate = FindZCandidate(muons);
            FillHist(leptonConfig + "masszcand", zCandidate.M(), weight, 400, 0., 400.);
        }
        if (leptonConfig == "3l_1m2e_" || leptonConfig == "3l_0m3e_"){
            TLorentzVector zCandidate = FindZCandidate(electrons);
            FillHist(leptonConfig + "masszcand", zCandidate.M(), weight, 400, 0., 400.);
        }

    }
    else if (leptonConfig != "NULL"){

        FillHist(leptonConfig + "nevents_weighted", 0., weight, 1, 0., 1.);
        FillHist(leptonConfig + "nevents_unweighted", 0., 1., 1, 0., 1.);

        FillHist(leptonConfig + "ptl1", leptons.at(0).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl2", leptons.at(1).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl3", leptons.at(2).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl4", leptons.at(3).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl", leptons.at(0).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl", leptons.at(1).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl", leptons.at(2).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "ptl", leptons.at(3).Pt(), weight, 400, 0., 400.);
        FillHist(leptonConfig + "etal1", leptons.at(0).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal2", leptons.at(1).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal3", leptons.at(2).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal4", leptons.at(3).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal", leptons.at(0).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal", leptons.at(1).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal", leptons.at(2).Eta(), weight, 40, -5., 5.);
        FillHist(leptonConfig + "etal", leptons.at(3).Eta(), weight, 40, -5., 5.);

    }






}

TLorentzVector ttbar_analysis::FindZCandidate(std::vector<Muon> leptons){

    TLorentzVector zCandidate;
    if (leptons.size() <= 1) zCandidate.SetPxPyPzE(0.,0.,0.,0.);
    else if (leptons.size() == 2) zCandidate = leptons.at(0) + leptons.at(1);
    else{
        zCandidate = leptons.at(0) + leptons.at(1);
        for (unsigned int i = 0; i < leptons.size(); i++){
            for (unsigned int j = i+1 ; j < leptons.size(); j++){
                TLorentzVector TMPzCandidate = leptons.at(i) + leptons.at(j);
                if (fabs(zCandidate.M() - 91.1876) > fabs(TMPzCandidate.M() - 91.1876)){
                    zCandidate = TMPzCandidate;
                }
            }
        }
    }
    return zCandidate;

}

TLorentzVector ttbar_analysis::FindZCandidate(std::vector<Electron> leptons){

    TLorentzVector zCandidate;
    if (leptons.size() <= 1) zCandidate.SetPxPyPzE(0.,0.,0.,0.);
    else if (leptons.size() == 2) zCandidate = leptons.at(0) + leptons.at(1);
    else{
        zCandidate = leptons.at(0) + leptons.at(1);
        for (unsigned int i = 0; i < leptons.size(); i++){
            for (unsigned int j = i+1 ; j < leptons.size(); j++){
                TLorentzVector TMPzCandidate = leptons.at(i) + leptons.at(j);
                if (fabs(zCandidate.M() - 91.1876) > fabs(TMPzCandidate.M() - 91.1876)){
                    zCandidate = TMPzCandidate;
                }
            }
        }
    }
    return zCandidate;

}

TLorentzVector ttbar_analysis::FindZCandidate(std::vector<Lepton> leptons){

    TLorentzVector zCandidate;
    if (leptons.size() <= 1) zCandidate.SetPxPyPzE(0.,0.,0.,0.);
    else if (leptons.size() == 2) zCandidate = leptons.at(0) + leptons.at(1);
    else{
        zCandidate = leptons.at(0) + leptons.at(1);
        for (unsigned int i = 0; i < leptons.size(); i++){
            for (unsigned int j = i+1 ; j < leptons.size(); j++){
                TLorentzVector TMPzCandidate = leptons.at(i) + leptons.at(j);
                if (fabs(zCandidate.M() - 91.1876) > fabs(TMPzCandidate.M() - 91.1876)){
                    zCandidate = TMPzCandidate;
                }
            }
        }
    }
    return zCandidate;

}
