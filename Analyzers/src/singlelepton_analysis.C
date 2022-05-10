#include "singlelepton_analysis.h"

singlelepton_analysis::singlelepton_analysis(){

}

void singlelepton_analysis::initializeAnalyzer(){

    muonTightIDs.clear();
    muonLooseIDs.clear();
    electronTightIDs.clear();
    electronLooseIDs.clear();
    jetIDs.clear();
    fatjetIDs.clear();

    muonTightIDs.push_back("POGHighPtWithLooseTrkIso");
    muonLooseIDs.push_back("POGHighPt");

//    muonTightIDs.push_back("POGTight");
//    muonLooseIDs.push_back("POGMedium");

    electronTightIDs.push_back("passTightID");
    electronLooseIDs.push_back("passMediumID");

    fatjetIDs.push_back("tight");
    jetIDs.push_back("tightLepVeto");

    muonTriggers.clear();
    electronTriggers.clear();

    muonPtCut = 60.; electronPtCut = 60.; jetPtCut = 30.; fatjetPtCut = 200.;

    if (DataYear == 2016){
        muonTriggers.push_back("HLT_Mu50_v");
        muonTriggers.push_back("HLT_TkMu50_v");
        electronTriggers.push_back("HLT_Ele27_WPTight_Gsf_v");
        electronTriggers.push_back("HLT_Photon175_v");
        electronTriggers.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
    }
    else if (DataYear == 2017){
        muonTriggers.push_back("HLT_Mu50_v");
        muonTriggers.push_back("HLT_OldMu100_v");
        muonTriggers.push_back("HLT_TkMu100_v");
        electronTriggers.push_back("HLT_Ele35_WPTight_Gsf_v");
        electronTriggers.push_back("HLT_Photon200_v");
        electronTriggers.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
    }
    else if(DataYear==2018){
        muonTriggers.push_back("HLT_Mu50_v");
        muonTriggers.push_back("HLT_OldMu100_v");
        muonTriggers.push_back("HLT_TkMu100_v");
        electronTriggers.push_back("HLT_Ele32_WPTight_Gsf_v");
        electronTriggers.push_back("HLT_Photon200_v");
        electronTriggers.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
    }


//    std::vector<JetTagging::Parameters> bTaggings bbTaggings;
    bTaggingLooseWP = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Loose, JetTagging::incl, JetTagging::comb);
    bTaggingMediumWP = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
    bTaggingTightWP = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Tight, JetTagging::incl, JetTagging::comb);

    bTaggings.push_back(bTaggingLooseWP);
    bTaggings.push_back(bTaggingMediumWP);
    bTaggings.push_back(bTaggingTightWP);

    mcCorr->SetJetTaggingParameters(bTaggings);

}


singlelepton_analysis::~singlelepton_analysis(){

}

void singlelepton_analysis::executeEvent(){

    allMuons = UseTunePMuon( GetAllMuons() );
    allElectrons = GetAllElectrons();
    allJets = GetAllJets();
    allFatJets = puppiCorr->Correct(GetAllFatJets());

    AnalyzerParameter param; param.Clear();

    param.syst_ = AnalyzerParameter::Central;

    param.Muon_Tight_ID = muonTightIDs.at(0);
    param.Muon_Loose_ID = muonLooseIDs.at(0);
    param.Electron_Tight_ID = electronTightIDs.at(0);
    param.Electron_Loose_ID = electronLooseIDs.at(0);
    param.Jet_ID = jetIDs.at(0);
    param.FatJet_ID = fatjetIDs.at(0);

    param.Muon_RECO_SF_Key = "HighPtMuonRecoSF";
    param.Muon_ID_SF_Key = "NUM_HighPtID_DEN_TrackerMuons";
    param.Muon_ISO_SF_Key = "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut";
    param.Muon_Trigger_SF_Key = "Mu50_POGHighPtLooseTrkIso";
    param.Electron_ID_SF_Key = "passTightID";
    param.Electron_Trigger_SF_Key = "";

    executeEventFromParameter(param);

}

void singlelepton_analysis::executeEventFromParameter(AnalyzerParameter param){

    if(!PassMETFilter()) return;

    Event event = GetEvent();

    if (!(event.PassTrigger(muonTriggers) || event.PassTrigger(electronTriggers))) return;

    // define leptons
    std::vector<Muon> muonsLoose = SelectMuons(allMuons, param.Muon_Loose_ID, 20., 2.1);
    std::vector<Muon> muons = SelectMuons(muonsLoose, param.Muon_Tight_ID, muonPtCut, 2.1);
    std::sort(muons.begin(), muons.end(), PtComparing);

    std::vector<Electron> electronsLoose = SelectElectrons(allElectrons, param.Electron_Loose_ID, 20, 2.1);
    std::vector<Electron> electrons = SelectElectrons(electronsLoose, param.Electron_Tight_ID, electronPtCut, 2.1);
    std::sort(electrons.begin(), electrons.end(), PtComparing);

    // define jets ak4
    std::vector<Jet> jetsLoose = SelectJets(allJets, param.Jet_ID, jetPtCut, 2.1);
    //std::vector<Jet> jets = SelectJetsPileupMVA(JetsVetoLeptonInside(jetsLoose, electronsLoose, muonsLoose, 0.4), "loose");
    std::vector<Jet> jets = JetsVetoLeptonInside(jetsLoose, electronsLoose, muonsLoose, 0.4);
    std::sort(jets.begin(), jets.end(), PtComparing);

    std::vector<Jet> nonbjets, bjets; // FIXME would it be the first two leading bjets the signal
    for (unsigned int i=0; i < jets.size(); i++){
        if( jets.at(i).GetTaggerResult(bTaggingMediumWP.j_Tagger) <= mcCorr->GetJetTaggingCutValue(bTaggingMediumWP.j_Tagger, bTaggingMediumWP.j_WP) ) bjets.push_back(jets.at(i));
        else nonbjets.push_back(jets.at(i));
    }

    // define jets ak8
    std::vector<FatJet> fatjetsLoose = SelectFatJets(allFatJets, param.FatJet_ID, fatjetPtCut, 2.1);
    std::vector<FatJet> fatjetsRaw = FatJetsVetoLeptonInside(fatjetsLoose, electronsLoose, muonsLoose, 0.4);
    std::sort(fatjetsRaw.begin(), fatjetsRaw.end(), PtComparing);

    std::vector<FatJet> nonbbfatjets, bbfatjets, fatjets;
    for (unsigned int i=0; i < fatjetsRaw.size(); i++){
        if(fatjetsRaw.at(i).SDMass() > 40.){
            if(mcCorr->IsBTagged_2a(bTaggingMediumWP, fatjetsRaw.at(i))) bbfatjets.push_back(fatjetsRaw.at(i));
            else nonbbfatjets.push_back(fatjetsRaw.at(i));
            fatjets.push_back(fatjetsRaw.at(i));
        }
    }

    // define missing et
    Particle missingEt = event.GetMETVector();

    std::vector<Lepton> leptons; leptons.clear();
    for (unsigned int i=0; i < muons.size(); i++) leptons.push_back(muons.at(i));
    for (unsigned int i=0; i < electrons.size(); i++) leptons.push_back(electrons.at(i));
    std::sort(leptons.begin(), leptons.end(), PtComparing);

    double weight = 1.;
    if (!IsDATA){
        weight = weight * event.MCweight();
        weight = weight * GetPrefireWeight(0);
        weight = weight * GetPileUpWeight(nPileUp,0); 
        weight = weight * weight_norm_1invpb  * event.GetTriggerLumi("Full");
        weight = weight * mcCorr->GetBTaggingReweight_1a(jets, bTaggingMediumWP);
        for (unsigned int i=0; i < electrons.size(); i++){
            weight = weight * mcCorr->ElectronReco_SF(electrons.at(i).scEta(), electrons.at(i).UncorrPt(), 0.);
            weight = weight * mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(i).scEta(), electrons.at(i).Pt(), 0.);
            // FIXME no trigger SF
        }
        for (unsigned int i=0; i < muons.size(); i++){
//            weight = weight * mcCorr->MuonReco_SF(param.Muon_RECO_SF_Key, muons.at(i).Eta(), sqrt(muons.at(i).MiniAODPt()*muons.at(i).MiniAODPt() + muons.at(i).Pz()*muons.at(i).Pz()), 0.);
            weight = weight * mcCorr->MuonID_SF(param.Muon_ID_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0.);
            weight = weight * mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt(), 0.);
        }
    }

    bool hasLepton = ((event.PassTrigger(muonTriggers) && muons.size()>= 1) || (event.PassTrigger(electronTriggers) && electrons.size() >= 1));
    bool hasOneLepton = (hasLepton && leptons.size() == 1 && (muonsLoose.size() + electronsLoose.size() == 1));
    bool hasTwoLepton =  (hasLepton && leptons.size() == 2 && (muonsLoose.size() + electronsLoose.size() == 2));
    bool hasOneMuon = (hasOneLepton && muons.size() == 1);
    bool hasTwoMuon = (hasTwoLepton && muons.size() == 2) ? (muons.at(0).Charge() != muons.at(1).Charge()) && ((leptons.at(0) + leptons.at(1)).M() > 70) : false;
    bool hasOneElectron = (hasOneLepton && electrons.size() == 1);
    bool hasTwoElectron = (hasTwoLepton && electrons.size() == 2) ? (electrons.at(0).Charge() != electrons.at(1).Charge()) && ((leptons.at(0) + leptons.at(1)).M() > 70) : false;
    bool hasAtLeastOneFatjet = (fatjets.size() >= 1);
    bool hasAtLeastTwoJet = (jets.size() >= 2);
    bool hasAtLeastTwoBJet = (bjets.size() >= 2);
    bool hasTwoBJet = (bjets.size() == 2);

    // Preselection
    if (hasLepton && (hasAtLeastOneFatjet || hasAtLeastTwoBJet)){
        DrawObjectPlots("Preselection_Inclusive", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);
        if (hasOneMuon) DrawObjectPlots("Preselection_InclusiveMu", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);
        if (hasOneElectron) DrawObjectPlots("Preselection_InclusiveEl", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);
        if (hasAtLeastOneFatjet) DrawObjectPlots("Preselection_Boosted", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);
        if (hasAtLeastTwoBJet) DrawObjectPlots("Preselection_Resolved", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);

        if (hasOneMuon && hasAtLeastOneFatjet) DrawObjectPlots("Preselection_BoostedMu", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);
        if (hasOneElectron && hasAtLeastOneFatjet) DrawObjectPlots("Preselection_BoostedEl", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);
        if (hasOneMuon && hasAtLeastTwoBJet) DrawObjectPlots("Preselection_ResolvedMu", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);
        if (hasOneElectron && hasAtLeastTwoBJet) DrawObjectPlots("Preselection_ResolvedEl", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);

    }

    if (hasTwoMuon){
        TString region = "Preselection_Inclusive2Mu";
        DrawObjectPlots("Preselection_Inclusive2Mu", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);
        FillHist(region + "_mass_mumu", (leptons.at(0) +leptons.at(1)).M(), weight, 500, 0., 500.);
    }
    if (hasTwoElectron){
        TString region = "Preselection_Inclusive2El";
        DrawObjectPlots("Preselection_Inclusive2El", weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);
        FillHist(region + "_mass_elel", (leptons.at(0) +leptons.at(1)).M(), weight, 500, 0., 500.);
    }
//        DrawSignalPlots(region, weight, lepton, missingEt, recoSecondaryBoson, recoHeavyNeutrino, recoPrimaryBoson);


    return;
}

void singlelepton_analysis::DrawSignalPlots(TString region, double weight, Lepton lepton, Particle missingEt, Particle recoSecondaryBoson, Particle recoHeavyNeutrino, Particle recoPrimaryBoson){

    FillHist(region + "_mass_secondaryboson", recoSecondaryBoson.M(), weight, 500, 0., 500.);
    FillHist(region + "_masst_heavyneutrino", recoHeavyNeutrino.Mt(), weight, 2500, 0., 2500.);
    FillHist(region + "_masst_primaryboson", recoPrimaryBoson.Mt(), weight, 4000, 0., 4000.);
    FillHist(region + "_dphi_leptonmet", lepton.DeltaPhi(missingEt), weight, 40, -5., 5.);
    FillHist(region + "_dphi_leptonsecondaryboson", lepton.DeltaPhi(recoSecondaryBoson), weight, 40, -5., 5.);
    FillHist(region + "_dphi_leptonheavyneutrino", lepton.DeltaPhi(recoHeavyNeutrino), weight, 40, -5., 5.);
    FillHist(region + "_dphi_metsecondaryboson", missingEt.DeltaPhi(recoSecondaryBoson), weight, 40, -5., 5.);

}

void singlelepton_analysis::DrawObjectPlots(TString region, double weight, std::vector<Muon> muons, std::vector<Electron> electrons, std::vector<Jet> jets, std::vector<Jet> nonbjets, std::vector<Jet> bjets, std::vector<FatJet> fatjets, std::vector<FatJet> nonbbfatjets, std::vector<FatJet> bbfatjets, Particle missingEt){

    double lt = 0., ht = 0., st = 0., met = missingEt.Pt();

    for (unsigned int i = 0; i < muons.size(); i++){
        Lepton lepton = muons.at(i);
        lt += lepton.Pt();
        FillHist(region + "_pt_mu", lepton.Pt(), weight, 500, 0., 500.);
        FillHist(region + "_eta_mu", lepton.Eta(), weight, 40, -5., 5.);
        FillHist(region + "_reliso_mu", lepton.RelIso(), weight, 40, 0., 1.);
        FillHist(region + "_minireliso_mu", lepton.MiniRelIso(), weight, 40, 0., 1.);
        FillHist(region + "_pt_lep", lepton.Pt(), weight, 500, 0., 500.);
        FillHist(region + "_eta_lep", lepton.Eta(), weight, 40, -5., 5.);
        FillHist(region + "_reliso_lep", lepton.RelIso(), weight, 40, 0., 1.);
        FillHist(region + "_minireliso_lep", lepton.MiniRelIso(), weight, 40, 0., 1.);
    }
    for (unsigned int i = 0; i < electrons.size(); i++){
        Lepton lepton = electrons.at(i);
        lt += lepton.Pt();
        FillHist(region + "_pt_el", lepton.Pt(), weight, 500, 0., 500.);
        FillHist(region + "_eta_el", lepton.Eta(), weight, 40, -5., 5.);
        FillHist(region + "_reliso_el", lepton.RelIso(), weight, 40, 0., 1.);
        FillHist(region + "_minireliso_el", lepton.MiniRelIso(), weight, 40, 0., 1.);
        FillHist(region + "_pt_lep", lepton.Pt(), weight, 500, 0., 500.);
        FillHist(region + "_eta_lep", lepton.Eta(), weight, 40, -5., 5.);
        FillHist(region + "_reliso_lep", lepton.RelIso(), weight, 40, 0., 1.);
        FillHist(region + "_minireliso_lep", lepton.MiniRelIso(), weight, 40, 0., 1.);
    }

    for (unsigned int i = 0; i < jets.size(); i++){
        Jet jet = jets.at(i);
        ht += jet.Pt();
        FillHist(region + "_pt_jet", jet.Pt(), weight, 500, 0., 500.);
        FillHist(region + "_eta_jet", jet.Eta(), weight, 40, -5., 5.);
    }
    for (unsigned int i = 0; i < bjets.size(); i++){
        Jet jet = bjets.at(i);
        FillHist(region + "_pt_bjet", jet.Pt(), weight, 500, 0., 500.);
        FillHist(region + "_eta_bjet", jet.Eta(), weight, 40, -5., 5.);
    }
    for (unsigned int i = 0; i < nonbjets.size(); i++){
        Jet jet = nonbjets.at(i);
        FillHist(region + "_pt_nonbjet", jet.Pt(), weight, 500, 0., 500.);
        FillHist(region + "_eta_nonbjet", jet.Eta(), weight, 40, -5., 5.);
    }

    for (unsigned int i = 0; i < fatjets.size(); i++){
        FatJet fatjet = fatjets.at(i);
        ht += fatjet.Pt();
        FillHist(region + "_pt_fatjet", fatjet.Pt(), weight, 500, 0., 500.);
        FillHist(region + "_eta_fatjet", fatjet.Eta(), weight, 40, -5., 5.);
        FillHist(region + "_sdmass_fatjet", fatjet.SDMass(), weight, 500, 0., 500.);
    }
    for (unsigned int i = 0; i < bbfatjets.size(); i++){
        FatJet fatjet = bbfatjets.at(i);
        FillHist(region + "_pt_bbfatjet", fatjet.Pt(), weight, 500, 0., 500.);
        FillHist(region + "_eta_bbfatjet", fatjet.Eta(), weight, 40, -5., 5.);
        FillHist(region + "_sdmass_bbfatjet", fatjet.SDMass(), weight, 500, 0., 500.);
    }
    for (unsigned int i = 0; i < nonbbfatjets.size(); i++){
        FatJet fatjet = nonbbfatjets.at(i);
        FillHist(region + "_pt_nonbbfatjet", fatjet.Pt(), weight, 500, 0., 500.);
        FillHist(region + "_eta_nonbbfatjet", fatjet.Eta(), weight, 40, -5., 5.);
        FillHist(region + "_sdmass_nonbbfatjet", fatjet.SDMass(), weight, 500, 0., 500.);
    }

    FillHist(region + "_n_events_weighted", 0, weight, 1, 0., 1.);
    FillHist(region + "_n_events_unweighted", 0, 1., 1, 0., 1.);
    FillHist(region + "_n_lep", (muons.size() + electrons.size()), weight, 5, 0., 5.);
    FillHist(region + "_n_mu", muons.size(), weight, 5, 0., 5.);
    FillHist(region + "_n_el", electrons.size(), weight, 5, 0., 5.);
    FillHist(region + "_n_jet", jets.size(), weight, 5, 0., 5.);
    FillHist(region + "_n_bjet", bjets.size(), weight, 5, 0., 5.);
    FillHist(region + "_n_nonbjet", nonbjets.size(), weight, 5, 0., 5.);
    FillHist(region + "_n_fatjet", fatjets.size(), weight, 5, 0., 5.);
    FillHist(region + "_n_bbfatjet", bbfatjets.size(), weight, 5, 0., 5.);
    FillHist(region + "_n_nonbbfatjet", nonbbfatjets.size(), weight, 5, 0., 5.);

   FillHist(region + "_lt", lt, weight, 1000, 0., 1000.);
   FillHist(region + "_ht", ht, weight, 1000, 0., 1000.); 
   FillHist(region + "_met", met, weight, 1000, 0., 1000.);
   FillHist(region + "_st", (lt + ht + met), weight, 1000, 0., 1000.);

}

TLorentzVector singlelepton_analysis::FindZCandidate(std::vector<Muon> leptons){

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

TLorentzVector singlelepton_analysis::FindZCandidate(std::vector<Electron> leptons){

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

TLorentzVector singlelepton_analysis::FindZCandidate(std::vector<Lepton> leptons){

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
