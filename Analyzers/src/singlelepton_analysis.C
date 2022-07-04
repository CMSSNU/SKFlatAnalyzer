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

    muonTightIDs.push_back("POGLoose");
    muonLooseIDs.push_back("POGMedium");

//    muonTightIDs.push_back("POGTight");
//    muonLooseIDs.push_back("POGMedium");

//    electronTightIDs.push_back("passMVAID_noIso_WP80");
//    electronLooseIDs.push_back("passMVAID_noIso_WP90");

    electronTightIDs.push_back("passMediumID");
    electronLooseIDs.push_back("passLooseID");

    fatjetIDs.push_back("tight");
    jetIDs.push_back("tightLepVeto");

    muonTriggers.clear();
    electronTriggers.clear();

    if (DataYear == 2016){
        muonTriggers.push_back("HLT_IsoMu24_v");
//        muonTriggers.push_back("HLT_Mu50_v");
//        muonTriggers.push_back("HLT_TkMu50_v");
        electronTriggers.push_back("HLT_Ele27_WPTight_Gsf_v");
        electronTriggers.push_back("HLT_Photon175_v");
        electronTriggers.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
        muonPtCut = 40.;
        electronPtCut = 40.;
        leptonPtCut = 20.;
    }
    else if (DataYear == 2017){
        muonTriggers.push_back("HLT_IsoMu27_v");
//        muonTriggers.push_back("HLT_Mu50_v");
//        muonTriggers.push_back("HLT_OldMu100_v");
//        muonTriggers.push_back("HLT_TkMu100_v");
        electronTriggers.push_back("HLT_Ele35_WPTight_Gsf_v");
        electronTriggers.push_back("HLT_Photon200_v");
        electronTriggers.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
        muonPtCut = 40.;
        electronPtCut = 40.;
        leptonPtCut = 20.;
    }
    else if(DataYear==2018){
        muonTriggers.push_back("HLT_IsoMu24_v");
//        muonTriggers.push_back("HLT_Mu50_v");
//        muonTriggers.push_back("HLT_OldMu100_v");
//        muonTriggers.push_back("HLT_TkMu100_v");
        electronTriggers.push_back("HLT_Ele32_WPTight_Gsf_v");
        electronTriggers.push_back("HLT_Photon200_v");
        electronTriggers.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
        muonPtCut = 40.;
        electronPtCut = 40.;
        leptonPtCut = 20.;
    }

    jetPtCut = 30.; fatjetPtCut = 200.;

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

    allMuons = GetAllMuons();
//    allMuons = UseTunePMuon( GetAllMuons() );
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

    param.Muon_ID_SF_Key = "NUM_MediumID_DEN_TrackerMuons";
    param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_MediumID";
    param.Muon_Trigger_SF_Key = "IsoMu27_POGTight";

    param.Electron_ID_SF_Key = "passMediumID";
    param.Electron_Trigger_SF_Key = "";

    executeEventFromParameter(param);

}

void singlelepton_analysis::executeEventFromParameter(AnalyzerParameter param){

    if(!PassMETFilter()) return;

    Event event = GetEvent();

    if (!(event.PassTrigger(muonTriggers) || event.PassTrigger(electronTriggers))) return;

    // define leptons
    std::vector<Muon> muonsLoose = SelectMuons(allMuons, param.Muon_Loose_ID, leptonPtCut, 2.1);
    std::vector<Muon> muonsNoIso = SelectMuons(muonsLoose, param.Muon_Tight_ID, leptonPtCut, 2.1);
    std::vector<Muon> muons;
    for (unsigned int i=0 ; i < muonsNoIso.size(); i++){
        if (muonsNoIso.at(i).RelIso() < 0.15){
            muons.push_back(muonsNoIso.at(i));
        }
    }
    std::sort(muons.begin(), muons.end(), PtComparing);

    std::vector<Electron> electronsLoose = SelectElectrons(allElectrons, param.Electron_Loose_ID, leptonPtCut, 2.1);
    std::vector<Electron> electrons = SelectElectrons(electronsLoose, param.Electron_Tight_ID, leptonPtCut, 2.1);
    std::sort(electrons.begin(), electrons.end(), PtComparing);

    // define jets ak8
    std::vector<FatJet> fatjetsNoSDMass = SelectFatJets(allFatJets, param.FatJet_ID, fatjetPtCut, 2.1);
    fatjetsNoSDMass = FatJetsVetoLeptonInside(fatjetsNoSDMass, electronsLoose, muonsLoose, 0.4);
    std::sort(fatjetsNoSDMass.begin(), fatjetsNoSDMass.end(), PtComparing);
    std::vector<FatJet> fatjets;

    std::vector<FatJet> nonbbfatjets, bbfatjets;
    for (unsigned int i=0; i < fatjetsNoSDMass.size(); i++){
        if(fatjetsNoSDMass.at(i).SDMass() > 50.){
//            if(mcCorr->IsBTagged_2a(bTaggingMediumWP, fatjets.at(i))) bbfatjets.push_back(fatjetsNoSDMass.at(i));
//            else nonbbfatjets.push_back(fatjetsNoSDMass.at(i));
            fatjets.push_back(fatjetsNoSDMass.at(i));
        }
    }

    // define jets ak4
    std::vector<Jet> jets = SelectJets(allJets, param.Jet_ID, jetPtCut, 2.1);
    jets = JetsAwayFromFatJet(jets, fatjets, 1.2);

    //std::vector<Jet> jets = SelectJetsPileupMVA(JetsVetoLeptonInside(jetsLoose, electronsLoose, muonsLoose, 0.4), "loose");
    jets = JetsVetoLeptonInside(jets, electronsLoose, muonsLoose, 0.4);
    std::sort(jets.begin(), jets.end(), PtComparing);

    std::vector<Jet> nonbjets, bjets; // FIXME would it be the first two leading bjets the signal
    for (unsigned int i=0; i < jets.size(); i++){
        if( jets.at(i).GetTaggerResult(bTaggingMediumWP.j_Tagger) >= mcCorr->GetJetTaggingCutValue(bTaggingMediumWP.j_Tagger, bTaggingMediumWP.j_WP) ) bjets.push_back(jets.at(i));
        else nonbjets.push_back(jets.at(i));
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
            double muon_pt_weight = muons.at(i).MiniAODPt() > 120. ? 119.9 :muons.at(i).MiniAODPt();
            weight = weight * mcCorr->MuonID_SF(param.Muon_ID_SF_Key, muons.at(i).Eta(), muon_pt_weight, 0.);
            weight = weight * mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muon_pt_weight, 0.);
        }

        if(MCSample.Contains("TT") && MCSample.Contains("powheg")){
            std::vector<Gen> gens = GetGens();
            double top_pt_weight = mcCorr->TopPtReweight(gens);
         //   weight = weight * top_pt_weight;
            FillHist("top_pt_weight", top_pt_weight, 1., 200, 0., 2.);
        }
    }

    // trigger settings
    bool passMuonTrigger = (event.PassTrigger(muonTriggers) && muons.size()>= 1) ? muons.at(0).Pt() > muonPtCut : false;
    bool passElectronTrigger = (event.PassTrigger(electronTriggers) && electrons.size() >= 1) ? electrons.at(0).Pt() > electronPtCut : false;

    // lepton related selections
    bool hasOneMuon = (muons.size() == 1 && muonsLoose.size() == 1);
    bool hasOneElectron = (electrons.size() == 1 && electronsLoose.size() == 1);
    bool hasZeroMuon = (muons.size() == 0 && muonsLoose.size() == 0);
    bool hasZeroElectron = (electrons.size() == 0 && electronsLoose.size() == 0);

    // jet related selections
    bool hasZeroFatjet = (fatjets.size() == 0);
    bool hasAtLeastOneFatjet = (fatjets.size() >= 1);
    bool hasZeroBJet = (bjets.size() == 0);
    bool hasAtLeastTwoJet = (jets.size() >= 2);

    // met related selections
    bool hasMetAbove50 = (missingEt.Pt() > 50.);
    bool hasMetAbove100 = (missingEt.Pt() > 100.);

    // signal region specific selections
    bool hasSignalLepton = ((hasOneMuon && hasZeroElectron) || (hasOneElectron && hasZeroMuon));
    bool hasSignalJet = (hasAtLeastOneFatjet || (hasZeroFatjet && hasAtLeastTwoJet));
    double signalSecondaryBosonMass = -1.;
    Particle signalSecondaryBoson;
    if (hasSignalJet){
        if (hasAtLeastOneFatjet){
            signalSecondaryBosonMass = fatjets.at(0).SDMass();
            signalSecondaryBoson = fatjets.at(0);
        }
        else{
            signalSecondaryBosonMass = (jets.at(0) + jets.at(1)).M();
            signalSecondaryBoson = jets.at(0) + jets.at(1);
        }
    }

    bool hasSignalSecondaryBosonMass = (signalSecondaryBosonMass >= 0.) ? ((signalSecondaryBosonMass >= 70.) && (signalSecondaryBosonMass <= 145.)) : false;
    double bTagWP = mcCorr->GetJetTaggingCutValue(bTaggingMediumWP.j_Tagger, bTaggingMediumWP.j_WP);
    bool hasSignalBJet = hasAtLeastTwoJet ? ((jets.at(0).GetTaggerResult(bTaggingMediumWP.j_Tagger) >= bTagWP) && (jets.at(1).GetTaggerResult(bTaggingMediumWP.j_Tagger) >= bTagWP)) : false;

    Lepton signalLepton;
    Particle signalNeutrino;
    int signalNeutrinoDet = -1.;
    Particle signalPrimaryBoson;
    if (leptons.size() > 0) {
        signalLepton = leptons.at(0);
        signalNeutrino = GetReconstructedNeutrino(signalLepton, missingEt);
        signalNeutrinoDet = GetReconstructedNeutrinoDet(signalLepton, missingEt);
        if (hasSignalJet){
            signalPrimaryBoson = signalLepton + missingEt + signalSecondaryBoson;
        }
    }
    bool hasMtLeptonMissingEtAbove150 = (leptons.size() > 0) ? ((signalLepton + missingEt).Mt() > 150) : false;
    bool hasSignalObject = hasSignalLepton && hasSignalJet;

    bool hasImaginarySolution = (signalNeutrinoDet == 0);
    bool hasSignalBoosted = hasSignalJet && (hasAtLeastOneFatjet && hasZeroBJet) && hasSignalSecondaryBosonMass;
    bool hasSignalResolved = hasSignalJet && (hasZeroFatjet && hasAtLeastTwoJet && hasSignalBJet) && hasSignalSecondaryBosonMass;

    std::map<TString, bool> eventRegions;

    eventRegions["Signal_MuonBoostedSR1"] = passMuonTrigger && hasOneMuon && hasZeroElectron && hasSignalBoosted && hasMtLeptonMissingEtAbove150;
    eventRegions["Signal_ElectronBoostedSR1"] = passElectronTrigger && hasZeroMuon && hasOneElectron && hasSignalBoosted && hasMtLeptonMissingEtAbove150;
    eventRegions["Signal_MuonResolvedSR1"] = passMuonTrigger && hasOneMuon && hasZeroElectron && hasSignalResolved && hasMtLeptonMissingEtAbove150;
    eventRegions["Signal_ElectronResolvedSR1"] = passElectronTrigger && hasZeroMuon && hasOneElectron && hasSignalResolved && hasMtLeptonMissingEtAbove150;

    eventRegions["Signal_MuonBoostedSR2"] = passMuonTrigger && hasOneMuon && hasZeroElectron && hasSignalBoosted && hasMtLeptonMissingEtAbove150 && hasImaginarySolution;
    eventRegions["Signal_ElectronBoostedSR2"] = passElectronTrigger && hasZeroMuon && hasOneElectron && hasSignalBoosted && hasMtLeptonMissingEtAbove150 && hasImaginarySolution;
    eventRegions["Signal_MuonResolvedSR2"] = passMuonTrigger && hasOneMuon && hasZeroElectron && hasSignalResolved && hasMtLeptonMissingEtAbove150 && hasImaginarySolution;
    eventRegions["Signal_ElectronResolvedSR2"] = passElectronTrigger && hasZeroMuon && hasOneElectron && hasSignalResolved && hasMtLeptonMissingEtAbove150 && hasImaginarySolution;

    eventRegions["Control_MuonBoostedCR1"] = passMuonTrigger && hasOneMuon && hasZeroElectron && (hasSignalJet && (hasAtLeastOneFatjet && !hasZeroBJet) && hasSignalSecondaryBosonMass) && hasMtLeptonMissingEtAbove150;
    eventRegions["Control_ElectronBoostedCR1"] = passElectronTrigger && hasZeroMuon && hasOneElectron && (hasSignalJet && (hasAtLeastOneFatjet && !hasZeroBJet) && hasSignalSecondaryBosonMass) && hasMtLeptonMissingEtAbove150;
    eventRegions["Control_MuonResolvedCR1"] = passMuonTrigger && hasOneMuon && hasZeroElectron && (hasSignalJet && (hasZeroFatjet && !hasAtLeastTwoJet && !hasSignalBJet) && hasSignalSecondaryBosonMass) && hasMtLeptonMissingEtAbove150;
    eventRegions["Control_ElectronResolvedCR1"] = passElectronTrigger && hasZeroMuon && hasOneElectron && (hasSignalJet && (hasZeroFatjet && !hasAtLeastTwoJet && !hasSignalBJet) && hasSignalSecondaryBosonMass) && hasMtLeptonMissingEtAbove150;

    eventRegions["Control_MuonBoostedCR2"] = passMuonTrigger && hasOneMuon && hasZeroElectron && (hasSignalJet && (hasAtLeastOneFatjet && hasZeroBJet) && !hasSignalSecondaryBosonMass) && hasMtLeptonMissingEtAbove150;
    eventRegions["Control_ElectronBoostedCR2"] = passElectronTrigger && hasZeroMuon && hasOneElectron && (hasSignalJet && (hasAtLeastOneFatjet && hasZeroBJet) && !hasSignalSecondaryBosonMass) && hasMtLeptonMissingEtAbove150;
    eventRegions["Control_MuonResolvedCR2"] = passMuonTrigger && hasOneMuon && hasZeroElectron && (hasSignalJet && (hasZeroFatjet && hasAtLeastTwoJet && hasSignalBJet) && !hasSignalSecondaryBosonMass) && hasMtLeptonMissingEtAbove150;
    eventRegions["Control_ElectronResolvedCR2"] = passElectronTrigger && hasZeroMuon && hasOneElectron && (hasSignalJet && (hasZeroFatjet && hasAtLeastTwoJet && hasSignalBJet) && !hasSignalSecondaryBosonMass) && hasMtLeptonMissingEtAbove150;

    eventRegions["Control_MuonBoostedCR3"] = passMuonTrigger && hasOneMuon && hasZeroElectron && (hasSignalJet && (hasAtLeastOneFatjet && hasZeroBJet) && hasSignalSecondaryBosonMass) && !hasMtLeptonMissingEtAbove150;
    eventRegions["Control_ElectronBoostedCR3"] = passElectronTrigger && hasZeroMuon && hasOneElectron && (hasSignalJet && (hasAtLeastOneFatjet && hasZeroBJet) && hasSignalSecondaryBosonMass) && !hasMtLeptonMissingEtAbove150;
    eventRegions["Control_MuonResolvedCR3"] = passMuonTrigger && hasOneMuon && hasZeroElectron && (hasSignalJet && (hasZeroFatjet && hasAtLeastTwoJet && hasSignalBJet) && hasSignalSecondaryBosonMass) && !hasMtLeptonMissingEtAbove150;
    eventRegions["Control_ElectronResolvedCR3"] = passElectronTrigger && hasZeroMuon && hasOneElectron && (hasSignalJet && (hasZeroFatjet && hasAtLeastTwoJet && hasSignalBJet) && hasSignalSecondaryBosonMass) && !hasMtLeptonMissingEtAbove150;


/*
    eventRegions["Control_MuonBoostedCR3"] = passMuonTrigger && hasOneMuon && hasOneElectron && isSignalBoosted && hasMtLeptonMissingEtAbove150 ? muons.at(0).Charge() != electrons.at(0).Charge() : false;
    eventRegions["Control_ElectronBoostedCR3"] = passElectronTrigger && hasOneMuon && hasOneElectron && isSignalBoosted && hasMtLeptonMissingEtAbove150 ? muons.at(0).Charge() != electrons.at(0).Charge() : false;
    eventRegions["Control_MuonResolvedCR3"] = passMuonTrigger && hasOneMuon && hasOneElectron && isSignalResolved && hasMtLeptonMissingEtAbove150 ? muons.at(0).Charge() != electrons.at(0).Charge() : false;
    eventRegions["Control_ElectronResolvedCR3"] = passElectronTrigger && hasOneMuon && hasOneElectron && isSignalResolved && hasMtLeptonMissingEtAbove150 ? muons.at(0).Charge() != electrons.at(0).Charge() : false;
*/

    std::map<TString, bool>::iterator itEventRegions;

    for (itEventRegions = eventRegions.begin(); itEventRegions != eventRegions.end(); ++itEventRegions){

        if (itEventRegions->second){
            TString eventRegion = itEventRegions->first;
            DrawObjectPlots(eventRegion, weight, muons, electrons, jets, nonbjets, bjets, fatjets, nonbbfatjets, bbfatjets, missingEt);

            if (hasSignalLepton){
                FillHist(eventRegion + "_masst_leptonmet", (signalLepton + missingEt).Mt(), weight, 5000, 0., 5000.);
                FillHist(eventRegion + "_dphi_leptonmet", signalLepton.DeltaPhi(missingEt), weight, 50, -5., 5.);
                FillHist(eventRegion + "_mass_recowboson", (signalLepton + signalNeutrino).M(), weight, 5000,0., 5000.);
                FillHist(eventRegion + "_pz_reconeutrino", signalNeutrino.Pz(), weight, 100, -100., 100.);
                FillHist(eventRegion + "_det_reconeutrino", signalNeutrinoDet, weight, 2, 0., 2.);
            }
            if (hasSignalObject){
                FillHist(eventRegion + "_mass_secondaryboson", signalSecondaryBosonMass, weight, 5000, 0., 5000.);
                FillHist(eventRegion + "_masst_primaryboson", signalPrimaryBoson.Mt(), weight, 5000, 0., 5000.);
            }
        }
    }


    return;
}

double singlelepton_analysis::GetReconstructedNeutrinoDet(Lepton lepton, Particle missingEt){

    double mW = 80.4;
    double A = mW*mW + 2*missingEt.Px()*lepton.Px() + 2*missingEt.Py()*lepton.Py();
    double a = 4*lepton.Pz()*lepton.Pz() - 4*lepton.E()*lepton.E();
    double b = 4*lepton.Pz()*A;
    double c = A*A - 4*lepton.E()*lepton.E()*missingEt.Pt()*missingEt.Pt();
    double det = b*b - 4*a*c;

    if (det < 0) return 0.;
    else return 1.;

}

Particle singlelepton_analysis::GetReconstructedNeutrino(Lepton lepton, Particle missingEt){

    double pz;
    double mW = 80.4;

    double A = mW*mW + 2*missingEt.Px()*lepton.Px() + 2*missingEt.Py()*lepton.Py();
    double a = 4*lepton.Pz()*lepton.Pz() - 4*lepton.E()*lepton.E();
    double b = 4*lepton.Pz()*A;
    double c = A*A - 4*lepton.E()*lepton.E()*missingEt.Pt()*missingEt.Pt();
    double det = b*b - 4*a*c;

    if (det < 0){
        pz = -b / (2*a);
    }
    else{
        if (TMath::Abs(-b + TMath::Sqrt(det)) < TMath::Abs(-b - TMath::Sqrt(det))) pz = (-b + TMath::Sqrt(det)) / (2*a);
        else pz = (-b - TMath::Sqrt(det)) / (2*a);
    }


    Particle neutrino;
    neutrino.SetPxPyPzE(missingEt.Px(), missingEt.Py(), pz, TMath::Sqrt(missingEt.E()*missingEt.E() + pz*pz));

    return neutrino;

}

void singlelepton_analysis::DrawObjectPlots(TString region, double weight, std::vector<Muon> muons, std::vector<Electron> electrons, std::vector<Jet> jets, std::vector<Jet> nonbjets, std::vector<Jet> bjets, std::vector<FatJet> fatjets, std::vector<FatJet> nonbbfatjets, std::vector<FatJet> bbfatjets, Particle missingEt){

    double lt = 0., ht = 0., st = 0., met = missingEt.Pt();

    for (unsigned int i = 0; i < muons.size(); i++){
        Lepton lepton = muons.at(i);
        lt += lepton.Pt();
        FillHist(region + "_pt_mu", lepton.Pt(), weight, 5000, 0., 5000.);
        FillHist(region + "_eta_mu", lepton.Eta(), weight, 50, -5., 5.);
        FillHist(region + "_reliso_mu", lepton.RelIso(), weight, 50, 0., 1.);
        FillHist(region + "_minireliso_mu", lepton.MiniRelIso(), weight, 50, 0., 1.);
        FillHist(region + "_pt_lep", lepton.Pt(), weight, 5000, 0., 5000.);
        FillHist(region + "_eta_lep", lepton.Eta(), weight, 50, -5., 5.);
        FillHist(region + "_reliso_lep", lepton.RelIso(), weight, 50, 0., 1.);
        FillHist(region + "_minireliso_lep", lepton.MiniRelIso(), weight, 50, 0., 1.);
    }
    for (unsigned int i = 0; i < electrons.size(); i++){
        Lepton lepton = electrons.at(i);
        lt += lepton.Pt();
        FillHist(region + "_pt_el", lepton.Pt(), weight, 5000, 0., 5000.);
        FillHist(region + "_eta_el", lepton.Eta(), weight, 50, -5., 5.);
        FillHist(region + "_reliso_el", lepton.RelIso(), weight, 50, 0., 1.);
        FillHist(region + "_minireliso_el", lepton.MiniRelIso(), weight, 50, 0., 1.);
        FillHist(region + "_pt_lep", lepton.Pt(), weight, 5000, 0., 5000.);
        FillHist(region + "_eta_lep", lepton.Eta(), weight, 50, -5., 5.);
        FillHist(region + "_reliso_lep", lepton.RelIso(), weight, 50, 0., 1.);
        FillHist(region + "_minireliso_lep", lepton.MiniRelIso(), weight, 50, 0., 1.);
    }

    for (unsigned int i = 0; i < jets.size(); i++){
        Jet jet = jets.at(i);
        ht += jet.Pt();
        FillHist(region + "_pt_jet", jet.Pt(), weight, 5000, 0., 5000.);
        FillHist(region + "_eta_jet", jet.Eta(), weight, 50, -5., 5.);
    }
    for (unsigned int i = 0; i < bjets.size(); i++){
        Jet jet = bjets.at(i);
        FillHist(region + "_pt_bjet", jet.Pt(), weight, 5000, 0., 5000.);
        FillHist(region + "_eta_bjet", jet.Eta(), weight, 50, -5., 5.);
    }
    for (unsigned int i = 0; i < nonbjets.size(); i++){
        Jet jet = nonbjets.at(i);
        FillHist(region + "_pt_nonbjet", jet.Pt(), weight, 5000, 0., 5000.);
        FillHist(region + "_eta_nonbjet", jet.Eta(), weight, 50, -5., 5.);
    }

    for (unsigned int i = 0; i < fatjets.size(); i++){
        FatJet fatjet = fatjets.at(i);
        ht += fatjet.Pt();
        FillHist(region + "_pt_fatjet", fatjet.Pt(), weight, 5000, 0., 5000.);
        FillHist(region + "_eta_fatjet", fatjet.Eta(), weight, 50, -5., 5.);
        FillHist(region + "_sdmass_fatjet", fatjet.SDMass(), weight, 5000, 0., 5000.);
    }
    for (unsigned int i = 0; i < bbfatjets.size(); i++){
        FatJet fatjet = bbfatjets.at(i);
        FillHist(region + "_pt_bbfatjet", fatjet.Pt(), weight, 5000, 0., 5000.);
        FillHist(region + "_eta_bbfatjet", fatjet.Eta(), weight, 50, -5., 5.);
        FillHist(region + "_sdmass_bbfatjet", fatjet.SDMass(), weight, 5000, 0., 5000.);
    }
    for (unsigned int i = 0; i < nonbbfatjets.size(); i++){
        FatJet fatjet = nonbbfatjets.at(i);
        FillHist(region + "_pt_nonbbfatjet", fatjet.Pt(), weight, 5000, 0., 5000.);
        FillHist(region + "_eta_nonbbfatjet", fatjet.Eta(), weight, 50, -5., 5.);
        FillHist(region + "_sdmass_nonbbfatjet", fatjet.SDMass(), weight, 5000, 0., 5000.);
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

    FillHist(region + "_ltmet", (lt + met), weight, 5000, 0., 5000.);
    FillHist(region + "_lt", lt, weight, 5000, 0., 5000.);
    FillHist(region + "_ht", ht, weight, 5000, 0., 5000.); 
    FillHist(region + "_met", met, weight, 5000, 0., 5000.);
    FillHist(region + "_st", (lt + ht + met), weight, 5000, 0., 5000.);

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
