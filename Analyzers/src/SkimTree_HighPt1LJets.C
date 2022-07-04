#include "SkimTree_HighPt1LJets.h"

void SkimTree_HighPt1LJets::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_HighPt1LJets::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  /*
  if(!IsDATA){
    newtree->SetBranchStatus("gen_*",0);
    newtree->SetBranchStatus("LHE_*",0);
    newtree->SetBranchStatus("gen_weight",1); // for MCweight()
  }*/

  triggers.clear();

  if(DataYear==2016){
    triggers = {
      "HLT_IsoMu24_v",
      "HLT_Mu50_v",
      "HLT_TkMu50_v",
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Photon175_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    };
  }
  else if(DataYear==2017){
    triggers = {
      "HLT_IsoMu27_v",
      "HLT_Mu50_v",
      "HLT_OldMu100_v",
      "HLT_TkMu100_v",
      "HLT_Ele35_WPTight_Gsf_v",
      "HLT_Photon200_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    };
  }
  else if(DataYear==2018){
    triggers = {
      "HLT_IsoMu24_v",
      "HLT_Mu50_v",
      "HLT_OldMu100_v",
      "HLT_TkMu100_v",
      "HLT_Ele32_WPTight_Gsf_v",
      "HLT_Photon200_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    };
  }
  else{
    cout << "[SkimTree_HighPt1LJets::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_HighPt1LJets::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<triggers.size(); i++){
    cout << "[SkimTree_HighPt1LJets::initializeAnalyzer]   " << triggers.at(i) << endl;
  }

  LeptonPtCut = 18.;
  TunePPtCut = 40.;
  AK4JetPtCut = 20.;
  AK8JetPtCut = 170.;

}

void SkimTree_HighPt1LJets::executeEvent(){

    Event ev;
    ev.SetTrigger(*HLT_TriggerName);

    //==== Skim 1 ) trigger
    if(! (ev.PassTrigger(triggers)) ) return;

    //==== Skim 2) at least two leptons (e or mu) with pt > "LeptonPtCut"
    vector<Muon> muons = GetMuons("NOCUT", LeptonPtCut, 2.4);
    vector<Muon> tunep_muons = UseTunePMuon( GetAllMuons() );
    int n_tunep_muons = 0.;
    for(unsigned int i=0; i<tunep_muons.size(); i++){
        if( tunep_muons.at(i).Pt() > TunePPtCut ) n_tunep_muons++;
    }
    vector<Electron> electrons = GetElectrons("NOCUT", LeptonPtCut, 2.5);

    if (n_tunep_muons + muons.size() + electrons.size() == 0) return;

    //==== Skim 3) Jets
    vector<FatJet> fatjets = puppiCorr->Correct( GetFatJets("tight", AK8JetPtCut, 2.4) ); //==== corret SDMass
    int n_fatjets = 0;
    for(unsigned int i=0; i<fatjets.size(); i++){
        if(fatjets.at(i).SDMass()>20.) n_fatjets++;
    }

    vector<Jet> jets = GetJets("tightLepVeto", AK4JetPtCut, 2.4);

    if (!((n_fatjets >= 1) || (jets.size() >= 2))) return;

    //=============================
    //==== If survived, fill tree
    //=============================

    newtree->Fill();

}

void SkimTree_HighPt1LJets::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_HighPt1LJets::SkimTree_HighPt1LJets(){

  newtree = NULL;

}

SkimTree_HighPt1LJets::~SkimTree_HighPt1LJets(){

}

void SkimTree_HighPt1LJets::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}


