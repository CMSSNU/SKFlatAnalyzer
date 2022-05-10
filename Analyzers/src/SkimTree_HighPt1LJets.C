#include "SkimTree_HighPt1LJets.h"

void SkimTree_HighPt1LJets::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_HighPt1LJets::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  if(!IsDATA){
    newtree->SetBranchStatus("gen_*",0);
    newtree->SetBranchStatus("LHE_*",0);
    newtree->SetBranchStatus("gen_weight",1); // for MCweight()
  }

  triggers.clear();
  if(DataYear==2016){
    triggers = {
      "HLT_Mu50_v",
      "HLT_TkMu50_v",
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Photon175_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    };
  }
  else if(DataYear==2017){
    triggers = {
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

  LeptonPtCut = 40.;
  AK4JetPtCut = 25.;
  AK8JetPtCut = 170.;

  cout << "[SkimTree_HighPt1LJets::initializeAnalyzer] LeptonPtCut = " << LeptonPtCut << endl;
  cout << "[SkimTree_HighPt1LJets::initializeAnalyzer] AK8JetPtCut = " << AK8JetPtCut << endl;

}

void SkimTree_HighPt1LJets::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  //==== Skim 1 ) trigger
  if(! (ev.PassTrigger(triggers)) ) return;

  //==== Skim 2) at least two leptons (e or mu) with pt > "LeptonPtCut"

  vector<Muon> allmuons = UseTunePMuon( GetAllMuons() );
  int Nmuon = 0;
  for(unsigned int i=0; i<allmuons.size(); i++){
    if( allmuons.at(i).Pt() > LeptonPtCut ) Nmuon++;
  }
  vector<Electron> electrons = GetElectrons("NOCUT", LeptonPtCut, 2.5);
  int Nelectron = electrons.size();

  if (Nelectron + Nmuon == 0) return;

  //==== Skim 3) Jets

  vector<FatJet> allfatjets = puppiCorr->Correct( GetFatJets("tight", AK8JetPtCut, 2.4) ); //==== corret SDMass
  int Nfatjet_SDMassCut = 0;
  for(unsigned int i=0; i<allfatjets.size(); i++){
    if(allfatjets.at(i).SDMass()>40.) Nfatjet_SDMassCut++;
  }

  vector<Jet> jets = GetJets("tightLepVeto", AK4JetPtCut, 2.4);
  int Njet = jets.size();

  bool JetReq = (Njet>=2) || (Nfatjet_SDMassCut>=1);

  if ( !JetReq ) return;

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


