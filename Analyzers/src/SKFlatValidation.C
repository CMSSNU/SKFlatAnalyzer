#include "SKFlatValidation.h"

void SKFlatValidation::initializeAnalyzer(){


  //==== Triggers

  Triggers_POG_Electron.clear();
  Triggers_POG_Muon.clear();
  Triggers_POGHighPt_Electron.clear();
  Triggers_POGHighPt_Muon.clear();

  if(DataYear==2016){

    //==== Normal POG ID
    Triggers_POG_Electron = {
      "HLT_Ele27_WPTight_Gsf_v",
    };
    Triggers_POG_Muon = {
      "HLT_IsoMu24_v",
      "HLT_IsoTkMu24_v",
    };
    TriggerNameForSF_POG_Electron = "Ele27_WPTight_Gsf";
    TriggerNameForSF_POG_Muon = "IsoMu24";
    TriggerSafePt_POG_Electron = 30.;
    TriggerSafePt_POG_Muon = 26.;

    //==== POG High pt ID
    Triggers_POGHighPt_Electron = {
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Photon175_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    };
    Triggers_POGHighPt_Muon = {
      "HLT_Mu50",
      "HLT_TkMu50_v",
    };
    TriggerNameForSF_POGHighPt_Electron = "WREGammaTrigger";
    TriggerNameForSF_POGHighPt_Muon = "Mu50";
    TriggerSafePt_POGHighPt_Electron = 30.;
    TriggerSafePt_POGHighPt_Muon = 52.;
    DoublePhotonSafePtCut = 65.;
  }
  else if(DataYear==2017){

    //==== Normal POG ID
    Triggers_POG_Electron = {
      "HLT_Ele35_WPTight_Gsf_v",
    };
    Triggers_POG_Muon = {
      "HLT_IsoMu27_v",
    };
    TriggerNameForSF_POG_Electron = "Ele35_WPTight_Gsf";
    TriggerNameForSF_POG_Muon = "IsoMu27";
    TriggerSafePt_POG_Electron = 38.;
    TriggerSafePt_POG_Muon = 29.;

    //==== POG High pt ID
    Triggers_POGHighPt_Electron = {
      "HLT_Ele35_WPTight_Gsf_v",
      "HLT_Photon200_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    };
    Triggers_POGHighPt_Muon = {
      "HLT_Mu50_v",
      "HLT_OldMu100_v",
      "HLT_TkMu100_v",
    };
    TriggerNameForSF_POGHighPt_Electron = "WREGammaTrigger";
    TriggerNameForSF_POGHighPt_Muon = "Mu50";
    TriggerSafePt_POGHighPt_Electron = 38.;
    TriggerSafePt_POGHighPt_Muon = 52.;
    DoublePhotonSafePtCut = 75.;
  }
  else if(DataYear==2018){

    //==== Normal POG ID
    Triggers_POG_Electron = {
      "HLT_Ele32_WPTight_Gsf_v",
    };
    Triggers_POG_Muon = {
      "HLT_IsoMu24_v",
    };
    TriggerNameForSF_POG_Electron = "Ele32_WPTight_Gsf";
    TriggerNameForSF_POG_Muon = "IsoMu24";
    TriggerSafePt_POG_Electron = 35.;
    TriggerSafePt_POG_Muon = 26.;

    //==== POG High pt ID
    Triggers_POGHighPt_Electron = {
      "HLT_Ele32_WPTight_Gsf_v",
      "HLT_Photon200_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    };
    Triggers_POGHighPt_Muon = {
      "HLT_Mu50_v",
      "HLT_OldMu100_v",
      "HLT_TkMu100_v",
    };
    TriggerNameForSF_POGHighPt_Electron = "WREGammaTrigger";
    TriggerNameForSF_POGHighPt_Muon = "Mu50";
    TriggerSafePt_POGHighPt_Electron = 35.;
    TriggerSafePt_POGHighPt_Muon = 52.;
    DoublePhotonSafePtCut = 75.;
  }

  else{

  }

  //==== B-Tagging

  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a and 2a method
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== If you want to use 1d, which is a reshaping method
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::iterativefit, JetTagging::iterativefit) );
  //==== now put temporary vector (jtps) to MCCorrection::jetTaggingPars
  //==== jetTaggingPars will be looped over and used to read csv files when initializeAnalyzerTools() is ran.
  mcCorr->SetJetTaggingParameters(jtps);

}

void SKFlatValidation::executeEvent(){

  //==== Found one strange event from /gv0/DATA/SKFlat/Run2Legacy_v4/2016/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/191229_034335/0000/SKFlatNtuple_2016_MC_541.root
  //==== (run,event,lumi) = (1,111148613,240842)
  //==== - This event has "nan" for all pfMET_* variables
  //==== - I found it contains a jet with "(Pt, Eta, Phi, M, Charge) = 162845	-2.38395	2.04845	250.096	0"
  //==== - Let's not use this event.. it gives stderr when calculating MT()
  //==== - The PR that fixes this issue was not included in CMSSW_10_2_18
  //====   - https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1969/1.html
  //====   - https://github.com/cms-sw/cmssw/pull/27988/files
  if(pfMET_Type1_pt!=pfMET_Type1_pt) return;

  //==== Prefire reweight

  weight_Prefire = GetPrefireWeight(0);

  //==== AnalyzerParameter

  AnalyzerParameter param;

  //==== POG IDs

  param.Name = "POG";

  param.Electron_Tight_ID = "passMediumID";
  param.Electron_ID_SF_Key = "passMediumID";
  param.Electron_Trigger_SF_Key = "Default";

  param.Muon_Tight_ID = "POGTightWithTightIso";
  param.Muon_ID_SF_Key = "NUM_TightID_DEN_genTracks";
  param.Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";
  param.Muon_Trigger_SF_Key = "POGTight";

  param.Jet_ID = "HN";

  executeEventFromParameter(param);

  //==== POG HighPt

  param.Clear();

  param.Name = "POGHighPt";

  param.Electron_Tight_ID = "passHEEPID";
  param.Electron_ID_SF_Key = "HEEP";
  param.Electron_Trigger_SF_Key = "HEEP";

  param.Muon_Tight_ID = "POGHighPtWithLooseTrkIso";
  param.Muon_ID_SF_Key = "NUM_HighPtID_DEN_genTracks";
  param.Muon_ISO_SF_Key = "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut";
  param.Muon_Trigger_SF_Key = "POGHighPtLooseTrkIso";
  param.Muon_UseTuneP = true;

  param.Jet_ID = "HN";

  executeEventFromParameter(param);

}

void SKFlatValidation::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

  bool PassSingleElectron, PassSingleMuon;
  TString TriggerNameForSF_Electron, TriggerNameForSF_Muon;
  int TriggerSafePt_Electron, TriggerSafePt_Muon;
  double MinLeptonPt = 20.;
  if(param.Name=="POG"){
    TriggerNameForSF_Electron = TriggerNameForSF_POG_Electron;
    TriggerNameForSF_Muon = TriggerNameForSF_POG_Muon;
    TriggerSafePt_Electron = TriggerSafePt_POG_Electron;
    TriggerSafePt_Muon = TriggerSafePt_POG_Muon;
    PassSingleElectron = ev.PassTrigger(Triggers_POG_Electron);
    PassSingleMuon = ev.PassTrigger(Triggers_POG_Muon);
  }
  else if(param.Name=="POGHighPt"){
    TriggerNameForSF_Electron = TriggerNameForSF_POGHighPt_Electron;
    TriggerNameForSF_Muon = TriggerNameForSF_POGHighPt_Muon;
    TriggerSafePt_Electron = TriggerSafePt_POGHighPt_Electron;
    TriggerSafePt_Muon = TriggerSafePt_POGHighPt_Muon;
    PassSingleElectron = ev.PassTrigger(Triggers_POGHighPt_Electron);
    PassSingleMuon = ev.PassTrigger(Triggers_POGHighPt_Muon);
  }
  else{
    cout << "[SKFlatValidation::executeEventFromParameter] Wrong param.Name = " << param.Name << endl;
    exit(EXIT_FAILURE);
  }

  std::vector<Muon> muons = GetMuons(param.Muon_Tight_ID, MinLeptonPt, 2.4);
  if(param.Muon_UseTuneP) muons = UseTunePMuon(muons);
  std::vector<Electron> electrons = GetElectrons(param.Electron_Tight_ID, MinLeptonPt, 2.5);

  std::vector<Jet> myjets = JetsVetoLeptonInside( GetJets("tight", 30., 2.4), electrons, muons);
  int NBJets_NoSF(0), NBJets_WithSF_2a(0);
  double HT=0;

  JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV,
                                                                     JetTagging::Medium,
                                                                     JetTagging::incl, JetTagging::comb);
  JetTagging::Parameters jtp_DeepCSV_Medium_Reshape = JetTagging::Parameters(JetTagging::DeepCSV,
                                                                             JetTagging::Medium,
                                                                             JetTagging::iterativefit, JetTagging::iterativefit);

  //==== b tag SF; method 1a
  double btagWeight_1a = mcCorr->GetBTaggingReweight_1a(myjets, jtp_DeepCSV_Medium);
  double btagWeight_1d = mcCorr->GetBTaggingReweight_1d(myjets, jtp_DeepCSV_Medium_Reshape);

  for(unsigned int i=0; i<myjets.size(); i++){
    Jet this_jet = myjets.at(i);
    HT += this_jet.Pt();

    double this_discr = this_jet.GetTaggerResult(JetTagging::DeepCSV);

    //==== No SF
    if( this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium) ) NBJets_NoSF++;

    //==== 2a
    if( mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, this_jet) ) NBJets_WithSF_2a++;

  }

  //==== Based on which trigger is fired
  std::vector< TString > Suffixs = {
    "SingleMuon",
    "SingleElectron",
  };
  std::vector< bool > PassTriggers = {
    PassSingleMuon && (muons.size()>=1) && (electrons.size()==0),
    PassSingleElectron && (electrons.size()>=1) && (muons.size()==0),
  };

  for(unsigned int i=0; i<Suffixs.size(); i++){

    TString Suffix = Suffixs.at(i);
    if( !PassTriggers.at(i) ) continue;

    if(Suffix.Contains("SingleMuon")){
      if( muons.at(0).Pt() < TriggerSafePt_Muon ) continue;
    }
    else if(Suffix.Contains("SingleElectron")){
      if( electrons.at(0).Pt() < TriggerSafePt_Electron ) continue;
    }
    else{

    }

    if(this->IsDATA){
      if(this->DataStream == "SingleMuon"){
        if(!( Suffix.Contains("SingleMuon") )) continue;
      }
      else if(this->DataStream == "SingleElectron"){
        if(!( Suffix.Contains("SingleElectron") )) continue;
      }
      else{

      }
    }
    
    std::vector<Lepton *> leps;
    if(Suffix.Contains("SingleMuon")){
      leps = MakeLeptonPointerVector(muons);
    }
    else if(Suffix.Contains("SingleElectron")){
      leps = MakeLeptonPointerVector(electrons);
    }
    else{

    }

    Particle METv = ev.GetMETVector();

    int n_lepton = leps.size();
    //==== DiLepton variables
    bool IsOS = false;
    Particle Z;
    if(n_lepton==2){
      IsOS = (leps[0]->Charge() != leps[1]->Charge());
      Z = (*leps[0])+(*leps[1]);
    }
    //==== SingleLepton variables
    double this_MT=-999.;
    if(n_lepton==1){
      this_MT = MT( METv, (*leps[0]) );
    }

    double weight = 1.;
    if(!IsDATA){
      //cout << "GetTriggerLumi = " << ev.GetTriggerLumi("Full") << endl;
      //cout << "MCweight = " << MCweight() << endl;
      //cout << "weight_Prefire = " << weight_Prefire << endl;
      weight *= ev.GetTriggerLumi("Full")*MCweight()*weight_Prefire;

      mcCorr->IgnoreNoHist = param.MCCorrrectionIgnoreNoHist;

      //==== FIXME add third lepton veto later
      if(Suffix.Contains("SingleMuon")){
        for(unsigned int i=0; i<muons.size(); i++){

          double this_pt = muons.at(i).MiniAODPt();
          double this_eta = muons.at(i).Eta();

          double this_idsf  = mcCorr->MuonID_SF (param.Muon_ID_SF_Key,  this_eta, this_pt);
          double this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, this_eta, this_pt);
          double this_trigsf = mcCorr->MuonTrigger_SF(param.Muon_Trigger_SF_Key, TriggerNameForSF_Muon, muons);

          weight *= this_idsf*this_isosf*this_trigsf;

          if(param.Name=="POGHighPt"){
            double MiniAODP = sqrt( muons.at(i).MiniAODPt() * muons.at(i).MiniAODPt() + muons.at(i).Pz() * muons.at(i).Pz() );
            double this_recosf = mcCorr->MuonReco_SF("HighPtMuonRecoSF", this_eta, MiniAODP);
            weight *= this_recosf;
          }

        }
      }
      if( Suffix.Contains("DiElectron") || Suffix.Contains("SingleElectron") || Suffix.Contains("DiPhoton") ){
        for(unsigned int i=0; i<electrons.size(); i++){
          double this_recosf = mcCorr->ElectronReco_SF(electrons.at(i).scEta(),electrons.at(i).Pt());
          double this_idsf = mcCorr->ElectronID_SF(param.Electron_ID_SF_Key, electrons.at(i).scEta(), electrons.at(i).Pt());
          double this_trigsf = mcCorr->ElectronTrigger_SF(param.Electron_Trigger_SF_Key, TriggerNameForSF_Electron, electrons);

          weight *= this_recosf*this_idsf*this_trigsf;

        }
      }

    }

    std::map<TString, bool> map_bool_To_Region;
    if(n_lepton==2){

      bool BaseDiLepSelection = (Z.M() > 15.);
      bool HigherDiLeptonPtCut = (leps[0]->Pt() > DoublePhotonSafePtCut) && (leps[1]->Pt() > DoublePhotonSafePtCut);

      if(BaseDiLepSelection){

        if( IsOS ){
          //==== generic two OS lepton
          map_bool_To_Region["OS"] = true;
          //==== OnZ event
          map_bool_To_Region["OnZ_OS"] = IsOnZ(Z.M(), 15.);
          //==== High ZMass event
          map_bool_To_Region["HigherDiLeptonPtCut_OS"] = HigherDiLeptonPtCut;
          //==== OffZ, at least two jets, MET > 30 for dilepton ttbar
          map_bool_To_Region["mllgt110_TwoJets_METgt30_OS"] = (Z.M() > 110) && (myjets.size()>2) && (METv.Pt()>30.);
          //==== With B-jet, MET > 30 for dilepton ttbar
          map_bool_To_Region["WithBJet_METgt30_OS"] = (NBJets_NoSF>0) && (METv.Pt()>30.);
        }
        else{
          //==== generic two SS lepton
          map_bool_To_Region["SS"] = true;
          //==== OnZ event
          map_bool_To_Region["OnZ_SS"] = IsOnZ(Z.M(), 15.);
          //==== High ZMass event
          map_bool_To_Region["HigherDiLeptonPtCut_SS"] = HigherDiLeptonPtCut;
          //==== OffZ, at least two jets, MET > 30 for dilepton ttbar
          map_bool_To_Region["mllgt110_TwoJets_METgt30_SS"] = (Z.M() > 110) && (myjets.size()>2) && (METv.Pt()>30.);
          //==== With B-jet, MET > 30 for dilepton ttbar
          map_bool_To_Region["WithBJet_METgt30_SS"] = (NBJets_NoSF>0) && (METv.Pt()>30.);
        }

      }

    }
    if(n_lepton==1){
      map_bool_To_Region["W_CR"] = ( METv.Pt() > 30. ) && ( this_MT > 30. );
    }

    int nPUForReweight = nPV;
    if(!IsDATA) nPUForReweight = nPileUp;

    for(std::map<TString, bool>::iterator it_map = map_bool_To_Region.begin(); it_map != map_bool_To_Region.end(); it_map++){

      TString this_region = it_map->first;
      this_region = (param.Name)+"_"+Suffix+"_"+this_region;

      if(it_map->second){

        FillHist(this_region+"/NEvent_"+this_region, 0., weight, 1, 0., 1.);

        FillHist(this_region+"/nPileUp_"+this_region, nPileUp, weight, 200., 0., 200.);
        FillHist(this_region+"/nPV_"+this_region, nPV, weight, 200., 0., 200.);
        FillHist(this_region+"/nPUForReweight_"+this_region, nPUForReweight, weight, 200., 0., 200.);

        FillHist(this_region+"/MET_"+this_region, METv.Pt(), weight, 500, 0., 500.);
        FillHist(this_region+"/METphi_"+this_region, METv.Phi(), weight, 60, -3., 3.);

        FillHist(this_region+"/HT_"+this_region, HT, weight, 1000, 0., 1000.);

        FillHist(this_region+"/Jet_Size_"+this_region, myjets.size(), weight, 10, 0., 10.);
        FillHist(this_region+"/NBJets_NoSF_"+this_region, NBJets_NoSF, weight, 10, 0., 10.);
        FillHist(this_region+"/NBJets_WithSF_1a_"+this_region, NBJets_NoSF, weight*btagWeight_1a, 10, 0., 10.);
        FillHist(this_region+"/NBJets_WithSF_1d_"+this_region, NBJets_NoSF, weight*btagWeight_1d, 10, 0., 10.);
        FillHist(this_region+"/NBJets_WithSF_2a_"+this_region, NBJets_WithSF_2a, weight, 10, 0., 10.);

        for(unsigned int ij=0; ij<myjets.size(); ij++){
          TString this_itoa = TString::Itoa(ij,10);

          double this_discr = myjets.at(ij).GetTaggerResult(JetTagging::DeepCSV);
          FillHist(this_region+"/Jet_"+this_itoa+"_DeepCSV_"+this_region, this_discr, weight, 120, 0., 1.2);
          FillHist(this_region+"/Jet_"+this_itoa+"_DeepCSV_Scaled_"+this_region, this_discr, weight*btagWeight_1d, 120, 0., 1.2);
        }

        FillLeptonPlots(leps, this_region, weight);

        if(n_lepton==1){
          FillHist(this_region+"/MT_"+this_region, this_MT, weight, 500, 0., 500.);
        }

        if(n_lepton>=2){

          FillHist(this_region+"/ZCand_Mass_"+this_region, Z.M(), weight, 7000, 0., 7000.);
          FillHist(this_region+"/ZCand_Pt_"+this_region, Z.Pt(), weight, 500, 0., 500.);
          FillHist(this_region+"/ZCand_Eta_"+this_region, Z.Eta(), weight, 60, -3., 3.);

        }



      } // END if(pass Region)

    } // END Loop mapRegionBool

  } // END Loop Suffixs

}
