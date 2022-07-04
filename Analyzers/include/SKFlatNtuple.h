#ifndef SKFlatNtuple_h
#define SKFlatNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

// Fixed size dimensions of array or collections stored in the TTree if any.

using namespace std;

class SKFlatNtuple {
public :

   SKFlatNtuple();
   virtual ~SKFlatNtuple();

   virtual Int_t GetEntry(Long64_t entry);

   virtual void SetTreeName(TString tname){
     fChain = new TChain(tname);
   }

   virtual int AddFile(TString filename){
     return fChain->Add(filename,-1);
   }

   Long64_t MaxEvent, NSkipEvent;
   int LogEvery;
   bool IsDATA;
   TString DataStream;
   TString MCSample;
   bool IsFastSim;
   int DataYear;
   TString DataEra;
   double xsec, sumW, sumSign;
   vector<TString> Userflags;

   virtual void Init();
   virtual void Loop();

   virtual void beginEvent(){};
   virtual void executeEvent(){};
   virtual void endEvent(){};

   virtual void SetEra(TString era){
     DataEra=era;
     DataYear=TString(era(0,4)).Atoi();
   }
   virtual TString GetEra() const { return DataEra; }
   virtual TString GetEraShort() const;
   virtual int GetYear() const { return DataYear; }

   std::string AddZeroToTime(int twodigit){
     if(twodigit<10){
       return "0"+std::to_string(twodigit);
     }
     else{
       return std::to_string(twodigit);
     }
   }

   std::string printcurrunttime(){

     std::stringstream out;
     TDatime datime;
     out << datime.GetYear()<<"-"<<AddZeroToTime(datime.GetMonth())<<"-"<<AddZeroToTime(datime.GetDay())<<" "<<AddZeroToTime(datime.GetHour())<<":"<<AddZeroToTime(datime.GetMinute())<<":"<<AddZeroToTime(datime.GetSecond());
     return out.str();

   }

   TChain *fChain=NULL;

   // Declaration of leaf types
   Bool_t          IsData;
   Int_t           nTotal;
   Int_t           run;
   ULong64_t       event;
   Int_t           lumi;
   Float_t         Rho;
   Float_t         RhoNC;
   Int_t           nPV;
   Float_t         L1PrefireReweight_Central;
   Float_t         L1PrefireReweight_Up;
   Float_t         L1PrefireReweight_Down;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadPFMuonDzFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Int_t           PVtrackSize;
   Float_t         PVchi2;
   Float_t         PVndof;
   Float_t         PVnormalizedChi2;
   Float_t         vertex_X;
   Float_t         vertex_Y;
   Float_t         vertex_Z;
   vector<string>  *HLT_TriggerName;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_charge;
   vector<float>   *jet_area;
   vector<int>     *jet_partonFlavour;
   vector<int>     *jet_hadronFlavour;
   vector<int>     *jet_GenHFHadronMatcher_flavour;
   vector<int>     *jet_GenHFHadronMatcher_origin;
   vector<float>   *jet_DeepCSV;
   vector<float>   *jet_DeepCSV_CvsL;
   vector<float>   *jet_DeepCSV_CvsB;
   vector<float>   *jet_DeepJet;
   vector<float>   *jet_DeepJet_CvsL;
   vector<float>   *jet_DeepJet_CvsB;
   vector<float>   *jet_chargedHadronEnergyFraction;
   vector<float>   *jet_neutralHadronEnergyFraction;
   vector<float>   *jet_neutralEmEnergyFraction;
   vector<float>   *jet_chargedEmEnergyFraction;
   vector<float>   *jet_muonEnergyFraction;
   vector<int>     *jet_chargedMultiplicity;
   vector<int>     *jet_neutralMultiplicity;
   vector<bool>    *jet_tightJetID;
   vector<bool>    *jet_tightLepVetoJetID;
   vector<int>     *jet_partonPdgId;
   vector<int>     *jet_vtxNtracks;
   vector<float>   *jet_m;
   vector<float>   *jet_energy;
   vector<float>   *jet_PileupJetId;
   vector<float>   *jet_shiftedEnUp;
   vector<float>   *jet_shiftedEnDown;
   vector<float>   *jet_smearedRes;
   vector<float>   *jet_smearedResUp;
   vector<float>   *jet_smearedResDown;
   vector<float>   *jet_bJetNN_corr;
   vector<float>   *jet_bJetNN_res;
   vector<float>   *jet_cJetNN_corr;
   vector<float>   *jet_cJetNN_res;
   vector<float>   *jet_JECL1FastJet;
   vector<float>   *jet_JECFull;
   vector<float>   *fatjet_pt;
   vector<float>   *fatjet_eta;
   vector<float>   *fatjet_phi;
   vector<float>   *fatjet_charge;
   vector<float>   *fatjet_area;
   vector<int>     *fatjet_partonFlavour;
   vector<int>     *fatjet_hadronFlavour;
   vector<float>   *fatjet_DeepCSV;
   vector<float>   *fatjet_DeepCSV_CvsL;
   vector<float>   *fatjet_DeepCSV_CvsB;
   vector<float>   *fatjet_particleNet_TvsQCD;
   vector<float>   *fatjet_particleNet_WvsQCD;
   vector<float>   *fatjet_particleNet_ZvsQCD;
   vector<float>   *fatjet_particleNet_HbbvsQCD;
   vector<float>   *fatjet_particleNet_HccvsQCD;
   vector<float>   *fatjet_particleNet_H4qvsQCD;
   vector<float>   *fatjet_particleNet_QCD;
   vector<float>   *fatjet_particleNetMD_Xbb;
   vector<float>   *fatjet_particleNetMD_Xcc;
   vector<float>   *fatjet_particleNetMD_Xqq;
   vector<float>   *fatjet_particleNetMD_QCD;
   vector<bool>    *fatjet_tightJetID;
   vector<bool>    *fatjet_tightLepVetoJetID;
   vector<int>     *fatjet_partonPdgId;
   vector<int>     *fatjet_vtxNtracks;
   vector<float>   *fatjet_m;
   vector<float>   *fatjet_energy;
   vector<float>   *fatjet_puppi_tau1;
   vector<float>   *fatjet_puppi_tau2;
   vector<float>   *fatjet_puppi_tau3;
   vector<float>   *fatjet_puppi_tau4;
   vector<float>   *fatjet_softdropmass;
   vector<float>   *fatjet_chargedHadronEnergyFraction;
   vector<float>   *fatjet_neutralHadronEnergyFraction;
   vector<float>   *fatjet_neutralEmEnergyFraction;
   vector<float>   *fatjet_chargedEmEnergyFraction;
   vector<float>   *fatjet_muonEnergyFraction;
   vector<int>     *fatjet_chargedMultiplicity;
   vector<int>     *fatjet_neutralMultiplicity;
   vector<float>   *fatjet_shiftedEnUp;
   vector<float>   *fatjet_shiftedEnDown;
   vector<float>   *fatjet_smearedRes;
   vector<float>   *fatjet_smearedResUp;
   vector<float>   *fatjet_smearedResDown;
   vector<float>   *fatjet_LSF;
   vector<float>   *fatjet_LSFlep_PID;
   vector<float>   *fatjet_LSFlep_Pt;
   vector<float>   *fatjet_LSFlep_Eta;
   vector<float>   *fatjet_LSFlep_Phi;
   vector<float>   *electron_MVAIso;
   vector<float>   *electron_MVANoIso;
   vector<float>   *electron_Energy;
   vector<float>   *electron_Energy_Scale_Up;
   vector<float>   *electron_Energy_Scale_Down;
   vector<float>   *electron_Energy_Smear_Up;
   vector<float>   *electron_Energy_Smear_Down;
   vector<float>   *electron_eta;
   vector<float>   *electron_phi;
   vector<int>     *electron_charge;
   vector<float>   *electron_gsfpt;
   vector<float>   *electron_gsfEta;
   vector<float>   *electron_gsfPhi;
   vector<int>     *electron_gsfCharge;
   vector<float>   *electron_scEta;
   vector<float>   *electron_scPhi;
   vector<float>   *electron_etaWidth;
   vector<float>   *electron_phiWidth;
   vector<float>   *electron_dEtaIn;
   vector<float>   *electron_dEtaInSeed;
   vector<float>   *electron_dPhiIn;
   vector<float>   *electron_sigmaIEtaIEta;
   vector<float>   *electron_Full5x5_SigmaIEtaIEta;
   vector<float>   *electron_e2x5OverE5x5;
   vector<float>   *electron_e1x5OverE5x5;
   vector<float>   *electron_HoverE;
   vector<float>   *electron_fbrem;
   vector<float>   *electron_eOverP;
   vector<float>   *electron_InvEminusInvP;
   vector<float>   *electron_dxyVTX;
   vector<float>   *electron_dxyerrVTX;
   vector<float>   *electron_dzVTX;
   vector<float>   *electron_dzerrVTX;
   vector<float>   *electron_3DIPVTX;
   vector<float>   *electron_3DIPerrVTX;
   vector<float>   *electron_dxy;
   vector<float>   *electron_sigdxy;
   vector<float>   *electron_dz;
   vector<float>   *electron_dxyBS;
   vector<float>   *electron_dzBS;
   vector<float>   *electron_chIso03;
   vector<float>   *electron_nhIso03;
   vector<float>   *electron_phIso03;
   vector<float>   *electron_puChIso03;
   vector<bool>    *electron_passConversionVeto;
   vector<bool>    *electron_isGsfCtfScPixChargeConsistent;
   vector<bool>    *electron_isGsfScPixChargeConsistent;
   vector<bool>    *electron_isGsfCtfChargeConsistent;
   vector<int>     *electron_mHits;
   vector<int>     *electron_ecalDriven;
   vector<float>   *electron_r9;
   vector<float>   *electron_l1et;
   vector<float>   *electron_scEnergy;
   vector<float>   *electron_scPreEnergy;
   vector<float>   *electron_scRawEnergy;
   vector<float>   *electron_scEt;
   vector<float>   *electron_E15;
   vector<float>   *electron_E25;
   vector<float>   *electron_E55;
   vector<float>   *electron_RelPFIso_dBeta;
   vector<float>   *electron_RelPFIso_Rho;
   vector<unsigned int> *electron_IDBit;
   vector<int>     *electron_IDCutBit;
   vector<float>   *electron_EnergyUnCorr;
   vector<float>   *electron_chMiniIso;
   vector<float>   *electron_nhMiniIso;
   vector<float>   *electron_phMiniIso;
   vector<float>   *electron_puChMiniIso;
   vector<float>   *electron_trackIso;
   vector<float>   *electron_dr03EcalRecHitSumEt;
   vector<float>   *electron_dr03HcalDepth1TowerSumEt;
   vector<float>   *electron_dr03HcalTowerSumEt;
   vector<float>   *electron_dr03TkSumPt;
   vector<float>   *electron_ecalPFClusterIso;
   vector<float>   *electron_hcalPFClusterIso;
   vector<ULong64_t> *electron_pathbits;
   vector<ULong64_t> *electron_filterbits;
   vector<float>   *muon_PfChargedHadronIsoR04;
   vector<float>   *muon_PfNeutralHadronIsoR04;
   vector<float>   *muon_PfGammaIsoR04;
   vector<float>   *muon_PFSumPUIsoR04;
   vector<float>   *muon_PfChargedHadronIsoR03;
   vector<float>   *muon_PfNeutralHadronIsoR03;
   vector<float>   *muon_PfGammaIsoR03;
   vector<float>   *muon_PFSumPUIsoR03;
   vector<unsigned int> *muon_TypeBit;
   vector<unsigned int> *muon_IDBit;
   vector<bool>    *muon_ishighpt;
   vector<bool>    *muon_ismedium_hip;
   vector<bool>    *muon_ismedium_nohip;
   vector<float>   *muon_dB;
   vector<float>   *muon_phi;
   vector<float>   *muon_eta;
   vector<float>   *muon_pt;
   vector<float>   *muon_mass;
   vector<float>   *muon_trkiso;
   vector<float>   *muon_hcaliso;
   vector<float>   *muon_ecaliso;
   vector<float>   *muon_trkisoR05;
   vector<float>   *muon_hcalisoR05;
   vector<float>   *muon_ecalisoR05;
   vector<int>     *muon_charge;
   vector<int>     *muon_nChambers;
   vector<int>     *muon_matchedstations;
   vector<int>     *muon_stationMask;
   vector<int>     *muon_nSegments;
   vector<float>   *muon_normchi;
   vector<int>     *muon_validhits;
   vector<int>     *muon_trackerHits;
   vector<int>     *muon_pixelHits;
   vector<int>     *muon_validmuonhits;
   vector<int>     *muon_trackerLayers;
   vector<float>   *muon_qoverp;
   vector<float>   *muon_theta;
   vector<float>   *muon_lambda;
   vector<float>   *muon_dxy;
   vector<float>   *muon_d0;
   vector<float>   *muon_dsz;
   vector<float>   *muon_dz;
   vector<float>   *muon_dxyBS;
   vector<float>   *muon_dzBS;
   vector<float>   *muon_dszBS;
   vector<float>   *muon_dxyVTX;
   vector<float>   *muon_dxyerrVTX;
   vector<float>   *muon_dzVTX;
   vector<float>   *muon_dzerrVTX;
   vector<float>   *muon_3DIPVTX;
   vector<float>   *muon_3DIPerrVTX;
   vector<float>   *muon_dszVTX;
   vector<float>   *muon_vx;
   vector<float>   *muon_vy;
   vector<float>   *muon_vz;
   vector<float>   *muon_Best_pt;
   vector<float>   *muon_Best_ptError;
   vector<float>   *muon_Best_eta;
   vector<float>   *muon_Best_phi;
   vector<float>   *muon_Inner_pt;
   vector<float>   *muon_Inner_ptError;
   vector<float>   *muon_Inner_eta;
   vector<float>   *muon_Inner_phi;
   vector<float>   *muon_Outer_pt;
   vector<float>   *muon_Outer_ptError;
   vector<float>   *muon_Outer_eta;
   vector<float>   *muon_Outer_phi;
   vector<float>   *muon_GLB_pt;
   vector<float>   *muon_GLB_ptError;
   vector<float>   *muon_GLB_eta;
   vector<float>   *muon_GLB_phi;
   vector<float>   *muon_TuneP_pt;
   vector<float>   *muon_TuneP_ptError;
   vector<float>   *muon_TuneP_eta;
   vector<float>   *muon_TuneP_phi;
   vector<float>   *muon_TuneP_charge;
   vector<float>   *muon_roch_sf;
   vector<float>   *muon_roch_sf_up;
   vector<float>   *muon_PfChargedHadronMiniIso;
   vector<float>   *muon_PfNeutralHadronMiniIso;
   vector<float>   *muon_PfGammaMiniIso;
   vector<float>   *muon_PFSumPUMiniIso;
   vector<float>   *muon_MVA;
   vector<float>   *muon_lowptMVA;
   vector<float>   *muon_softMVA;
   vector<float>   *muon_jetPtRatio;
   vector<float>   *muon_jetPtRel;
   vector<int>     *muon_simType;
   vector<int>     *muon_simExtType;
   vector<int>     *muon_simFlavour;
   vector<int>     *muon_simHeaviestMotherFlavour;
   vector<int>     *muon_simPdgId;
   vector<int>     *muon_simMotherPdgId;
   vector<float>   *muon_simMatchQuality;
   vector<ULong64_t> *muon_pathbits;
   vector<ULong64_t> *muon_filterbits;
   vector<float>   *LHE_Px;
   vector<float>   *LHE_Py;
   vector<float>   *LHE_Pz;
   vector<float>   *LHE_E;
   vector<int>     *LHE_Status;
   vector<int>     *LHE_ID;
   vector<float>   *weight_AlphaS;
   vector<float>   *weight_PDF;
   vector<float>   *weight_Scale;
   vector<float>   *weight_PSSyst;
   vector<float>   *weight_alpsfact;
   vector<float>   *weight_largeptscales;
   vector<float>   *weight_q0;
   vector<float>   *weight_sthw2;
   vector<float>   *weight_kHtt;
   vector<float>   *gen_phi;
   vector<float>   *gen_eta;
   vector<float>   *gen_pt;
   vector<float>   *gen_mass;
   vector<float>   *gen_charge;
   vector<int>     *gen_mother_index;
   vector<int>     *gen_status;
   vector<int>     *gen_PID;
   vector<bool>    *gen_isPrompt;
   vector<bool>    *gen_isPromptFinalState;
   vector<bool>    *gen_isTauDecayProduct;
   vector<bool>    *gen_isPromptTauDecayProduct;
   vector<bool>    *gen_isDirectPromptTauDecayProductFinalState;
   vector<bool>    *gen_isHardProcess;
   vector<bool>    *gen_isLastCopy;
   vector<bool>    *gen_isLastCopyBeforeFSR;
   vector<bool>    *gen_isPromptDecayed;
   vector<bool>    *gen_isDecayedLeptonHadron;
   vector<bool>    *gen_fromHardProcessBeforeFSR;
   vector<bool>    *gen_fromHardProcessDecayed;
   vector<bool>    *gen_fromHardProcessFinalState;
   vector<bool>    *gen_isMostlyLikePythia6Status3;
   Float_t         gen_weight;
   Float_t         genWeight_Q;
   Float_t         genWeight_X1;
   Float_t         genWeight_X2;
   Int_t           genWeight_id1;
   Int_t           genWeight_id2;
   Float_t         genWeight_alphaQCD;
   Float_t         genWeight_alphaQED;
   vector<float>   *photon_Energy;
   vector<float>   *photon_EnergyUnCorr;
   vector<float>   *photon_eta;
   vector<float>   *photon_phi;
   vector<float>   *photon_scEta;
   vector<float>   *photon_scPhi;
   vector<float>   *photon_HoverE;
   vector<bool>    *photon_hasPixelSeed;
   vector<float>   *photon_Full5x5_SigmaIEtaIEta;
   vector<float>   *photon_ChIso;
   vector<float>   *photon_NhIso;
   vector<float>   *photon_PhIso;
   vector<float>   *photon_ChIsoWithEA;
   vector<float>   *photon_NhIsoWithEA;
   vector<float>   *photon_PhIsoWithEA;
   vector<bool>    *photon_passMVAID_WP80;
   vector<bool>    *photon_passMVAID_WP90;
   vector<bool>    *photon_passLooseID;
   vector<bool>    *photon_passMediumID;
   vector<bool>    *photon_passTightID;
   Int_t           nPileUp;
   Float_t         pfMET_pt;
   Float_t         pfMET_phi;
   Float_t         pfMET_SumEt;
   Float_t         pfMET_Type1_pt;
   Float_t         pfMET_Type1_phi;
   Float_t         pfMET_Type1_SumEt;
   Float_t         pfMET_Type1_PhiCor_pt;
   Float_t         pfMET_Type1_PhiCor_phi;
   Float_t         pfMET_Type1_PhiCor_SumEt;
   vector<float>   *pfMET_Type1_pt_shifts;
   vector<float>   *pfMET_Type1_phi_shifts;
   vector<float>   *pfMET_Type1_SumEt_shifts;
   vector<float>   *pfMET_Type1_PhiCor_pt_shifts;
   vector<float>   *pfMET_Type1_PhiCor_phi_shifts;
   vector<float>   *pfMET_Type1_PhiCor_SumEt_shifts;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_SumEt;
   Float_t         PuppiMET_Type1_pt;
   Float_t         PuppiMET_Type1_phi;
   Float_t         PuppiMET_Type1_SumEt;
   Float_t         PuppiMET_Type1_PhiCor_pt;
   Float_t         PuppiMET_Type1_PhiCor_phi;
   Float_t         PuppiMET_Type1_PhiCor_SumEt;
   vector<float>   *PuppiMET_Type1_pt_shifts;
   vector<float>   *PuppiMET_Type1_phi_shifts;
   vector<float>   *PuppiMET_Type1_SumEt_shifts;
   vector<float>   *tau_phi;
   vector<float>   *tau_eta;
   vector<float>   *tau_pt;
   vector<float>   *tau_mass;
   vector<float>   *tau_dz;
   vector<float>   *tau_dxy;
   vector<int>     *tau_decaymode;
   vector<int>     *tau_charge;
   vector<unsigned int> *tau_IDBit;
   vector<bool>    *tau_idDecayModeNewDMs;

   // List of branches
   TBranch        *b_IsData;   //!
   TBranch        *b_nTotal;   //!
   TBranch        *b_runNum;   //!
   TBranch        *b_evtNum;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_RhoNC;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_L1PrefireReweight_Central;   //!
   TBranch        *b_L1PrefireReweight_Up;   //!
   TBranch        *b_L1PrefireReweight_Down;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_PVtrackSize;   //!
   TBranch        *b_PVchi2;   //!
   TBranch        *b_PVndof;   //!
   TBranch        *b_PVnormalizedChi2;   //!
   TBranch        *b_PVx;   //!
   TBranch        *b_PVy;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_HLT_TriggerName;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_charge;   //!
   TBranch        *b_jet_area;   //!
   TBranch        *b_jet_partonFlavour;   //!
   TBranch        *b_jet_hadronFlavour;   //!
   TBranch        *b_jet_GenHFHadronMatcher_flavour;   //!
   TBranch        *b_jet_GenHFHadronMatcher_origin;   //!
   TBranch        *b_jet_DeepCSV;   //!
   TBranch        *b_jet_DeepCSV_CvsL;   //!
   TBranch        *b_jet_DeepCSV_CvsB;   //!
   TBranch        *b_jet_DeepJet;   //!
   TBranch        *b_jet_DeepJet_CvsL;   //!
   TBranch        *b_jet_DeepJet_CvsB;   //!
   TBranch        *b_jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_jet_neutralHadronEnergyFraction;   //!
   TBranch        *b_jet_neutralEmEnergyFraction;   //!
   TBranch        *b_jet_chargedEmEnergyFraction;   //!
   TBranch        *b_jet_muonEnergyFraction;   //!
   TBranch        *b_jet_chargedMultiplicity;   //!
   TBranch        *b_jet_neutralMultiplicity;   //!
   TBranch        *b_jet_tightJetID;   //!
   TBranch        *b_jet_tightLepVetoJetID;   //!
   TBranch        *b_jet_partonPdgId;   //!
   TBranch        *b_jet_vtxNtracks;   //!
   TBranch        *b_jet_m;   //!
   TBranch        *b_jet_energy;   //!
   TBranch        *b_jet_PileupJetId;   //!
   TBranch        *b_jet_shiftedEnUp;   //!
   TBranch        *b_jet_shiftedEnDown;   //!
   TBranch        *b_jet_smearedRes;   //!
   TBranch        *b_jet_smearedResUp;   //!
   TBranch        *b_jet_smearedResDown;   //!
   TBranch        *b_jet_bJetNN_corr;   //!
   TBranch        *b_jet_bJetNN_res;   //!
   TBranch        *b_jet_cJetNN_corr;   //!
   TBranch        *b_jet_cJetNN_res;   //!
   TBranch        *b_jet_JECL1FastJet;   //!
   TBranch        *b_jet_JECFull;   //!
   TBranch        *b_fatjet_pt;   //!
   TBranch        *b_fatjet_eta;   //!
   TBranch        *b_fatjet_phi;   //!
   TBranch        *b_fatjet_charge;   //!
   TBranch        *b_fatjet_area;   //!
   TBranch        *b_fatjet_partonFlavour;   //!
   TBranch        *b_fatjet_hadronFlavour;   //!
   TBranch        *b_fatjet_DeepCSV;   //!
   TBranch        *b_fatjet_DeepCSV_CvsL;   //!
   TBranch        *b_fatjet_DeepCSV_CvsB;   //!
   TBranch        *b_fatjet_particleNet_TvsQCD;   //!
   TBranch        *b_fatjet_particleNet_WvsQCD;   //!
   TBranch        *b_fatjet_particleNet_ZvsQCD;   //!
   TBranch        *b_fatjet_particleNet_HbbvsQCD;   //!
   TBranch        *b_fatjet_particleNet_HccvsQCD;   //!
   TBranch        *b_fatjet_particleNet_H4qvsQCD;   //!
   TBranch        *b_fatjet_particleNet_QCD;   //!
   TBranch        *b_fatjet_particleNetMD_Xbb;   //!
   TBranch        *b_fatjet_particleNetMD_Xcc;   //!
   TBranch        *b_fatjet_particleNetMD_Xqq;   //!
   TBranch        *b_fatjet_particleNetMD_QCD;   //!
   TBranch        *b_fatjet_tightJetID;   //!
   TBranch        *b_fatjet_tightLepVetoJetID;   //!
   TBranch        *b_fatjet_partonPdgId;   //!
   TBranch        *b_fatjet_vtxNtracks;   //!
   TBranch        *b_fatjet_m;   //!
   TBranch        *b_fatjet_energy;   //!
   TBranch        *b_fatjet_puppi_tau1;   //!
   TBranch        *b_fatjet_puppi_tau2;   //!
   TBranch        *b_fatjet_puppi_tau3;   //!
   TBranch        *b_fatjet_puppi_tau4;   //!
   TBranch        *b_fatjet_softdropmass;   //!
   TBranch        *b_fatjet_chargedHadronEnergyFraction;   //!
   TBranch        *b_fatjet_neutralHadronEnergyFraction;   //!
   TBranch        *b_fatjet_neutralEmEnergyFraction;   //!
   TBranch        *b_fatjet_chargedEmEnergyFraction;   //!
   TBranch        *b_fatjet_muonEnergyFraction;   //!
   TBranch        *b_fatjet_chargedMultiplicity;   //!
   TBranch        *b_fatjet_neutralMultiplicity;   //!
   TBranch        *b_fatjet_shiftedEnUp;   //!
   TBranch        *b_fatjet_shiftedEnDown;   //!
   TBranch        *b_fatjet_smearedRes;   //!
   TBranch        *b_fatjet_smearedResUp;   //!
   TBranch        *b_fatjet_smearedResDown;   //!
   TBranch        *b_fatjet_LSF;   //!
   TBranch        *b_fatjet_LSFlep_PID;   //!
   TBranch        *b_fatjet_LSFlep_Pt;   //!
   TBranch        *b_fatjet_LSFlep_Eta;   //!
   TBranch        *b_fatjet_LSFlep_Phi;   //!
   TBranch        *b_electron_MVAIso;   //!
   TBranch        *b_electron_MVANoIso;   //!
   TBranch        *b_electron_Energy;   //!
   TBranch        *b_electron_Energy_Scale_Up;   //!
   TBranch        *b_electron_Energy_Scale_Down;   //!
   TBranch        *b_electron_Energy_Smear_Up;   //!
   TBranch        *b_electron_Energy_Smear_Down;   //!
   TBranch        *b_electron_eta;   //!
   TBranch        *b_electron_phi;   //!
   TBranch        *b_electron_charge;   //!
   TBranch        *b_electron_gsfpt;   //!
   TBranch        *b_electron_gsfEta;   //!
   TBranch        *b_electron_gsfPhi;   //!
   TBranch        *b_electron_gsfCharge;   //!
   TBranch        *b_electron_scEta;   //!
   TBranch        *b_electron_scPhi;   //!
   TBranch        *b_electron_etaWidth;   //!
   TBranch        *b_electron_phiWidth;   //!
   TBranch        *b_electron_dEtaIn;   //!
   TBranch        *b_electron_dEtaInSeed;   //!
   TBranch        *b_electron_dPhiIn;   //!
   TBranch        *b_electron_sigmaIEtaIEta;   //!
   TBranch        *b_electron_Full5x5_SigmaIEtaIEta;   //!
   TBranch        *b_electron_e2x5OverE5x5;   //!
   TBranch        *b_electron_e1x5OverE5x5;   //!
   TBranch        *b_electron_HoverE;   //!
   TBranch        *b_electron_fbrem;   //!
   TBranch        *b_electron_eOverP;   //!
   TBranch        *b_electron_InvEminusInvP;   //!
   TBranch        *b_electron_dxyVTX;   //!
   TBranch        *b_electron_dxyerrVTX;   //!
   TBranch        *b_electron_dzVTX;   //!
   TBranch        *b_electron_dzerrVTX;   //!
   TBranch        *b_electron_3DIPVTX;   //!
   TBranch        *b_electron_3DIPerrVTX;   //!
   TBranch        *b_electron_dxy;   //!
   TBranch        *b_electron_sigdxy;   //!
   TBranch        *b_electron_dz;   //!
   TBranch        *b_electron_dxyBS;   //!
   TBranch        *b_electron_dzBS;   //!
   TBranch        *b_electron_chIso03;   //!
   TBranch        *b_electron_nhIso03;   //!
   TBranch        *b_electron_phIso03;   //!
   TBranch        *b_electron_puChIso03;   //!
   TBranch        *b_electron_passConversionVeto;   //!
   TBranch        *b_electron_isGsfCtfScPixChargeConsistent;   //!
   TBranch        *b_electron_isGsfScPixChargeConsistent;   //!
   TBranch        *b_electron_isGsfCtfChargeConsistent;   //!
   TBranch        *b_electron_mHits;   //!
   TBranch        *b_electron_ecalDriven;   //!
   TBranch        *b_electron_r9;   //!
   TBranch        *b_electron_l1et;   //!
   TBranch        *b_electron_scEnergy;   //!
   TBranch        *b_electron_scPreEnergy;   //!
   TBranch        *b_electron_scRawEnergy;   //!
   TBranch        *b_electron_scEt;   //!
   TBranch        *b_electron_E15;   //!
   TBranch        *b_electron_E25;   //!
   TBranch        *b_electron_E55;   //!
   TBranch        *b_electron_RelPFIso_dBeta;   //!
   TBranch        *b_electron_RelPFIso_Rho;   //!
   TBranch        *b_electron_IDBit;   //!
   TBranch        *b_electron_IDCutBit;   //!
   TBranch        *b_electron_EnergyUnCorr;   //!
   TBranch        *b_electron_chMiniIso;   //!
   TBranch        *b_electron_nhMiniIso;   //!
   TBranch        *b_electron_phMiniIso;   //!
   TBranch        *b_electron_puChMiniIso;   //!
   TBranch        *b_electron_trackIso;   //!
   TBranch        *b_electron_dr03EcalRecHitSumEt;   //!
   TBranch        *b_electron_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_electron_dr03HcalTowerSumEt;   //!
   TBranch        *b_electron_dr03TkSumPt;   //!
   TBranch        *b_electron_ecalPFClusterIso;   //!
   TBranch        *b_electron_hcalPFClusterIso;   //!
   TBranch        *b_electron_pathbits;   //!
   TBranch        *b_electron_filterbits;   //!
   TBranch        *b_muon_PfChargedHadronIsoR04;   //!
   TBranch        *b_muon_PfNeutralHadronIsoR04;   //!
   TBranch        *b_muon_PfGammaIsoR04;   //!
   TBranch        *b_muon_PFSumPUIsoR04;   //!
   TBranch        *b_muon_PfChargedHadronIsoR03;   //!
   TBranch        *b_muon_PfNeutralHadronIsoR03;   //!
   TBranch        *b_muon_PfGammaIsoR03;   //!
   TBranch        *b_muon_PFSumPUIsoR03;   //!
   TBranch        *b_muon_TypeBit;   //!
   TBranch        *b_muon_IDBit;   //!
   TBranch        *b_muon_ishighpt;   //!
   TBranch        *b_muon_ismedium_hip;   //!
   TBranch        *b_muon_ismedium_nohip;   //!
   TBranch        *b_muon_dB;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_mass;   //!
   TBranch        *b_muon_trkiso;   //!
   TBranch        *b_muon_hcaliso;   //!
   TBranch        *b_muon_ecaliso;   //!
   TBranch        *b_muon_trkisoR05;   //!
   TBranch        *b_muon_hcalisoR05;   //!
   TBranch        *b_muon_ecalisoR05;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_nChambers;   //!
   TBranch        *b_muon_matchedstations;   //!
   TBranch        *b_muon_stationMask;   //!
   TBranch        *b_muon_nSegments;   //!
   TBranch        *b_muon_normchi;   //!
   TBranch        *b_muon_validhits;   //!
   TBranch        *b_muon_trackerHits;   //!
   TBranch        *b_muon_pixelHits;   //!
   TBranch        *b_muon_validmuonhits;   //!
   TBranch        *b_muon_trackerLayers;   //!
   TBranch        *b_muon_qoverp;   //!
   TBranch        *b_muon_theta;   //!
   TBranch        *b_muon_lambda;   //!
   TBranch        *b_muon_dxy;   //!
   TBranch        *b_muon_d0;   //!
   TBranch        *b_muon_dsz;   //!
   TBranch        *b_muon_dz;   //!
   TBranch        *b_muon_dxyBS;   //!
   TBranch        *b_muon_dzBS;   //!
   TBranch        *b_muon_dszBS;   //!
   TBranch        *b_muon_dxyVTX;   //!
   TBranch        *b_muon_dxyerrVTX;   //!
   TBranch        *b_muon_dzVTX;   //!
   TBranch        *b_muon_dzerrVTX;   //!
   TBranch        *b_muon_3DIPVTX;   //!
   TBranch        *b_muon_3DIPerrVTX;   //!
   TBranch        *b_muon_dszVTX;   //!
   TBranch        *b_muon_vx;   //!
   TBranch        *b_muon_vy;   //!
   TBranch        *b_muon_vz;   //!
   TBranch        *b_muon_Best_pt;   //!
   TBranch        *b_muon_Best_ptError;   //!
   TBranch        *b_muon_Best_eta;   //!
   TBranch        *b_muon_Best_phi;   //!
   TBranch        *b_muon_Inner_pt;   //!
   TBranch        *b_muon_Inner_ptError;   //!
   TBranch        *b_muon_Inner_eta;   //!
   TBranch        *b_muon_Inner_phi;   //!
   TBranch        *b_muon_Outer_pt;   //!
   TBranch        *b_muon_Outer_ptError;   //!
   TBranch        *b_muon_Outer_eta;   //!
   TBranch        *b_muon_Outer_phi;   //!
   TBranch        *b_muon_GLB_pt;   //!
   TBranch        *b_muon_GLB_ptError;   //!
   TBranch        *b_muon_GLB_eta;   //!
   TBranch        *b_muon_GLB_phi;   //!
   TBranch        *b_muon_TuneP_pt;   //!
   TBranch        *b_muon_TuneP_ptError;   //!
   TBranch        *b_muon_TuneP_eta;   //!
   TBranch        *b_muon_TuneP_phi;   //!
   TBranch        *b_muon_TuneP_charge;   //!
   TBranch        *b_muon_roch_sf;   //!
   TBranch        *b_muon_roch_sf_up;   //!
   TBranch        *b_muon_PfChargedHadronMiniIso;   //!
   TBranch        *b_muon_PfNeutralHadronMiniIso;   //!
   TBranch        *b_muon_PfGammaMiniIso;   //!
   TBranch        *b_muon_PFSumPUMiniIso;   //!
   TBranch        *b_muon_MVA;   //!
   TBranch        *b_muon_lowptMVA;   //!
   TBranch        *b_muon_softMVA;   //!
   TBranch        *b_muon_jetPtRatio;   //!
   TBranch        *b_muon_jetPtRel;   //!
   TBranch        *b_muon_simType;   //!
   TBranch        *b_muon_simExtType;   //!
   TBranch        *b_muon_simFlavour;   //!
   TBranch        *b_muon_simHeaviestMotherFlavour;   //!
   TBranch        *b_muon_simPdgId;   //!
   TBranch        *b_muon_simMotherPdgId;   //!
   TBranch        *b_muon_simMatchQuality;   //!
   TBranch        *b_muon_pathbits;   //!
   TBranch        *b_muon_filterbits;   //!
   TBranch        *b_LHE_Px;   //!
   TBranch        *b_LHE_Py;   //!
   TBranch        *b_LHE_Pz;   //!
   TBranch        *b_LHE_E;   //!
   TBranch        *b_LHE_Status;   //!
   TBranch        *b_LHE_ID;   //!
   TBranch        *b_weight_AlphaS;   //!
   TBranch        *b_weight_PDF;   //!
   TBranch        *b_weight_Scale;   //!
   TBranch        *b_weight_PSSyst;   //!
   TBranch        *b_weight_alpsfact;   //!
   TBranch        *b_weight_largeptscales;   //!
   TBranch        *b_weight_q0;   //!
   TBranch        *b_weight_sthw2;   //!
   TBranch        *b_weight_kHtt;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_mass;   //!
   TBranch        *b_gen_charge;   //!
   TBranch        *b_gen_mother_index;   //!
   TBranch        *b_gen_status;   //!
   TBranch        *b_gen_PID;   //!
   TBranch        *b_gen_isPrompt;   //!
   TBranch        *b_gen_isPromptFinalState;   //!
   TBranch        *b_gen_isTauDecayProduct;   //!
   TBranch        *b_gen_isPromptTauDecayProduct;   //!
   TBranch        *b_gen_isDirectPromptTauDecayProductFinalState;   //!
   TBranch        *b_gen_isHardProcess;   //!
   TBranch        *b_gen_isLastCopy;   //!
   TBranch        *b_gen_isLastCopyBeforeFSR;   //!
   TBranch        *b_gen_isPromptDecayed;   //!
   TBranch        *b_gen_isDecayedLeptonHadron;   //!
   TBranch        *b_gen_fromHardProcessBeforeFSR;   //!
   TBranch        *b_gen_fromHardProcessDecayed;   //!
   TBranch        *b_gen_fromHardProcessFinalState;   //!
   TBranch        *b_gen_isMostlyLikePythia6Status3;   //!
   TBranch        *b_gen_weight;   //!
   TBranch        *b_genWeight_Q;   //!
   TBranch        *b_genWeight_X1;   //!
   TBranch        *b_genWeight_X2;   //!
   TBranch        *b_genWeight_id1;   //!
   TBranch        *b_genWeight_id2;   //!
   TBranch        *b_genWeight_alphaQCD;   //!
   TBranch        *b_genWeight_alphaQED;   //!
   TBranch        *b_photon_Energy;   //!
   TBranch        *b_photon_EnergyUnCorr;   //!
   TBranch        *b_photon_eta;   //!
   TBranch        *b_photon_phi;   //!
   TBranch        *b_photon_scEta;   //!
   TBranch        *b_photon_scPhi;   //!
   TBranch        *b_photon_HoverE;   //!
   TBranch        *b_photon_hasPixelSeed;   //!
   TBranch        *b_photon_Full5x5_SigmaIEtaIEta;   //!
   TBranch        *b_photon_ChIso;   //!
   TBranch        *b_photon_NhIso;   //!
   TBranch        *b_photon_PhIso;   //!
   TBranch        *b_photon_ChIsoWithEA;   //!
   TBranch        *b_photon_NhIsoWithEA;   //!
   TBranch        *b_photon_PhIsoWithEA;   //!
   TBranch        *b_photon_passMVAID_WP80;   //!
   TBranch        *b_photon_passMVAID_WP90;   //!
   TBranch        *b_photon_passLooseID;   //!
   TBranch        *b_photon_passMediumID;   //!
   TBranch        *b_photon_passTightID;   //!
   TBranch        *b_nPileUp;   //!
   TBranch        *b_pfMET_pt;   //!
   TBranch        *b_pfMET_phi;   //!
   TBranch        *b_pfMET_SumEt;   //!
   TBranch        *b_pfMET_Type1_pt;   //!
   TBranch        *b_pfMET_Type1_phi;   //!
   TBranch        *b_pfMET_Type1_SumEt;   //!
   TBranch        *b_pfMET_Type1_PhiCor_pt;   //!
   TBranch        *b_pfMET_Type1_PhiCor_phi;   //!
   TBranch        *b_pfMET_Type1_PhiCor_SumEt;   //!
   TBranch        *b_pfMET_Type1_pt_shifts;   //!
   TBranch        *b_pfMET_Type1_phi_shifts;   //!
   TBranch        *b_pfMET_Type1_SumEt_shifts;   //!
   TBranch        *b_pfMET_Type1_PhiCor_pt_shifts;   //!
   TBranch        *b_pfMET_Type1_PhiCor_phi_shifts;   //!
   TBranch        *b_pfMET_Type1_PhiCor_SumEt_shifts;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_SumEt;   //!
   TBranch        *b_PuppiMET_Type1_pt;   //!
   TBranch        *b_PuppiMET_Type1_phi;   //!
   TBranch        *b_PuppiMET_Type1_SumEt;   //!
   TBranch        *b_PuppiMET_Type1_PhiCor_pt;   //!
   TBranch        *b_PuppiMET_Type1_PhiCor_phi;   //!
   TBranch        *b_PuppiMET_Type1_PhiCor_SumEt;   //!
   TBranch        *b_PuppiMET_Type1_pt_shifts;   //!
   TBranch        *b_PuppiMET_Type1_phi_shifts;   //!
   TBranch        *b_PuppiMET_Type1_SumEt_shifts;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_mass;   //!
   TBranch        *b_tau_dz;   //!
   TBranch        *b_tau_dxy;   //!
   TBranch        *b_tau_decaymode;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_tau_IDBit;   //!
   TBranch        *b_tau_idDecayModeNewDMs;   //!

};

#endif

