#include "FatJet.h"

ClassImp(FatJet)

FatJet::FatJet() : Particle() {
  j_area=-999.;
  j_partonFlavour=-999;
  j_hadronFlavour=-999;
  j_DeepCSV=-999.;
  j_DeepCSV_CvsL=-999.;
  j_DeepCSV_CvsB=-999.;
  j_particleNet_TvsQCD=-999.;
  j_particleNet_WvsQCD=-999.;
  j_particleNet_ZvsQCD=-999.;
  j_particleNet_HbbvsQCD=-999.;
  j_particleNet_HccvsQCD=-999.;
  j_particleNet_H4qvsQCD=-999.;
  j_particleNet_QCD=-999.;
  j_particleNetMD_Xbb=-999.;
  j_particleNetMD_Xcc=-999.;
  j_particleNetMD_Xqq=-999.;
  j_particleNetMD_QCD=-999.;
  j_chargedHadronEnergyFraction=-999.;
  j_neutralHadronEnergyFraction=-999.;
  j_neutralEmEnergyFraction=-999.;
  j_chargedEmEnergyFraction=-999.;
  j_muonEnergyFraction=-999.;
  j_chargedMultiplicity=-999;
  j_neutralMultiplicity=-999;
  j_lsf=-999.;
  j_lsf_pid=-999;
  j_En_up=1.;
  j_En_down=1.;;
  j_Res_up = 1.;
  j_Res_down = 1.;
  j_tightJetID=false;
  j_tightLepVetoJetID=false;
  j_puppi_tau1 = -999.;
  j_puppi_tau2 = -999.;
  j_puppi_tau3 = -999.;
  j_puppi_tau4 = -999.;
  j_SDMass = -999.;
}

FatJet::~FatJet(){

}

void FatJet::SetArea(double area){
  j_area = area;
}
void FatJet::SetGenFlavours(double pf, double hf){
  j_partonFlavour = pf;
  j_hadronFlavour = hf;
}
void FatJet::SetTaggerResults(std::vector<double> ds){
  j_DeepCSV           = ds.at(0);
  j_DeepCSV_CvsL      = ds.at(1);
  j_DeepCSV_CvsB      = ds.at(2);
  j_particleNet_TvsQCD   = ds.at(3);
  j_particleNet_WvsQCD   = ds.at(4);
  j_particleNet_ZvsQCD   = ds.at(5);
  j_particleNet_HbbvsQCD = ds.at(6);
  j_particleNet_HccvsQCD = ds.at(7);
  j_particleNet_H4qvsQCD = ds.at(8);
  j_particleNet_QCD      = ds.at(9);
  j_particleNetMD_Xbb    = ds.at(10);
  j_particleNetMD_Xcc    = ds.at(11);
  j_particleNetMD_Xqq    = ds.at(12);
  j_particleNetMD_QCD    = ds.at(13);
}
void FatJet::SetEnergyFractions(double cH, double nH, double nEM, double cEM, double muE){
  j_chargedHadronEnergyFraction = cH;
  j_neutralHadronEnergyFraction = nH;
  j_neutralEmEnergyFraction = nEM;
  j_chargedEmEnergyFraction = cEM;
  j_muonEnergyFraction = muE;
}
void FatJet::SetMultiplicities(double cM, double nM){
  j_chargedMultiplicity = cM;
  j_neutralMultiplicity = nM;
}

void FatJet::SetLSF(double lsf, int lsf_PID){
  j_lsf = lsf;
  j_lsf_pid = lsf_PID;
}
double FatJet::LSF() const{
  return j_lsf;
}
double FatJet::LSF_PID() const{
  return j_lsf_pid;
}

void FatJet::SetEnShift(double en_up, double en_down){
  j_En_up = en_up;
  j_En_down = en_down;
}

void FatJet::SetResShift(double res_up, double res_down){
  j_Res_up = res_up;
  j_Res_down = res_down;
}

void FatJet::SetTightJetID(double b){
  j_tightJetID = b;
}
void FatJet::SetTightLepVetoJetID(double b){
  j_tightLepVetoJetID = b;
}

bool FatJet::PassID(TString ID) const {

  if(ID=="tight") return Pass_tightJetID();
  if(ID=="tightLepVeto") return Pass_tightLepVetoJetID();

  cout << "[FatJet::PassID] No id : " << ID << endl;
  exit(ENODATA);

  return false;

}

double FatJet::GetTaggerResult(JetTagging::Tagger tg) const{
  if(tg==JetTagging::DeepCSV) return j_DeepCSV;
  else if(tg==JetTagging::DeepCSV_CvsL) return j_DeepCSV_CvsL;
  else if(tg==JetTagging::DeepCSV_CvsB) return j_DeepCSV_CvsB;
  else if(tg==JetTagging::particleNet_TvsQCD) return j_particleNet_TvsQCD;
  else if(tg==JetTagging::particleNet_WvsQCD) return j_particleNet_WvsQCD;
  else if(tg==JetTagging::particleNet_ZvsQCD) return j_particleNet_ZvsQCD;
  else if(tg==JetTagging::particleNet_HbbvsQCD) return j_particleNet_HbbvsQCD;
  else if(tg==JetTagging::particleNet_HccvsQCD) return j_particleNet_HccvsQCD;
  else if(tg==JetTagging::particleNet_H4qvsQCD) return j_particleNet_H4qvsQCD;
  else if(tg==JetTagging::particleNet_QCD) return j_particleNet_QCD;
  else if(tg==JetTagging::particleNetMD_Xbb) return j_particleNetMD_Xbb;
  else if(tg==JetTagging::particleNetMD_Xcc) return j_particleNetMD_Xcc;
  else if(tg==JetTagging::particleNetMD_Xqq) return j_particleNetMD_Xqq;
  else if(tg==JetTagging::particleNetMD_QCD) return j_particleNetMD_QCD;
  else{
    cout << "[FatJet::GetTaggerResult] ERROR; Wrong tagger" << endl;
    return -999;
  }
}

void FatJet::SetPuppiTaus(double t1, double t2, double t3, double t4){
  j_puppi_tau1 = t1;
  j_puppi_tau2 = t2;
  j_puppi_tau3 = t3;
  j_puppi_tau4 = t4;
}

void FatJet::SetSDMass(double m){
  j_SDMass = m;
}

