#include "Event.h"

ClassImp(Event)

Event::Event(){
  j_nPV = -999;
  j_DataYear = -999;
  j_DataEra = "";
}

Event::~Event(){

}

void Event::SetnPV(double n){
  j_nPV = n;
}

void Event::SetTrigger(vector<string> HLT_TriggerName){
  j_HLT_TriggerName = HLT_TriggerName;
}

bool Event::PassTrigger(TString trig){
  std::vector<TString> tmp_vec;
  tmp_vec.push_back(trig);
  return PassTrigger(tmp_vec);
}

bool Event::PassTrigger(std::vector<TString> trigs){

  for(unsigned int i=0; i<trigs.size(); i++){
    TString this_check_trig = trigs.at(i);

    //cout << this_check_trig << endl;
    for(unsigned int j=0; j<j_HLT_TriggerName.size(); j++){

      //cout << j_HLT_TriggerName.at(j) << endl;
      if( j_HLT_TriggerName.at(j).find(this_check_trig)!=std::string::npos ){
        return true;
      }
    }

  }
  return false;
}

double Event::GetTriggerLumi(TString trig){
  if(j_DataEra=="2016preVFP"){
    if(trig=="Full") return 19517.523849863;
    else if(trig=="HLT_Mu3_PFJet40_v") return 3.911703648;
    else if(trig=="HLT_Mu8_v") return 3.319291722;
    else if(trig=="HLT_Mu17_v") return 259.029439505;
    else if(trig=="HLT_Mu20_v") return 109.124528366;
    else if(trig=="HLT_Mu27_v") return 149.654826484;
    else if(trig=="HLT_Mu50_v") return 19517.523849863;
    else if(trig=="HLT_Mu55_v") return 19517.523849863;
    else if(trig=="HLT_DoublePhoton60_v") return 19517.523849863;
    else if(trig=="HLT_DoublePhoton85_v") return 19517.523849863;
    else if(trig=="HLT_Photon50_v") return 18.734556169;
    else if(trig=="HLT_Photon75_v") return 83.007650452;
    else if(trig=="HLT_Photon90_v") return 163.439249664;
    else if(trig=="HLT_Photon120_v") return 338.327429849;
    else if(trig=="HLT_Photon175_v") return 19517.523849863;
    else if(trig=="HLT_IsoMu24_v") return 19517.523849863;
    else if(trig=="HLT_IsoMu27_v") return 19517.523849863;
    else if(trig=="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 4.127188014;
    else if(trig=="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 11.033890219;
    else if(trig=="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 52.790026620;
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") return 19517.523849863;
  }else if(j_DataEra=="2016postVFP"){
    if(trig=="Full") return 16812.151722482;
    else if(trig=="HLT_Mu3_PFJet40_v") return 3.580526616;
    else if(trig=="HLT_Mu8_v") return 0.659286323;
    else if(trig=="HLT_Mu17_v") return 26.872793405;
    else if(trig=="HLT_Mu20_v") return 32.291949992;
    else if(trig=="HLT_Mu27_v") return 103.565936026;
    else if(trig=="HLT_Mu50_v") return 16812.151722482;
    else if(trig=="HLT_Mu55_v") return 16812.151722482;
    else if(trig=="HLT_DoublePhoton60_v") return 16812.151722482;
    else if(trig=="HLT_DoublePhoton85_v") return 16812.151722482;
    else if(trig=="HLT_Photon50_v") return 12.643834933;
    else if(trig=="HLT_Photon75_v") return 53.053657641;
    else if(trig=="HLT_Photon90_v") return 103.565936026;
    else if(trig=="HLT_Photon120_v") return 204.590057915;
    else if(trig=="HLT_Photon175_v") return 16812.151722482;
    else if(trig=="HLT_IsoMu24_v") return 16812.151722482;
    else if(trig=="HLT_IsoMu27_v") return 16812.151722482;
    else if(trig=="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 2.938006053;
    else if(trig=="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 3.980060594;
    else if(trig=="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 10.699483420;
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") return 16812.151722482;
  }else if(j_DataEra=="2017"){
    if(trig=="Full") return 41477.877400009;
    else if(trig=="HLT_Mu3_PFJet40_v") return 4.607782551;
    else if(trig=="HLT_Mu8_v"){
      return 2.602078443;
      //FIXME what is norm?? (JS's code)
      //==== 180829 Norm Checked : sf = 1.41036
      //==== 181206 Norm Checked : sf = 1.33461
      //return 2.605;;
      //return 2.605*1.33461;
    }
    else if(trig=="HLT_Mu17_v") return 70.005238765;
    else if(trig=="HLT_Mu19_v") return 46.204786104;
    else if(trig=="HLT_Mu20_v") return 575.585018236;
    else if(trig=="HLT_Mu27_v") return 184.843129035;
    else if(trig=="HLT_Mu50_v") return 41477.877400009;
    else if(trig=="HLT_Mu55_v") return 27121.940558082;
    else if(trig=="HLT_DoublePhoton70_v") return 41477.877400009;
    else if(trig=="HLT_DoublePhoton85_v") return 41477.877400009;
    else if(trig=="HLT_Photon25_v") return 1.339434797;
    else if(trig=="HLT_Photon33_v") return 14.319591410;
    else if(trig=="HLT_Photon50_v") return 22.182993392;
    else if(trig=="HLT_Photon75_v") return 101.406554985;
    else if(trig=="HLT_Photon90_v") return 122.130341981;
    else if(trig=="HLT_Photon120_v") return 490.132351457;
    else if(trig=="HLT_Photon150_v") return 689.862207927;
    else if(trig=="HLT_Photon175_v") return 1245.403311711;
    else if(trig=="HLT_Photon200_v") return 41477.877400009;
    else if(trig=="HLT_IsoMu24_v") return 37997.004705584;
    else if(trig=="HLT_IsoMu27_v") return 41477.877400009;
    else if(trig=="HLT_IsoMu30_v") return 27121.940558082;
    else if(trig=="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 3.970041451;
    else if(trig=="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 27.683553584;
    else if(trig=="HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 0.000600725;
    else if(trig=="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 43.453223350;
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") return 7455.394595245;
    else if(trig=="HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v") return 2613.476499000;
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v") return 41477.877400009;
    else if(trig=="HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v") return 27121.940558082;
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") return 36674.511073834;
    else if(trig=="HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v") return 27121.940558082;
  }else if(j_DataEra=="2018"){
    if(trig=="Full") return 59827.879502925;
    else if(trig=="HLT_Mu3_PFJet40_v") return 2.704239929;
    else if(trig=="HLT_Mu8_v") return 8.566603827;
    else if(trig=="HLT_Mu17_v") return 45.852032056;
    else if(trig=="HLT_Mu19_v") return 45.852032056;
    else if(trig=="HLT_Mu20_v") return 55.360539057;
    else if(trig=="HLT_Mu27_v") return 125.967123844;
    else if(trig=="HLT_Mu50_v") return 59827.879502925;
    else if(trig=="HLT_Mu55_v") return 59827.879502925;
    else if(trig=="HLT_DoublePhoton70_v") return 59827.879502925;
    else if(trig=="HLT_DoublePhoton85_v") return 59827.879502925;
    else if(trig=="HLT_Photon33_v") return 2.349799925;
    else if(trig=="HLT_Photon50_v") return 12.960071823;
    else if(trig=="HLT_Photon75_v") return 117.946263489;
    else if(trig=="HLT_Photon90_v") return 117.622232475;
    else if(trig=="HLT_Photon120_v") return 470.863308303;
    else if(trig=="HLT_Photon150_v") return 478.335294777;
    else if(trig=="HLT_Photon175_v") return 952.895805882;
    else if(trig=="HLT_Photon200_v") return 59827.879502925;
    else if(trig=="HLT_IsoMu24_v") return 59819.714470919;
    else if(trig=="HLT_IsoMu27_v") return 59827.879502925;
    else if(trig=="HLT_IsoMu30_v") return 59827.879502925;
    else if(trig=="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 6.424543223;
    else if(trig=="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 38.917235485;
    else if(trig=="HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 0.001834965;
    else if(trig=="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") return 38.973910595;
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") return 3008.003592622;
    else if(trig=="HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v") return 3135.331113947;
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v") return 59827.879502925;
    else if(trig=="HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v") return 59827.879502925;
    else if(trig=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") return 59827.879502925;
    else if(trig=="HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v") return 59827.879502925;
  }

  cout << "[Event::GetTriggerLumi] Wrong trigger (DataEra=" << j_DataEra << ", trig=" << trig << ")" << endl;
  exit(ENODATA);
  return -1.;

}

bool Event::IsPDForTrigger(TString trig, TString PD){

  if(j_DataYear==2017){
    if(trig=="HLT_Mu3_PFJet40_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Mu8_v"){
      if(PD=="DoubleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Mu17_v"){
      if(PD=="DoubleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Mu20_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Mu27_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Mu50_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_IsoMu24_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_IsoMu27_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_IsoMu30_v"){
      if(PD=="SingleMuon") return true;
      else return false;
    }
    else if(trig=="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      if(PD=="SingleElectron") return true;
      else return false;
    }
    else if(trig=="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      if(PD=="SingleElectron") return true;
      else return false;
    }
    else if(trig=="HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      if(PD=="SingleElectron") return true;
      else return false;
    }
    else if(trig=="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"){
      if(PD=="SingleElectron") return true;
      else return false;
    }
    else if(trig=="HLT_Photon25_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon33_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon50_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon75_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon90_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon120_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon150_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon175_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else if(trig=="HLT_Photon200_v"){
      if(PD=="SinglePhoton") return true;
      else return false;
    }
    else{
      cout << "[Event::IsPDForTrigger] trig = " << trig << endl;
      cout << "[Event::IsPDForTrigger] PD = " << PD << endl;
      exit(ENODATA);
      return false;
    }
  }

  cout << "[Event::IsPDForTrigger] trig = " << trig << endl;
  cout << "[Event::IsPDForTrigger] PD = " << PD << endl;
  exit(ENODATA);
  return false;

}

void Event::SetMET(double pt, double phi){
  j_METVector.SetPtEtaPhiE(pt, 0, phi, pt);
}

