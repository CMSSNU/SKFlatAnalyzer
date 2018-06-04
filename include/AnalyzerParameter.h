#include "TString.h"

class AnalyzerParameter{

public:

  TString Name;

  bool MCCorrrectionIgnoreNoHist;

  TString Electron_Tight_ID, Electron_Loose_ID, Electron_Veto_ID;
  TString Electron_ID_SF_Key;
  double Electron_Tight_RelIso, Electron_Loose_RelIso, Electron_Veto_RelIso;

  TString Muon_Tight_ID, Muon_Loose_ID, Muon_Veto_ID;
  TString Muon_ID_SF_Key, Muon_ISO_SF_Key, Muon_Trigger_SF_Key;
  double Muon_Tight_RelIso, Muon_Loose_RelIso, Muon_Veto_RelIso;

  TString Jet_ID, FatJet_ID;

  void Clear();

  AnalyzerParameter();
  ~AnalyzerParameter();

};

void AnalyzerParameter::Clear(){

  Name = "";

  MCCorrrectionIgnoreNoHist = false;

  Electron_Tight_ID = "";
  Electron_Loose_ID = "";
  Electron_Veto_ID = "";
  Electron_ID_SF_Key = "";
  Electron_Tight_RelIso = 999.;
  Electron_Loose_RelIso = 999.;
  Electron_Veto_RelIso = 999.;

  Muon_Tight_ID = "";

  Muon_Loose_ID = "";
  Muon_Veto_ID = "";
  Muon_ID_SF_Key = "";
  Muon_ISO_SF_Key = "";
  Muon_Trigger_SF_Key = "";
  Muon_Tight_RelIso = 999.;
  Muon_Loose_RelIso = 999.;
  Muon_Veto_RelIso = 999.;

  Jet_ID = "";
  FatJet_ID = "";

}

AnalyzerParameter::AnalyzerParameter(){

  Name = "Default";

  MCCorrrectionIgnoreNoHist = false;

  Electron_Tight_ID = "passTightID";
  Electron_Loose_ID = "passLooseID";
  Electron_Veto_ID = "passVetoID";
  Electron_ID_SF_Key = "passTightID";

  Muon_Tight_ID = "POGTightWithTightIso";
  Muon_Loose_ID = "POGLoose";
  Muon_Veto_ID = "POGLoose";
  Muon_ID_SF_Key = "NUM_TightID_DEN_genTracks";
  Muon_ISO_SF_Key = "NUM_TightRelIso_DEN_TightIDandIPCut";
  Muon_Trigger_SF_Key = "POGTight";

  Jet_ID = "HN";
  FatJet_ID = "HN";

}

AnalyzerParameter::~AnalyzerParameter(){
}
