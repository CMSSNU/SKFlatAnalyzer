R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so)


void run_6(){

  ExampleRun m;

  m.SetTreeName("recoTree/SKFlat");

  m.LogEvery = 1000;
  m.MCSample = "DYJets";
  m.IsDATA = false;
  m.xsec = 6077.22;
  m.sumSign = 132322003.0;
  m.sumW = 3.342410636064e+12;
  m.IsFastSim = false;
  m.SetEra("2018");
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_745.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_746.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_747.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_748.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_749.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_750.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_751.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_752.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_753.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_754.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_755.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_756.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_757.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_758.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_759.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_760.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_761.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_762.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_763.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_764.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_765.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_766.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_767.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_768.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_769.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_770.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_771.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_772.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_773.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_774.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_775.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_776.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_777.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_778.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_779.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_780.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_781.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_782.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_783.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_784.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_785.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_786.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_787.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_788.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_789.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_790.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_791.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_792.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_793.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_794.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_795.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_796.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_797.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_798.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_799.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_800.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_801.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_802.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_803.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_804.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_805.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_806.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_807.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_808.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_809.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_810.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_811.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_812.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_813.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_814.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_815.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_816.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_817.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_818.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_819.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_820.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_821.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_822.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_823.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_824.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_825.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_826.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_827.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_828.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_829.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_830.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_831.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_832.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_833.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_834.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_835.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_836.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_837.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_838.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_839.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_840.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_841.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_842.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_843.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_844.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_845.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_846.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_847.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_848.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_849.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_850.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_851.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_852.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_853.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_854.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_855.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_856.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_857.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_858.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_859.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_860.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_861.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_862.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_863.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_864.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_865.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_866.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_867.root")) exit(EIO);
  if(!m.AddFile("/gv0/DATA/SKFlat/Run2UltraLegacy_v3/2018/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220621_003647/0000/SKFlatNtuple_2018_MC_868.root")) exit(EIO);
  m.SetOutfilePath("hists.root");
  m.MaxEvent=m.fChain->GetEntries()/100.0;
  m.Init();
  m.initializeAnalyzer();
  m.initializeAnalyzerTools();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}
