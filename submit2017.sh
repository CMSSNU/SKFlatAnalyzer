SKFlat.py -a singlelepton_analysis --skim SkimTree_HighPt1LJets -e 2017 -i SingleMuon -n 80 &
SKFlat.py -a singlelepton_analysis --skim SkimTree_HighPt1LJets -e 2017 -i SingleElectron -n 80 &
SKFlat.py -a singlelepton_analysis --skim SkimTree_HighPt1LJets -e 2017 -l SampleList/80cores.dat -n 80 &
SKFlat.py -a singlelepton_analysis --skim SkimTree_HighPt1LJets -e 2017 -l SampleList/30cores.dat -n 30 &
SKFlat.py -a singlelepton_analysis --skim SkimTree_HighPt1LJets -e 2017 -l SampleList/10cores.dat -n 10 &
SKFlat.py -a singlelepton_analysis --skim SkimTree_HighPt1LJets -e 2017 -l SampleList/5cores.dat -n 5 &
SKFlat.py -a singlelepton_analysis --skim SkimTree_HighPt1LJets -e 2017 -l SampleList/1cores.dat -n 1 &
SKFlat.py -a singlelepton_analysis --skim SkimTree_HighPt1LJets -e 2017 -l SampleList/signal.dat -n 1 &
