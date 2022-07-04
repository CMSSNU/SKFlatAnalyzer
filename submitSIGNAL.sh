SKFlat.py -a signal_study --skim SkimTree_HighPt1LJets -e 2017 -l SampleList/QCD.dat -n 10 &
SKFlat.py -a signal_study -e 2017 -l SampleList/signal.dat -n 1 &
SKFlat.py -a signal_study --skim SkimTree_HighPt1LJets -e 2017 -i WJets_MG -n 40 &
