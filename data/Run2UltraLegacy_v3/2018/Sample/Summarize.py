import os

os.system('ls -1 CommonSampleInfo/*.txt &> temp.txt');
lines = open('temp.txt').readlines()
os.system('rm temp.txt')

out_MC = open('SampleSummary_MC.txt','w')
out_Sig_PairN = open('SampleSummary_Signal_PairN.txt','w')
out_Sig_HNWR = open('SampleSummary_Signal_HNWR.txt','w')

FirstLine = '# alias PD xsec nmc sumw\n'
out_MC.write(FirstLine)
out_Sig_PairN.write(FirstLine)
out_Sig_HNWR.write(FirstLine)

for line in lines:

  line = line.strip('\n')
  alias = line.replace('CommonSampleInfo/','').replace('.txt','')

  infoline = open(line).readlines()[1]

  if alias.startswith("HNPairToJJJJ"):
    out_Sig_PairN.write(infoline)
  elif alias.startswith("WR_"):
    out_Sig_HNWR.write(infoline)
  else:
    out_MC.write(infoline)

out_MC.close()
out_Sig_PairN.close()
out_Sig_HNWR.close()
