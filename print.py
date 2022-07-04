import os

skim = "SkimTree_HighPt1LJets"
year = "2017"

datasets = ["SingleElectron", "SingleMuon"]

for dataset in datasets:
    periods = os.listdir("/gv0/DATA/SKFlat/Run2UltraLegacy_v2/" + year + "/DATA_" + skim + "/" + dataset)
    for period in periods:
        data_path = "/gv0/DATA/SKFlat/Run2UltraLegacy_v2/" + year + "/DATA_" + skim + "/" + dataset + "/" + period + "/*/*root"
        os.system("ls -1 " + data_path + " > data/Run2UltraLegacy_v2/" + year + "/Sample/ForSNU/" + skim + "_" + dataset + "_" + period.replace("period", "") + ".txt")

mcsets = os.listdir("/gv0/DATA/SKFlat/Run2UltraLegacy_v2/" + year + "/MC_" + skim +"/")
for mcset in mcsets:
    mc_path = "/gv0/DATA/SKFlat/Run2UltraLegacy_v2/" + year + "/MC_" + skim + "/" + mcset + "/*/*root"
    os.system("grep -r '" + mcset + "' data/Run2UltraLegacy_v2/" + year + "/Sample/CommonSampleInfo/ > temp.txt")
    with open("temp.txt") as temp:
        line = temp.read()
        alias = line.strip().split(":")[1].split(" ")[0].split("\t")[0]
    
    os.system("ls -1 " + mc_path + " > data/Run2UltraLegacy_v2/" + year + "/Sample/ForSNU/" + skim + "_" + alias + ".txt")
