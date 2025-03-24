# SKFlatAnalyzer

## Manual

https://jskim.web.cern.ch/jskim/SKFlat/Manual/Manual_SKFlat.pdf

## Where to put the analyzer
TAMSA 1/2 : /data6/Users/$USER/

KISTI : ~/ # home directory

KNU :  ~/scartch/

## First time setup
```
#### When first time gie clone, use the option "--recursive" to initiate the submodules
git clone --recursive git@github.com:CMSSNU/SKFlatAnalyzer.git
cd SKFlatAnalyzer
#### add your remote repo
git remote add <nickname> git@github.com:<gitaccount>/SKFlatAnalyzer.git
git checkout <your working branch>

#### First time setup script
source bin/FirstTimeSetup.sh 
source setup.sh

#### You have to edit user info
#### First, copy the temply using the command below
cp $SKFlat_WD/python/UserInfo_template.py $SKFlat_WD/python/UserInfo_${USER}.py 
#### Then, edit $SKFlat_WD/python/UserInfo_${USER}.py

#### Now, run setup script.
#### This should be done for every new shell
source setup.sh
```
Compile
```
make clean
make
```

## Test job
```bash
SKFlat.py -a ExampleRun -i DYJets -n 50 -e 2018 &
```

## Making a new Ananlyzer
```bash
cd script/MakeCycleSkeleton/
```
Then, run
```bash
python MakeCycleSkeleton.py NewAnalyzer # <- put new analyzer name
```
It will print below lines (execute the lines) :
```bash
mv NewAnalyzer.h $SKFlat_WD/Analyzers/include/
mv NewAnalyzer.C $SKFlat_WD/Analyzers/src/
```

Then, add
```bash
#pragma link C++ class NewAnalyzer+;
```
in Analyzers/include/Analyzers_LinkDef.h

## Detailed descriptions

Look Analyzers/src/ExampleRun.C

## Adding samples
To add a sample, you should add two files.  
* data/$SKFlatV/$ERA/Sample/ForSNU/$ALIAS.txt: list of file paths.  
* data/$SKFlatV/$ERA/Sample/CommonSampleInfo/$ALIAS.txt: alias, DAS name, cross section, nevent, sum(sign) and sum(weight).  

And one file should be edited  
* data/$SKFlatV/$ERA/Sample/SampleSummary_*.txt: CommonSampleInfo files in one file. This file is not actually used by SKFlatAnalyzer. It is just a summary for users.  

You can do it manually, or use scripts as below. The scripts use SampleSummary files for all SKFlat versions to find alias and cross section. If the sample is never used before the scripts will ask you alias and cross section. If it is annoying you can make temporal SampleSummary file like data/$SKFlatV/$ERA/Sample/SampleSummary_temp.txt which containing alias, DAS name and cross section (other information is not needed).
1. Make the file list file using bin/UpdateSampleForSNU.sh script.
```
./bin/UpdateSampleForSNU.sh $SAMPLEDIRECTORY
```
2. Run GetEffLumi analyzer to get nevent, sum(sign) and sum(weight)
```
SKFlat.py -a GetEffLumi -e $ERA -n 20 -i $SAMPLEALIAS
```
3. Make CommonSampleInfo file using bin/UpdateCommonSampleInfo.sh script.
```
./bin/UpdateCommonSampleInfo.sh
```
4. Update SampleSummary file using Summarize.py
```
cd data/$SKFlatV/$ERA/Sample
python Summarize.py
```

# Tips

## Making PR

Start from the CMSSNU's master branch of CMSSNU when making pull request to prevent anoying conflicts.

