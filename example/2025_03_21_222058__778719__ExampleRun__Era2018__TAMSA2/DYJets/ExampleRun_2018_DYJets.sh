#!/bin/bash
SECTION=`printf $1`
WORKDIR=`pwd`
Trial=0

#### make sure use C locale
export LC_ALL=C

#### Don't make root history
export ROOT_HIST=0

#### use cvmfs for root ####
export CMS_PATH=/cvmfs/cms.cern.ch
source $CMS_PATH/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
export cmsswrel=cmssw/CMSSW_15_0_1
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/$cmsswrel/src
echo "@@@@ SCRAM_ARCH = "$SCRAM_ARCH
echo "@@@@ cmsswrel = "$cmsswrel
echo "@@@@ scram..."
eval `scramv1 runtime -sh`
cd -
#source /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/$cmsswrel/external/$SCRAM_ARCH/bin/thisroot.sh

### modifying LD_LIBRARY_PATH to use libraries in base_rundir
export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH|sed 's@'$SKFlat_WD'/lib@/data6/Users/choij/SKFlatRunlog//2025_03_21_222058__778719__ExampleRun__Era2018__TAMSA2//lib@')
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$SKFlat_WD/DataFormats/include:$SKFlat_WD/AnalyzerTools/include:$SKFlat_WD/Analyzers/include

while [ "$Trial" -lt 3 ]; do
  echo "#### running ####"
  echo "root -l -b -q /data6/Users/choij/SKFlatRunlog//2025_03_21_222058__778719__ExampleRun__Era2018__TAMSA2/DYJets//run_${SECTION}.C"
  root -l -b -q /data6/Users/choij/SKFlatRunlog//2025_03_21_222058__778719__ExampleRun__Era2018__TAMSA2/DYJets//run_${SECTION}.C 2> err.log 
  EXITCODE=$?
  if [ "$EXITCODE" -eq 5 ]; then
    echo "IO error occured.. running again in 300 seconds.."
    Trial=$((Trial+=1))
    sleep 300
  else
    break
  fi
done

if [ "$EXITCODE" -ne 0 ]; then
  echo "ERROR errno=$EXITCODE" >> err.log
fi

cat err.log >&2
exit $EXITCODE

