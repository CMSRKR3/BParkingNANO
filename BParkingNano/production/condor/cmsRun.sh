#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a bash script, use .sh instead of .csh
export X509_USER_PROXY=${7}
voms-proxy-info -all
voms-proxy-info -all -file ${7}
### for case 1. EOS have the following line, otherwise remove this line in case 2.
#Arguments 1)eosDir 2)jobName 3)rel 4)iJob 5)jobCfg 6)inputFile 7)fileName

# You need to move transfered files to the working area
cp ${1}${2}.tgz .
tar -xf ${2}.tgz
rm ${2}.tgz
mv ${5} ${3}/src/PhysicsTools/BParkingNano
mv filelist_${2}.tgz ${3}/src/PhysicsTools/BParkingNano

# export $SCRAM_ARCH=slc7_amd64_gcc10
cd ${3}/src/
scramv1 b ProjectRename -j 8
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

cd PhysicsTools/BParkingNano
tar -xf filelist_${2}.tgz
ls ./
touch log.txt
cmsRun ${5} inputFiles_load=${6} 2>&1 | tee log.txt
cp log.txt ${1}/log_${4}.txt
cp BParkingNANO_Run3_mc_124X.root ${1}/myNano_${4}.root
rm BParkingNANO_Run3_mc_124X.root
cd ${_CONDOR_SCRATCH_DIR}
rm -rf ${3}