#!/usr/bin/env python
import os, re, sys, commands, math, time, calendar

print '\nSTART\n'
ts = calendar.timegm(time.gmtime())

# joblist = ["BuToKEE","BdToK0starEE_postEE","BdToK0starEE_postEE_ext1","BdToK0starEE_preEE","BdToK0starEE_preEE_ext1"]
joblist = ["DataF"]

jobdic = {
           "DataF":"ParkingDoubleElectronLowMass0_Run2022F.list",
           "BuToKEE":"BuToKEE_SoftQCD_TuneCP5_13p6TeV_pythia8-evtgen.list",
           "BdToK0starEE_postEE":"BdToK0starEE_K0starToKPi_SoftQCD_TuneCP5_13p6TeV_pythia8-evtgen_postEE.list",
           "BdToK0starEE_postEE_ext1":"BdToK0starEE_K0starToKPi_SoftQCD_TuneCP5_13p6TeV_pythia8-evtgen_postEE_ext1.list",
           "BdToK0starEE_preEE":"BdToK0starEE_K0starToKPi_SoftQCD_TuneCP5_13p6TeV_pythia8-evtgen_preEE.list",
           "BdToK0starEE_preEE_ext1":"BdToK0starEE_K0starToKPi_SoftQCD_TuneCP5_13p6TeV_pythia8-evtgen_preEE_ext1.list",
}

for jobName in joblist:
    print(jobName, jobdic[jobName])
    jobCfg = "run_nano_cfg.py"
    jobScript = "cmsRun.sh"
    rel = "CMSSW_13_1_0"
    eosDir = "/eos/home-x/xuyan/RKProj/RKAnalyzer/KStarllStudy/Data/" + jobName + "_" + str(ts) + "/"
    isMC = "True"
    if "Data" in jobName:
        isMC = "False"
    print("Running on MC: ",isMC)

    rootDir = os.environ["CMSSW_BASE"] + "/src/PhysicsTools/BParkingNano/production/condor/"
    jobDir = rootDir + jobName + "_" + str(ts) + "/"
    files_batch = 10

    fileList = rootDir + jobdic[jobName]

    ret = os.system("rm filelist_tmp_*.list")
    ret = os.system("cp ../../test/" + jobCfg + " ./")
    ret = os.system("mkdir -p " + jobDir)
    ret = os.system("mkdir -p " + eosDir)
    ret = os.system("mkdir -p " + jobDir + "out/")
    ret = os.system("mkdir -p " + jobDir + "err/")
    ret = os.system("mkdir -p " + jobDir + "log/")
    ret = os.chdir(os.environ["CMSSW_BASE"]+"/../")
    print('Tarballing ' + rel + "/ into " + jobName + ".tgz...")
    ret = os.system("tar --exclude='*.root' --exclude='ignore' --exclude='.git' " + "-zcf " + jobName + ".tgz " + rel)
    print 'Done!'
    ret = os.system("mv " + jobName + ".tgz " + eosDir) 
    ret = os.chdir(rootDir)
    ret = os.system("cp /tmp/x509up_u93529 /afs/cern.ch/user/x/xuyan/private/x509up/x509up_u93529")
    proxy_path = "/afs/cern.ch/user/x/xuyan/private/x509up/x509up_u93529"

    file1 = open(fileList, 'r')
    file_content = file1.readlines()
    chunks = len(file_content) / files_batch + 1

    ret = os.system("rm filelist_tmp_*.list")
    for ichunk in range(chunks):
        filelist_name = "filelist_tmp_%i.list" %(ichunk)
        filelist_tmp = open(filelist_name, 'w')
        for file in file_content[ichunk*files_batch:(ichunk+1)*files_batch]:
            filelist_tmp.write(file)
        filelist_tmp.close()
    ret = os.system("tar -cvf filelist_" + jobName + ".tgz filelist_tmp_*")

    with open(jobDir + jobName + '.jdl', 'w') as jdl:
        jdl.write("universe = vanilla\n")
        jdl.write("Executable = " + jobScript + "\n")
        jdl.write("Should_Transfer_Files = YES\n")
        jdl.write("WhenToTransferOutput = ON_EXIT\n")
        jdl.write("Transfer_Input_Files = " + jobScript + ", " + jobCfg + ", filelist_" + jobName + ".tgz" + "\n")
        jdl.write("Output = "    + jobDir + "out/$(ProcId).out\n")
        jdl.write("Error = "     + jobDir + "err/$(ProcId).err\n")
        jdl.write("Log = "       + jobDir + "log/$(ProcId).log\n")
        jdl.write("Arguments = " + eosDir + " " + jobName + " " + rel + " $(ProcId) " + jobCfg + " filelist_tmp_$(ProcId).list " + proxy_path + " " + isMC +"\n")
        jdl.write("MY.WantOS = \"el7\"\n")
        jdl.write("+JobFlavour = " + "\"tomorrow\"" + "\n")
        # jdl.write("+MaxRuntime = 28800\n")
        # jdl.write("on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n")
        # jdl.write("max_retries = 3\n")
        # jdl.write("requirements = Machine =!= LastRemoteHost\n")
        jdl.write("Queue " + str(chunks) + "\n")      

    os.system("condor_submit " + jobDir + jobName + ".jdl")
    print str(chunks) + " jobs submitted."
    os.system("condor_q")