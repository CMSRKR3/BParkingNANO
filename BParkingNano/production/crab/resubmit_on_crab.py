import os

for job in os.listdir("/afs/cern.ch/work/x/xuyan/RKProj/BParkingNanoOfficial/CMSSW_13_1_0/src/PhysicsTools/BParkingNano/production/crab/BParkingNANO_2024May02"):
    print(job)
    os.system(f"crab resubmit BParkingNANO_2024May02/{job}")
