import csv
import os
import subprocess
import sys
from multiprocessing import Pool
import time

import matplotlib.pyplot as plt
import numpy as np

# Set number of runs per evaluation
nRuns=100

# Set number of simultaneous processes
nProcesses = 100

# Function that create the run folder and executes the ns3 command
def start_simulation(runParams):
    (
        nUes,
        deployment,
        d,
        relayDensity,
        nPaths,
        pathFindingAlgo,
        distanceConstraint,
        rsrpConstraint,
        maxNumRelays,
        allowReverse,
        allowRelayEndNode,
        packetSizeBe,
        dataRate,
        bidirectional,
        mu, 
        maxNumTx, 
        dynSch, 
        sensing,
        rngRun,
        outputDir
    ) = runParams

    runDir= ("%s/Run%d" % (outputDir, rngRun)) 
    print(f"Run dir: %s" % (runDir))
    try:
      os.mkdir(runDir)
    except FileExistsError:
      print("folder %s already exists, overwriting..." % runDir)

    #Run the simulation
    run_command = f"./ns3 run --cwd={runDir} nr-prose-u2u-multihop -- --RngRun={rngRun} --nUes={nUes} --deployment={deployment} --d={d} --relayDensity={relayDensity} --nPaths={nPaths}  --pathFindingAlgo={pathFindingAlgo} --distanceConstraint={distanceConstraint}  --rsrpConstraint={rsrpConstraint} --maxNumRelays={maxNumRelays} --allowReverse={allowReverse} --allowRelayEndNode={allowRelayEndNode} --packetSizeBe={packetSizeBe} --dataRate={dataRate}  --bidirectional={bidirectional} --mu={mu} --maxNumTx={maxNumTx} --dynSch={dynSch} --sensing={sensing}"
    print (run_command)
    with open(runDir + "/output.txt", "w") as f:
        output = subprocess.run(run_command.split(), stdout=f, stderr=subprocess.STDOUT)
        f.close()
    #Run the processing script
    subprocess.run(["python3", "nr-prose-u2u-multihop_plotSimStats.py", runDir], check=True) 

    
# Function that creates the evaluation folder and calls start_simulation for every run
def start_evaluation(params):
    (
        campaignName,
        nUes,
        deployment,
        d,
        relayDensity,
        nPaths,
        pathFindingAlgo,
        distanceConstraint,
        rsrpConstraint,
        maxNumRelays,
        allowReverse,
        allowRelayEndNode,
        packetSizeBe,
        dataRate,
        bidirectional,
        mu, 
        maxNumTx, 
        dynSch, 
        sensing
    ) = params

    print (params)
    outputDir = (
        "%s_nUes-%d_depl-%s_d-%d_relDen-%d_nPaths-%d_pAlgo-%s_distC-%d_rsrpC-%s_maxNRel-%d_reverse-%s_relEndN-%s_pktSz-%d_dRate-%d_biD-%s_mu-%d_maxNTx-%d_dynSch-%s_sens-%s"
        % (
            campaignName,
            nUes,
            deployment,
            d,
            relayDensity*100,
            nPaths,
            pathFindingAlgo,
            distanceConstraint,
            rsrpConstraint,
            maxNumRelays,
            allowReverse,
            allowRelayEndNode,
            packetSizeBe,
            dataRate,
            bidirectional,
            mu, 
            maxNumTx, 
            dynSch, 
            sensing
        )
    )
    print(f"Output dir: %s" % (outputDir))
    try:
        os.mkdir(outputDir)
    except FileExistsError:
        print("folder %s already exists, overwriting..." % outputDir)
        #Run simulations 
        
    # Create list with all simulations to run. Each element is a tuple of simulation parameters values
    allSims = []
    for rngRun in range(1, nRuns+1, 1):
        runParams = (
            nUes,
            deployment,
            d,
            relayDensity,
            nPaths,
            pathFindingAlgo,
            distanceConstraint,
            rsrpConstraint,
            maxNumRelays,
            allowReverse,
            allowRelayEndNode,
            packetSizeBe,
            dataRate,
            bidirectional,
            mu, 
            maxNumTx, 
            dynSch, 
            sensing,
            rngRun, #For the simulation
            outputDir #For the simulation
       )
        allSims.append(runParams)

    #  Create a Pool object with "nProcesses" number of simultaneous processes/simulations
    pool = Pool(processes=nProcesses)

    #  Run the simulations
    print("Running simulations...")
    for _ in pool.imap_unordered(start_simulation, allSims):
        pass
    pool.close()
    pool.join()

    return outputDir

#######################################################################################################
########################################### Main script ###############################################
#######################################################################################################

start_time = time.time()

# Inital compilation
print("Running ./ns3 build...")
output = subprocess.run(["./ns3", "show", "config"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if output.returncode:
    print("Error:  Is the project configured?  Run ./ns3 configure ... first")
    print(output.stdout.decode("utf-8"))
    print(output.stderr.decode("utf-8"))
    sys.exit(1)
output = subprocess.run(["./ns3", "build"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if output.returncode:
    print("Error:  The build failed; fix it to continue")
    print(output.stdout.decode("utf-8"))
    print(output.stderr.decode("utf-8"))
    sys.exit(1)


#Baseline parameters 
nUes=10
deployment="SqrRand"
d=1500
relayDensity=1.0
nPaths=1
pathFindingAlgo="DiscoveryShortest"
distanceConstraint=0
rsrpConstraint=3.16e-12
maxNumRelays=1
allowReverse=False
allowRelayEndNode=True
packetSizeBe=60
dataRate=24
bidirectional=False
mu=0
maxNumTx=1
dynSch=False
sensing=True

#################################
##### DOD 1MB configuration #####
#nUes=20
#d=800
#packetSizeBe=2500
#dataRate=1000
#rsrpConstraint=1e-10
#MCS=14 in the scenario
#20MHz bandwidht in the scenario
#################################



campaignName=""
for i in [1]:

    if i == 1:
        campaignName=f"z_1_BL"

    # Start evaluations for the campaign

    #Eval 1
    maxNumRelays=0
    params = (campaignName, nUes, deployment, d, relayDensity, nPaths, pathFindingAlgo, distanceConstraint, rsrpConstraint, maxNumRelays, allowReverse, allowRelayEndNode, packetSizeBe, dataRate, bidirectional, mu, maxNumTx, dynSch, sensing)
    evalFolder_1=start_evaluation (params)
    print(f"Eval folder: %s" % (evalFolder_1))

    #Eval 2
    maxNumRelays=1
    params = (campaignName, nUes, deployment, d, relayDensity, nPaths, pathFindingAlgo, distanceConstraint, rsrpConstraint, maxNumRelays, allowReverse, allowRelayEndNode, packetSizeBe, dataRate, bidirectional, mu, maxNumTx, dynSch, sensing)
    evalFolder_2=start_evaluation (params)
    print(f"Eval folder: %s" % (evalFolder_2))

    #Eval 3
    maxNumRelays=2
    params = (campaignName, nUes, deployment, d, relayDensity, nPaths, pathFindingAlgo, distanceConstraint, rsrpConstraint, maxNumRelays, allowReverse, allowRelayEndNode, packetSizeBe, dataRate, bidirectional, mu, maxNumTx, dynSch, sensing)

    evalFolder_3=start_evaluation (params)
    print(f"Eval folder: %s" % (evalFolder_3))

    #Eval 4
    maxNumRelays=3
    params = (campaignName, nUes, deployment, d, relayDensity, nPaths, pathFindingAlgo, distanceConstraint, rsrpConstraint, maxNumRelays, allowReverse, allowRelayEndNode, packetSizeBe, dataRate, bidirectional, mu, maxNumTx, dynSch, sensing)

    evalFolder_4=start_evaluation (params)
    print(f"Eval folder: %s" % (evalFolder_4)) 

    #Eval 5
    maxNumRelays=4
    params = (campaignName, nUes, deployment, d, relayDensity, nPaths, pathFindingAlgo, distanceConstraint, rsrpConstraint, maxNumRelays, allowReverse, allowRelayEndNode, packetSizeBe, dataRate, bidirectional, mu, maxNumTx, dynSch, sensing)

    evalFolder_5=start_evaluation (params)
    print(f"Eval folder: %s" % (evalFolder_5))

    # Create campaign plots with all evaluations
    scriptParam = [
        "python3",
        "nr-prose-u2u-multihop_plotEvalStats.py",
        "--campaignName", campaignName,
        "--campaignLabel", "maxNRelays",
        "--eval", evalFolder_1, "0",
        "--eval", evalFolder_2, "1",
        "--eval", evalFolder_3, "2",
        "--eval", evalFolder_4, "3",
        "--eval", evalFolder_5, "4"
    ]

    subprocess.run(scriptParam, check=True)

end_time = time.time()
duration_s = end_time - start_time
duration_m = duration_s / 60.0
print(f"The script took {duration_m} minutes to run.")
