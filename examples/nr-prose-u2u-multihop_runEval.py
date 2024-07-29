import csv
import os
import subprocess
import sys
from multiprocessing import Pool
import time

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

# Set number of runs per evaluation
nRuns=20

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


#Function used for the higher level postprocesing
def plot_with_error_bars(csv_file, x_column, y_column, set_column, campaignPrefix):
    # Read data from CSV
    with open(csv_file, mode='r') as file:
        reader = csv.DictReader(file)
        
        # Initialize data containers
        data_by_set = {}

        for row in reader:
            set_value = row[set_column]
            x_value = float(row[x_column])
            mean, ci = map(float, row[y_column].strip("()").split(','))
            nUEs = int(row['nUEs'])
            interUeDistance = float(row['interUeDistance'])
            d_value = (nUEs - 1) * interUeDistance
            
            if set_value not in data_by_set:
                data_by_set[set_value] = {'x': [], 'mean': [], 'ci': [], 'd': []}

            data_by_set[set_value]['x'].append(x_value)
            data_by_set[set_value]['mean'].append(mean)
            data_by_set[set_value]['ci'].append(ci)
            data_by_set[set_value]['d'].append(d_value)
            
    # Sort sets by numeric order if possible
    def sort_key(item):
        try:
            return float(item)
        except ValueError:
            return item

    sorted_sets = sorted(data_by_set.keys(), key=sort_key)
    

    # Set default font sizes
    plt.rcParams.update({'font.size': 12, 'legend.fontsize': 'medium'})

    # Create subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)
    
    # Colors for each set
    colors = plt.cm.tab10(np.linspace(0, 1, len(sorted_sets)))
    set_colors = {set_value: color for set_value, color in zip(sorted_sets, colors)}
    
    # Plot data for each sorted set
    for set_value in sorted_sets:
        data = data_by_set[set_value]
        sorted_indices = np.argsort(data['x'])
        x_sorted = np.array(data['x'])[sorted_indices]
        mean_sorted = np.array(data['mean'])[sorted_indices]
        ci_sorted = np.array(data['ci'])[sorted_indices]
        d_sorted = np.array(data['d'])[sorted_indices]
        
        # Plot first subplot
        ax1.errorbar(x_sorted, mean_sorted, yerr=ci_sorted, fmt='o-', label=set_value, color=set_colors[set_value])
        
        # Plot second subplot
        ax2.plot(x_sorted, d_sorted, 'x-', label=f"{set_value}", color=set_colors[set_value])

    # Additional settings for avgFlowLossRatio
    # Check if we are plotting avgFlowLossRatio and need to add horizontal lines
    if y_column == 'avgFlowLossRatio':
        ax1.set_yscale('log')
        ax1.set_ylim([0.0001,1])
        lines = []  # List to keep track of the line objects for the legend
        for y_line, style in zip([0.02, 0.05, 0.1], ['--', '-.', ':']):
            lines.append(ax1.axhline(y=y_line, color='gray', linestyle=style))

    # After all plotting is done, handle the legend
    handles, labels = ax1.get_legend_handles_labels()
    if y_column == 'avgFlowLossRatio':
        # Add the horizontal line legend entries at the end
        for i, y_line in enumerate([0.02, 0.05, 0.1]):
            # Create custom handles for the horizontal lines
            handles.append(lines[i])
            labels.append(f'y={y_line}')

    # Now set the updated handles and labels for the legend of ax1
    ax1.legend(handles, labels, title=set_column, bbox_to_anchor=(1.02, 1), loc='upper left')

#    ax1.set_title(f"{y_column} by {x_column} for each {set_column}")
    ax1.set_xlabel(x_column)
    ax1.set_ylabel(y_column)
    ax1.grid(True, which='major', color='#cccccc', linestyle='-', linewidth=0.5)
    ax1.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)

#    ax2.set_title(f'Calculated d = (nUEs-1) * interUeDistance by {set_column}')
    ax2.set_xlabel(x_column)
    ax2.set_ylabel('d (m)')
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(1000))
    ax2.legend(title=set_column, bbox_to_anchor=(1.02, 1), loc='upper left')
    ax2.grid(True, which='major', color='#cccccc', linestyle='-', linewidth=0.5)

    plt.tight_layout()
    plt.tight_layout()
    png_file_name = f"{campaignPrefix}_X-{x_column}_Y-{y_column}_Set-{set_column}.png"
    plt.savefig(png_file_name)
    plt.close()  # Close the plot



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


# Baseline campaign varying maxNumRelays  
""" 
#TODO:Postprocessing not working for this one. DEBUG!
campaignName=""
for i in [1]:

    evalFolders = {}  # Dictionary to store maxNumRelays and corresponding eval folder

    if i == 1:
        campaignName=f"z_0_BL"

    # Start evaluations for the campaign

    for maxNumRelays in [0,1,2,3,4]: 
        params = (campaignName, nUes, deployment, d, relayDensity, nPaths, pathFindingAlgo, distanceConstraint, rsrpConstraint, maxNumRelays, allowReverse, allowRelayEndNode, packetSizeBe, dataRate, bidirectional, mu, maxNumTx, dynSch, sensing)
        evalFolder = start_evaluation(params)
        evalFolders[maxNumRelays] = evalFolder
        print(f"Eval folder for maxNumRelays {maxNumRelays}: {evalFolder}")

    # Now, prepare the parameters for the script
    scriptParam = ["python3", "nr-prose-u2u-multihop_plotEvalStats.py", "--campaignName", campaignName, "--campaignLabel", "maxNRelays"]

    for maxNumRelays, folder in evalFolders.items():
        scriptParam.extend(["--eval", folder, str(maxNumRelays)])

    subprocess.run(scriptParam, check=True) """
 



""" # NPSTC campaign varying maxNumRelays  
#packetSizeBe=320 #1 pcket every 10 ms
packetSizeBe=640 #Larger packet less often (20ms)
dataRate=256
dynSch=True
sensing=False
campaignName=""
rsrpConstraint=1e-11 # 110 dBm
maxRel=4

for i in [3, 6]:

    evalFolders = {}  # Dictionary to store maxNumRelays and corresponding eval folder

    if i == 1:
        nUes=46
        d=1200
        maxRel=4
        nPaths=1
        campaignName=f"z_30_NPSTC_nUes-{nUes}_d-{d}_nPaths-{nPaths}"

    if i == 2:
        #nUes=184
        nUes=46
        d=2400
        maxRel=8
        nPaths=1
        campaignName=f"z_31_NPSTC_nUes-{nUes}_d-{d}_nPaths-{nPaths}"

    if i == 3:
 #       nUes=750
        nUes=184
        d=4750
        maxRel=16
        nPaths=1
        campaignName=f"z_36_NPSTC_nUes-{nUes}_d-{d}_nPaths-{nPaths}"

    if i == 4:
        nUes=46
        d=1200
        maxRel=4
        nPaths=3
        campaignName=f"z_33_NPSTC_nUes-{nUes}_d-{d}_nPaths-{nPaths}"

    if i == 5:
        #nUes=184
        nUes=46
        d=2400
        maxRel=8
        nPaths=3
        campaignName=f"z_34_NPSTC_nUes-{nUes}_d-{d}_nPaths-{nPaths}"

    if i == 6:
#       nUes=750
        nUes=184
        d=4750
        maxRel=16
        nPaths=3
        campaignName=f"z_37_NPSTC_nUes-{nUes}_d-{d}_nPaths-{nPaths}" 



    # Start evaluations for the campaign

    for maxNumRelays in range (maxRel +1): #Iterate from 0 to maxRel
        params = (campaignName, nUes, deployment, d, relayDensity, nPaths, pathFindingAlgo, distanceConstraint, rsrpConstraint, maxNumRelays, allowReverse, allowRelayEndNode, packetSizeBe, dataRate, bidirectional, mu, maxNumTx, dynSch, sensing)
        evalFolder = start_evaluation(params)
        evalFolders[maxNumRelays] = evalFolder
        print(f"Eval folder for maxNumRelays {maxNumRelays}: {evalFolder}")

    # Now, prepare the parameters for the script
    scriptParam = ["python3", "nr-prose-u2u-multihop_plotEvalStats.py", "--campaignName", campaignName, "--campaignLabel", "maxNRelays"]

    for maxNumRelays, folder in evalFolders.items():
        scriptParam.extend(["--eval", folder, str(maxNumRelays)])

    subprocess.run(scriptParam, check=True)"""





#Campaing for seeing delay depending on inter UE distance
deployment="Linear"
relayDensity=1.0
nPaths=1
pathFindingAlgo="DistShortest"
mu=0
dynSch=False

campaignPrefix="z_2"
campaignName=""
#for nUes in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
for nUes in [10]:

    evalFolders = {}  # Dictionary to store maxNumRelays and corresponding eval folder
    campaignName=f"{campaignPrefix}_nUes-{nUes}"

    # Start evaluations for the campaign
#    for iud in [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]: 
    for iud in [500, 800]: 
#    for iud in [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]: 
        d=iud * (nUes -1) 
        distanceConstraint=iud+1
        maxNumRelays=nUes

        params = (campaignName, nUes, deployment, d, relayDensity, nPaths, pathFindingAlgo, distanceConstraint, rsrpConstraint, maxNumRelays, allowReverse, allowRelayEndNode, packetSizeBe, dataRate, bidirectional, mu, maxNumTx, dynSch, sensing)
        evalFolder = start_evaluation(params)
        evalFolders[iud] = evalFolder
        print(f"Eval folder for iud {iud}: {evalFolder}")

    # Now, prepare the parameters for the script
    scriptParam = ["python3", "nr-prose-u2u-multihop_plotEvalStats.py", "--campaignName", campaignName, "--campaignLabel", "interUeDistance"]

    for iud, folder in evalFolders.items():
        scriptParam.extend(["--eval", folder, str(iud)])

    subprocess.run(scriptParam, check=True)   


#The data for the following plots has been curated by hand:
# campaignPrefix="z_20"
# plot_with_error_bars('z_20_nUes-All_evalSimStats_meanCI.csv', 'interUeDistance', 'avgFlowMeanDelay', 'nUEs', campaignPrefix)
# plot_with_error_bars('z_20_nUes-All_evalSimStats_meanCI.csv', 'interUeDistance', 'avgFlowLossRatio', 'nUEs', campaignPrefix)
# plot_with_error_bars('z_20_nUes-All_evalSimStats_meanCI.csv', 'nUEs', 'avgFlowMeanDelay', 'interUeDistance', campaignPrefix)
# plot_with_error_bars('z_20_nUes-All_evalSimStats_meanCI.csv', 'nUEs', 'avgFlowLossRatio', 'interUeDistance', campaignPrefix)
# plot_with_error_bars('z_20_nUes-All_evalSimPhyStatsGlobal_meanCI.csv', 'interUeDistance', 'nHdRxData', 'nUEs', campaignPrefix)
# plot_with_error_bars('z_20_nUes-All_evalSimPhyStatsGlobal_meanCI.csv', 'interUeDistance', 'nCorruptRxData', 'nUEs', campaignPrefix)
# plot_with_error_bars('z_20_nUes-All_evalSimPhyStatsGlobal_meanCI.csv', 'interUeDistance', 'nIgnoredData', 'nUEs', campaignPrefix)


end_time = time.time()
duration_s = end_time - start_time
duration_m = duration_s / 60.0
print(f"The script took {duration_m} minutes to run.")
