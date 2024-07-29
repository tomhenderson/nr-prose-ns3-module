import csv
import os
import subprocess
import sys
from multiprocessing import Pool
import time
import signal
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import psutil  

# Set number of runs per evaluation
nRuns=1000

# Set number of simultaneous processes
nProcesses = 50

#Set maximum execution time per simulation
maxExecTime = 10*60 # s

def kill_proc_tree(pid, including_parent=True):  
    parent = psutil.Process(pid)
    children = parent.children(recursive=True)
    for child in children:
        child.kill()
    psutil.wait_procs(children, timeout=5)
    if including_parent:
        parent.kill()
        parent.wait(5)

# Function that create the run folder and executes the ns3 command
def start_simulation(runParams):
    (
        nUes,
        d,
        nPaths,
        routingAlgo,
        maxNumTx, 
        sensing,
        harqType,
        nDisrupt,
        interval,
        rngRun,
        outputDir
    ) = runParams

    runDir= ("%s/Run%d" % (outputDir, rngRun)) 
    #print(f"Run dir: %s" % (runDir))
    try:
      os.mkdir(runDir)
    except FileExistsError:
      print("folder %s already exists, overwriting..." % runDir)

     #Run the simulation
    run_command = f"./ns3 run --cwd={runDir} nr-prose-u2u-multihop-routing -- --RngRun={rngRun} --nUes={nUes} --d={d} --nPaths={nPaths}  --routingAlgo={routingAlgo} --maxNumTx={maxNumTx} --sensing={sensing} --harqType={harqType} --nDisrupt={nDisrupt} --interval={interval}"
    print (run_command)
    try:
        with open(runDir + "/output.txt", "w") as f:
            process = subprocess.Popen(run_command.split(), stdout=f, stderr=subprocess.STDOUT)
            start = time.time()
            while True:
                if process.poll() is not None:
                    break
                if time.time() - start > maxExecTime:
                    kill_proc_tree(process.pid)  # Kill the process tree
                    process.wait()  # Wait for process to terminate
                    raise subprocess.TimeoutExpired(run_command, maxExecTime)
                time.sleep(0.5)  # Sleep a short time to prevent busy waiting
    except subprocess.TimeoutExpired:
        print(f"Simulation in {runDir} exceeded the maximum execution time of {maxExecTime} seconds and was terminated.")
        return {'rngRun': rngRun, 'status': 'killed'}


    #Run the processing script
    subprocess.run(["python3", "routing2024_plotSimStats.py", runDir], check=True)
    return {'rngRun': rngRun, 'status': 'completed'}

    
# Function that creates the evaluation folder and calls start_simulation for every run
def start_evaluation(params):
    (
        campaignName,
        nUes,
        d,
        nPaths,
        routingAlgo,
        maxNumTx, 
        sensing,
        harqType,
        nDisrupt,
        interval
    ) = params

    outputDir = (
        "%s_nUes-%d_d-%d_nPaths-%d_rAlgo-%s_maxNTx-%d_sens-%s_harqT-%s_nDisr-%d_int-%d"
        % (
            campaignName,
            nUes,
            d,
            nPaths,
            routingAlgo,
            maxNumTx, 
            sensing,
            harqType,
            nDisrupt,
            interval
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
            d,
            nPaths,
            routingAlgo,
            maxNumTx, 
            sensing,
            harqType,
            nDisrupt,
            interval,
            rngRun, #For the simulation
            outputDir #For the simulation
       )
        allSims.append(runParams)

    #  Create a Pool object with "nProcesses" number of simultaneous processes/simulations
    pool = Pool(processes=nProcesses)

    #  Run the simulations
    print("Running simulations...")
    results = pool.imap_unordered(start_simulation, allSims)
    pool.close()
    pool.join()

    killed_sims = [result for result in results if result['status'] == 'killed']
    print(f"Number of killed simulations: {len(killed_sims)}")
    if killed_sims:
        print("Details of killed simulations:")
        for sim in killed_sims:
            print(f"Simulation {sim['rngRun']} was terminated.")

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
            
            if set_value not in data_by_set:
                data_by_set[set_value] = {'x': [], 'mean': [], 'ci': [], 'd': []}

            data_by_set[set_value]['x'].append(x_value)
            data_by_set[set_value]['mean'].append(mean)
            data_by_set[set_value]['ci'].append(ci)
            
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
    fig, ax1, = plt.subplots(figsize=(5, 5))

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
        
        # Plot 
        ax1.errorbar(x_sorted, mean_sorted, yerr=ci_sorted, fmt='o-', label=set_value, color=set_colors[set_value])
        
    ax1.set_ylim(bottom=0)

    # Additional settings for avgFlowLossRatio
    # Check if we are plotting avgFlowLossRatio and need to add horizontal lines
    if y_column == 'avgFlowLossRatio':
        ax1.set_ylim(bottom=0,top=0.4)


        #ax1.set_yscale('log')
        #ax1.set_ylim([0.0001,1])
        
#        lines = []  # List to keep track of the line objects for the legend
#        for y_line, style in zip([0.02, 0.05, 0.1], ['--', '-.', ':']):
#            lines.append(ax1.axhline(y=y_line, color='gray', linestyle=style))

    # After all plotting is done, handle the legend
    handles, labels = ax1.get_legend_handles_labels()
#    if y_column == 'avgFlowLossRatio':
        # Add the horizontal line legend entries at the end
#        for i, y_line in enumerate([0.02, 0.05, 0.1]):
            # Create custom handles for the horizontal lines
#            handles.append(lines[i])
 #           labels.append(f'y={y_line}')

    # Now set the updated handles and labels for the legend of ax1
    ax1.legend(handles, labels, bbox_to_anchor=(0.5, -0.15), loc='upper center')

#    ax1.set_title(f"{y_column} by {x_column} for each {set_column}")
    ax1.set_xlabel(x_column)
    ax1.set_ylabel(y_column)
    ax1.grid(True, which='major', color='#cccccc', linestyle='-', linewidth=0.5)
    ax1.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)


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
d=2000
nPaths=1
routingAlgo="OLSR"
maxNumTx=1 
sensing="True"
harqType="No"
nDisrupt=0
interval=4.0

campaignPrefix="z_02_sensing-True_HP25"
campaignName=""

#for routingAlgo in ["OLSR","BATMAN"]:
for routingAlgo in ["BATMAN"]:

    evalFolders = {}  # Dictionary to store maxNumRelays and corresponding eval folder
    campaignName=f"{campaignPrefix}-{routingAlgo}"
 
    for d in [1600, 2000, 2400, 2800, 3200]: 
        params = (campaignName, nUes, d, nPaths, routingAlgo, maxNumTx, sensing, harqType, nDisrupt, interval)
        evalFolder = start_evaluation(params)
        evalFolders[d] = evalFolder
        print(f"Eval folder for d {d}: {evalFolder}")    
    scriptParam = ["python3", "routing2024_plotEvalStats.py", "--campaignName", campaignName, "--campaignLabel", "d"]
    for d, folder in evalFolders.items():
        scriptParam.extend(["--eval", folder, str(d)])
    subprocess.run(scriptParam, check=True)   
    
#     for interval in [0.5, 1.0, 2.0, 4.0, 8.0]: 
# #    for interval in [1.0]:
#         params = (campaignName, nUes, d, nPaths, routingAlgo, maxNumTx, sensing, harqType, nDisrupt, interval)
#         evalFolder = start_evaluation(params)
#         evalFolders[interval] = evalFolder
#         print(f"Eval folder for d {interval}: {evalFolder}")    
#     scriptParam = ["python3", "routing2024_plotEvalStats.py", "--campaignName", campaignName, "--campaignLabel", "interval"]
#     for interval, folder in evalFolders.items():
#         scriptParam.extend(["--eval", folder, str(interval)])
#     subprocess.run(scriptParam, check=True)  



# The below csv files were curated by hand after all simulations finished.
# plot_with_error_bars('z_01_ALL_nDisrupt-0_evalSimStats_meanCI_HP25.csv', 'interval', 'avgFlowLossRatio', 'Campaign', 'z_01_nDisrupt-0_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-1_evalSimStats_meanCI_HP25.csv', 'interval', 'avgFlowLossRatio', 'Campaign', 'z_01_nDisrupt-1_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-2_evalSimStats_meanCI.csv', 'interval', 'avgFlowLossRatio', 'Campaign', 'z_01_nDisrupt-2')

# plot_with_error_bars('z_01_ALL_nDisrupt-0_evalSimStats_meanCI_HP25.csv', 'interval', 'avgFlowMeanDelay', 'Campaign', 'z_01_nDisrupt-0_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-1_evalSimStats_meanCI_HP25.csv', 'interval', 'avgFlowMeanDelay', 'Campaign', 'z_01_nDisrupt-1_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-2_evalSimStats_meanCI.csv', 'interval', 'avgFlowMeanDelay', 'Campaign', 'z_01_nDisrupt-2')

# plot_with_error_bars('z_01_ALL_nDisrupt-0_evalSimStats_meanCI_HP25.csv', 'interval', 'sysNRouteChanges', 'Campaign', 'z_01_nDisrupt-0_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-1_evalSimStats_meanCI_HP25.csv', 'interval', 'sysNRouteChanges', 'Campaign', 'z_01_nDisrupt-1_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-2_evalSimStats_meanCI.csv', 'interval', 'sysNRouteChanges', 'Campaign', 'z_01_nDisrupt-2')

# plot_with_error_bars('z_01_ALL_nDisrupt-0_evalSimStats_meanCI_HP25.csv', 'interval', 'avgMaxNHops', 'Campaign', 'z_01_nDisrupt-0_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-1_evalSimStats_meanCI_HP25.csv', 'interval', 'avgMaxNHops', 'Campaign', 'z_01_nDisrupt-1_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-2_evalSimStats_meanCI.csv', 'interval', 'avgMaxNHops', 'Campaign', 'z_01_nDisrupt-2')


# plot_with_error_bars('z_01_ALL_nDisrupt-0_evalSimStats_meanCI_HP25.csv', 'interval', 'sysRoutingOhBytes', 'Campaign', 'z_01_nDisrupt-0_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-1_evalSimStats_meanCI_HP25.csv', 'interval', 'sysRoutingOhBytes', 'Campaign', 'z_01_nDisrupt-1_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-2_evalSimStats_meanCI.csv', 'interval', 'sysRoutingOhBytes', 'Campaign', 'z_01_nDisrupt-2')

# plot_with_error_bars('z_01_ALL_nDisrupt-0_evalSimStats_meanCI_HP25.csv', 'interval', 'ratioPathsFound', 'Campaign', 'z_01_nDisrupt-0_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-1_evalSimStats_meanCI_HP25.csv', 'interval', 'ratioPathsFound', 'Campaign', 'z_01_nDisrupt-1_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-2_evalSimStats_meanCI.csv', 'interval', 'ratioPathsFound', 'Campaign', 'z_01_nDisrupt-2')

# plot_with_error_bars('z_01_ALL_nDisrupt-0_evalSimStats_meanCI_HP25.csv', 'interval', 'avgFlowThroughput', 'Campaign', 'z_01_nDisrupt-0_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-1_evalSimStats_meanCI_HP25.csv', 'interval', 'avgFlowThroughput', 'Campaign', 'z_01_nDisrupt-1_HP25')

# plot_with_error_bars('z_01_ALL_nDisrupt-0_evalSimStats_meanCI_HP25.csv', 'interval', 'avgNDisruptions', 'Campaign', 'z_01_nDisrupt-0_HP25')
# plot_with_error_bars('z_01_ALL_nDisrupt-1_evalSimStats_meanCI_HP25.csv', 'interval', 'avgNDisruptions', 'Campaign', 'z_01_nDisrupt-1_HP25')


plot_with_error_bars('z_02_ALL_evalSimStats_meanCI.csv', 'd', 'avgFlowLossRatio', 'Campaign', 'z_02')
plot_with_error_bars('z_02_ALL_evalSimStats_meanCI.csv', 'd', 'avgFlowMeanDelay', 'Campaign', 'z_02')
plot_with_error_bars('z_02_ALL_evalSimStats_meanCI.csv', 'd', 'sysNRouteChanges', 'Campaign', 'z_02')
plot_with_error_bars('z_02_ALL_evalSimStats_meanCI.csv', 'd', 'avgMaxNHops', 'Campaign', 'z_02')
plot_with_error_bars('z_02_ALL_evalSimStats_meanCI.csv', 'd', 'sysRoutingOhBytes', 'Campaign', 'z_02')
plot_with_error_bars('z_02_ALL_evalSimStats_meanCI.csv', 'd', 'ratioPathsFound', 'Campaign', 'z_02')


end_time = time.time()
duration_s = end_time - start_time
duration_m = duration_s / 60.0
print(f"The script took {duration_m} minutes to run.")
