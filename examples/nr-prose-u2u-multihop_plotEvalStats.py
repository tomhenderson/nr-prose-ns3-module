import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import argparse

def merge_csv_files_simStats(evalOutputPath):
    all_rows = []
    header_saved = False

    for run_dir in os.listdir(evalOutputPath):
        run_path = os.path.join(evalOutputPath, run_dir)
        if os.path.isdir(run_path):
            csv_file = os.path.join(run_path, "simStats.csv")
            if os.path.exists(csv_file):
                with open(csv_file, 'r') as file:
                    csv_reader = csv.reader(file)
                    header = next(csv_reader)
                    if not header_saved:
                        all_rows.append(header)
                        header_saved = True
                    all_rows.extend(csv_reader)

    # Save the merged rows to a CSV file in the evalOutputPath
    csv_file_path = os.path.join(evalOutputPath, 'evalSimStats.csv')
    with open(csv_file_path, 'w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerows(all_rows)

    return all_rows, header

def plot_metrics_per_run(rows, metrics, colors, output_path, header):
    data = {metric: [] for metric in metrics}
    rng_run = []

    for row in rows[1:]:  # Skip header
        if int(row[1]) not in rng_run:
            rng_run.append(int(row[1]))  # RngRun is the second column
        for metric in metrics:
            data[metric].append(float(row[header.index(metric)]))

    # Creating subplots
    fig, axs = plt.subplots(2, 5, figsize=(20, 10))

    for i, metric in enumerate(metrics):
        # Calculate position in the grid
        row_pos = i // 5
        col_pos = i % 5

        # Bar plot for each metric with different color
        axs[row_pos, col_pos].bar(rng_run, data[metric], color=colors[i % len(colors)])
        axs[row_pos, col_pos].set_title(metric)
        axs[row_pos, col_pos].set_xlabel('RngRun')
        axs[row_pos, col_pos].set_ylabel(metric)
        axs[row_pos, col_pos].set_xticks(rng_run)

    plt.tight_layout()
    plot_path = os.path.join(output_path, "evalSimStats.png")
    plt.savefig(plot_path)
    plt.close()  # Close the plot


def plot_metrics_CI(evalOutputPath, metrics, colors):
    
    # Reading the merged data
    filepath = os.path.join(evalOutputPath, 'evalSimStats.csv')
    data = np.genfromtxt(filepath, delimiter=',', names=True, dtype=None, encoding='utf-8')

    # Prepare the figure for plotting
    fig, axs = plt.subplots(2, 5, figsize=(20, 10))

    for i, metric in enumerate(metrics):
        metric_data = data[metric]  # Get data for each metric using the header name
        
        # Filter out NaN values
        metric_data = metric_data[~np.isnan(metric_data)]

        # Check if there's enough data left
        if len(metric_data) > 1:
            mean = np.mean(metric_data)
            se = np.std(metric_data, ddof=1) / np.sqrt(len(metric_data))
            ci95 = 1.96 * se

            row_pos = i // 5
            col_pos = i % 5

            # Plotting the mean with the 95% CI as error bars
            axs[row_pos, col_pos].bar(1, mean, color=colors[i % len(colors)])
            axs[row_pos, col_pos].errorbar(1, mean, yerr=ci95, fmt='o', color='black')
            axs[row_pos, col_pos].set_title(metric)
            axs[row_pos, col_pos].set_xticks([])
            axs[row_pos, col_pos].set_ylim([0, mean + ci95 * 3])  # Adjust y-limit to fit the CI
        else:
            # Handle the case where there's not enough data
            axs[i // 5, i % 5].text(0.5, 0.5, 'Insufficient data for\n' + metric, ha='center')

    plt.tight_layout()
    plt.savefig(os.path.join(evalOutputPath, "evalSimStats_meanCI.png"))
    plt.close()  # Close the plot


def gather_flows_stats(evalOutputPath):
    flows_data = []
    for run_dir in os.listdir(evalOutputPath):
        run_path = os.path.join(evalOutputPath, run_dir)
        if os.path.isdir(run_path):
            csv_file = os.path.join(run_path, "flowsStats.csv")
            if os.path.exists(csv_file):
                with open(csv_file, 'r') as file:
                    csv_reader = csv.DictReader(file)
                    for row in csv_reader:
                        flows_data.append(row)
    return flows_data

def calculate_flow_statistics(flows_data, eval_name):
    # Organize data by nHops and calculate statistics
    metrics_by_hops = {}
    for row in flows_data:
        nHops = int(row['nHops'])
        if nHops not in metrics_by_hops:
            metrics_by_hops[nHops] = {'LossRatio': [], 'MeanDelay': [], 'MeanJitter': []}

        metrics_by_hops[nHops]['LossRatio'].append(float(row['LossRatio']))
        metrics_by_hops[nHops]['MeanDelay'].append(float(row['MeanDelay']))
        metrics_by_hops[nHops]['MeanJitter'].append(float(row['MeanJitter']))

    stats_by_hops = {}
    for nHops, metrics in metrics_by_hops.items():
        stats_by_hops[nHops] = {}
        for metric, values in metrics.items():
            mean = np.mean(values)
            se = np.std(values, ddof=1) / np.sqrt(len(values))
            ci95 = 1.96 * se
            stats_by_hops[nHops][metric] = (mean, ci95)

    return stats_by_hops, eval_name

def plot_flow_statistics(all_stats_by_hops, eval_names, output_path, output_file_name):
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    metrics = ['LossRatio', 'MeanDelay', 'MeanJitter']
    colors = plt.cm.tab10.colors

    for i, metric in enumerate(metrics):
        # Collect all unique nHops values across all evaluations
        nHops_sorted = sorted(set(hops for stats_by_hops in all_stats_by_hops for hops in stats_by_hops.keys()))

        for j, nHops in enumerate(nHops_sorted):
            # For each nHops, plot a bar for each evaluation
            for k, stats_by_hops in enumerate(all_stats_by_hops):
                mean, ci95 = stats_by_hops.get(nHops, {}).get(metric, (0, 0))
                bar_position = j - 0.15 + k * 0.1
                axs[i].bar(bar_position, mean, width=0.1, yerr=ci95, color=colors[k % len(colors)], label=eval_names[k] if j == 0 else "")

        axs[i].set_title(metric)
        axs[i].set_xlabel('nHops')
        axs[i].set_xticks(range(len(nHops_sorted)))
        axs[i].set_xticklabels(nHops_sorted)
        if i == 0:
            axs[i].legend()

    plt.tight_layout()
    png_file_name = f"{output_file_name}.png"

    plt.savefig(os.path.join(output_path, png_file_name))
    plt.close()  # Close the plot


def plot_combined_metrics(metrics, all_eval_data, eval_names, output_file_name, campaignLabel):
    fig, axs = plt.subplots(2, 5, figsize=(20, 10))
    
    colors = plt.cm.tab10.colors

    for i, metric in enumerate(metrics):
        row_pos = i // 5
        col_pos = i % 5

        # Find unique RngRuns across all evaluations
        rng_runs = sorted(set(run for eval_data in all_eval_data for run in eval_data['RngRun']))

        for j, run in enumerate(rng_runs):
            # For each RngRun, plot a bar for each evaluation
            for k, eval_data in enumerate(all_eval_data):
                data = eval_data[metric]
                run_index = eval_data['RngRun'].index(run) if run in eval_data['RngRun'] else -1
                value = data[run_index] if run_index != -1 else 0

                # Calculate position of the bar
                bar_position = j + k * 0.1 - (len(all_eval_data) * 0.1) / 2
                axs[row_pos, col_pos].bar(bar_position, value, width=0.1, color=colors[k % len(colors)])

        axs[row_pos, col_pos].set_title(metric+" by "+campaignLabel)
        axs[row_pos, col_pos].set_xticks(range(len(rng_runs)))
        axs[row_pos, col_pos].set_xticklabels(rng_runs)
        axs[row_pos, col_pos].legend(eval_names)


    
    plt.tight_layout()
    png_file_name = f"{output_file_name}.png"
    plt.savefig(png_file_name)
    plt.close()  # Close the plot



def calculate_and_plot_combined_statistics(metrics, all_eval_data, eval_names, output_file_name, campaignLabel, campaignName):
    fig, axs = plt.subplots(2, 5, figsize=(20, 10))

    colors = plt.cm.tab10.colors

    # Initialize a dictionary to hold the statistics for each campaignLabel
    stats_dict = {metric: [] for metric in metrics}

    for i, metric in enumerate(metrics):
        row_pos = i // 5
        col_pos = i % 5
        max_mean = 0  # Initialize the maximum mean value for the current metric

        for j, eval_data in enumerate(all_eval_data):
            data = eval_data[metric]
            data = np.array(eval_data[metric])  

            # Filter out NaN values
            filtered_data = data[~np.isnan(data)]

            # Check if there's enough data left
            if len(filtered_data) > 1:
                mean = np.mean(filtered_data)
                se = np.std(filtered_data, ddof=1) / np.sqrt(len(filtered_data))
                ci95 = 1.96 * se

                if mean > max_mean:  # Update max_mean if current mean is higher
                    max_mean = mean
            else:
                #Handle the case with insufficient data
                mean=0
                ci95=0
            
             # Store the statistics in the dictionary
            stats_dict[metric].append((mean, ci95))

            axs[row_pos, col_pos].bar(j, mean, color=colors[j % len(colors)], label=eval_names[j])
            axs[row_pos, col_pos].errorbar(j, mean, yerr=ci95, fmt='o', color='black')
            axs[row_pos, col_pos].text(j-0.25, mean, f'{mean:.2f}', ha='center', va='bottom', rotation=90)
 
        if metric == "ratioPathsFound" or metric ==  "sysLossRatio" or metric == "avgFlowLossRatio":
            axs[row_pos, col_pos].set_ylim([0, 1.0])  
        else:
            # Increase ylim to 20% above the highest mean for legibility
            axs[row_pos, col_pos].set_ylim([0, max_mean * 1.2])

        axs[row_pos, col_pos].set_title(metric)
        axs[row_pos, col_pos].set_xticks(range(len(eval_names)))
        axs[row_pos, col_pos].set_xticklabels(eval_names)
        axs[row_pos, col_pos].set_xlabel(campaignLabel)


    plt.tight_layout()
    png_file_name = f"{output_file_name}.png"
    plt.savefig(png_file_name)
    plt.close()  # Close the plot

    # Write the statistics to a CSV file
    csv_file_name = f"{output_file_name}.csv"
    with open(csv_file_name, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Campaign", campaignLabel, *metrics])

        for i in range(len(eval_names)):
            row = [campaignName, eval_names[i]]
            for metric in metrics:
                mean_ci = stats_dict[metric][i]
                row.append(f"({mean_ci[0]:.6f}, {mean_ci[1]:.6f})")
            writer.writerow(row)


def merge_csv_files_simPhyStats(evalOutputPath):
    all_rows = []
    header_saved = False

    for run_dir in os.listdir(evalOutputPath):
        run_path = os.path.join(evalOutputPath, run_dir)
        if os.path.isdir(run_path):
            csv_file = os.path.join(run_path, "simPhyStats.csv")
            if os.path.exists(csv_file):
                with open(csv_file, 'r') as file:
                    csv_reader = csv.reader(file)
                    header = next(csv_reader)
                    if not header_saved:
                        all_rows.append(header)
                        header_saved = True
                    all_rows.extend(csv_reader)

    # Save the merged rows to a CSV file in the evalOutputPath
    csv_file_path = os.path.join(evalOutputPath, 'evalSimPhyStats.csv')
    with open(csv_file_path, 'w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerows(all_rows)

    return all_rows, header


def plot_phy_stats_metrics_per_run(rows, metrics, colors, output_path, header):
    data = {metric: [] for metric in metrics}
    rng_run = []

    for row in rows[1:]:  # Skip header
        if int(row[1]) not in rng_run:
            rng_run.append(int(row[1]))  # RngRun is the second column
        for metric in metrics:
            data[metric].append(float(row[header.index(metric)]))

    # Creating subplots
    fig, axs = plt.subplots(2, 5, figsize=(20, 10))

    for i, metric in enumerate(metrics):
        # Calculate position in the grid
        row_pos = i // 5
        col_pos = i % 5

        # Bar plot for each metric with different color
        axs[row_pos, col_pos].bar(rng_run, data[metric], color=colors[i % len(colors)])
        axs[row_pos, col_pos].set_title(metric)
        axs[row_pos, col_pos].set_xlabel('RngRun')
        axs[row_pos, col_pos].set_ylabel(metric)
        axs[row_pos, col_pos].set_xticks(rng_run)

    plt.tight_layout()
    plot_path = os.path.join(output_path, "evalSimPhyStatsPercents.png")
    plt.savefig(plot_path)
    plt.close()  # Close the plot


def plot_phy_stats_metrics_CI(evalOutputPath, metrics, colors):
    
    # Reading the merged data
    filepath = os.path.join(evalOutputPath, 'evalSimPhyStatsPercents.csv')
    data = np.genfromtxt(filepath, delimiter=',', names=True, dtype=None, encoding='utf-8')

    # Prepare the figure for plotting
    fig, axs = plt.subplots(2, 5, figsize=(20, 10))

    for i, metric in enumerate(metrics):
        metric_data = data[metric]  # Get data for each metric using the header name
        
        # Filter out NaN values
        metric_data = metric_data[~np.isnan(metric_data)]

        # Check if there's enough data left
        if len(metric_data) > 1:
            mean = np.mean(metric_data)
            se = np.std(metric_data, ddof=1) / np.sqrt(len(metric_data))
            ci95 = 1.96 * se

            row_pos = i // 5
            col_pos = i % 5

            # Plotting the mean with the 95% CI as error bars
            axs[row_pos, col_pos].bar(1, mean, color=colors[i % len(colors)])
            axs[row_pos, col_pos].errorbar(1, mean, yerr=ci95, fmt='o', color='black')
            axs[row_pos, col_pos].set_title(metric)
            axs[row_pos, col_pos].set_xticks([])
            axs[row_pos, col_pos].set_ylim([0, mean + ci95 * 3])  # Adjust y-limit to fit the CI
        else:
            # Handle the case where there's not enough data
            axs[i // 5, i % 5].text(0.5, 0.5, 'Insufficient data for\n' + metric, ha='center')

    plt.tight_layout()
    plt.savefig(os.path.join(evalOutputPath, "evalSimPhyStatsPercents_meanCI.png"))
    plt.close()  # Close the plot



# Function to calculate phy stats percentages
def calculate_percentages(row):
    # Convert string values to integers
    for key in row:
        row[key] = int(row[key])
    
    # Calculations
    results = {
        "RngSeed": row["RngSeed"],
        "RngRun": row["RngRun"],
        "percentSuccessRxCtrl": (row["totalnRxCtrl"] / row["totalnTxCtrl"] * 100) if row["totalnTxCtrl"] > 0 else 0,
        "percentSuccessRxData": (row["totalnRxData"] / row["totalnTxData"] * 100) if row["totalnTxData"] > 0 else 0,
        "percentCorruptRxCtrl": (row["totalnCorruptRxCtrl"] / row["totalnTxCtrl"] * 100) if row["totalnTxCtrl"] > 0 else 0,
        "percentCorruptRxData": (row["totalnCorruptRxData"] / row["totalnTxData"] * 100) if row["totalnTxData"] > 0 else 0,
        "percentHdRxCtrl": (row["totalnHdRxCtrl"] / row["totalnTxCtrl"] * 100) if row["totalnTxCtrl"] > 0 else 0,
        "percentHdRxData": (row["totalnHdRxData"] / row["totalnTxData"] * 100) if row["totalnTxData"] > 0 else 0,
        "percentIgnoredData": (row["totalnIgnoredData"] / row["totalnTxData"] * 100) if row["totalnTxData"] > 0 else 0,
        "percentUnaccountedCtrl": 100 - ((row["totalnRxCtrl"] + row["totalnCorruptRxCtrl"] + row["totalnHdRxCtrl"]) / row["totalnTxCtrl"] * 100) if row["totalnTxCtrl"] > 0 else 0,
        "percentUnaccountedData": 100 - ((row["totalnRxData"] + row["totalnCorruptRxData"] + row["totalnHdRxData"] + row["totalnIgnoredData"]) / row["totalnTxData"] * 100) if row["totalnTxData"] > 0 else 0,
    }

    return results

def process_csv_files_simPhyStats (evalOutputPath):

    input_filename = os.path.join(evalOutputPath, 'evalSimPhyStats.csv')
    output_filename = os.path.join(evalOutputPath, 'evalSimPhyStatsPercents.csv')
    all_rows = []  # List to store all processed rows for return

    # Open the input CSV file
    with open(input_filename, mode='r', newline='') as infile:
        reader = csv.DictReader(infile)

        # Prepare the output CSV file
        fieldnames = ['RngSeed', 'RngRun', 'percentSuccessRxCtrl', 'percentSuccessRxData', 'percentCorruptRxCtrl', 'percentCorruptRxData', 'percentHdRxCtrl', 'percentHdRxData', 'percentIgnoredData', 'percentUnaccountedCtrl', 'percentUnaccountedData']
        with open(output_filename, mode='w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()

            # Process each row from the input file and write to the output file
            for row in reader:
                if int(row["totalnTxCtrl"]) > 0 and int(row["totalnTxData"]) > 0:
                    percentages = calculate_percentages(row)
                    # Check if success rates for control or data are less than 100 #TODO WHY? REmove?
 #                   if percentages["percentSuccessRxCtrl"] < 100 or percentages["percentSuccessRxData"] < 100:
                    writer.writerow(percentages)
                    # Convert the processed dictionary back to a CSV row format
                    processed_row = [percentages[field] for field in fieldnames]
                    all_rows.append(processed_row)  # Add the processed row to the list

    return all_rows, fieldnames  # Return all processed rows and the header



def merge_csv_files_simPhyStats_global_data(evalOutputPath):
    all_rows = []
    header_saved = False

    for run_dir in os.listdir(evalOutputPath):
        run_path = os.path.join(evalOutputPath, run_dir)
        if os.path.isdir(run_path):
            csv_file = os.path.join(run_path, "simPhyStats_PerNode_Data_total.csv")
            if os.path.exists(csv_file):
                with open(csv_file, 'r') as file:
                    csv_reader = csv.reader(file)
                    header = next(csv_reader)
                    if not header_saved:
                        all_rows.append(header)
                        header_saved = True
                    all_rows.extend(csv_reader)

    # Save the merged rows to a CSV file in the evalOutputPath
    csv_file_path = os.path.join(evalOutputPath, 'evalSimPhyStats_global.csv')
    with open(csv_file_path, 'w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerows(all_rows)

    return all_rows, header


def plot_phy_stats_global_metrics_CI(evalOutputPath, metrics, colors):
    
    # Reading the merged data
    filepath = os.path.join(evalOutputPath, 'evalSimPhyStats_global.csv')
    data = np.genfromtxt(filepath, delimiter=',', names=True, dtype=None, encoding='utf-8')

    # Prepare the figure for plotting
    fig, axs = plt.subplots(2, 5, figsize=(20, 10))

    for i, metric in enumerate(metrics):
        metric_data = data[metric]  # Get data for each metric using the header name
        
        # Filter out NaN values
        metric_data = metric_data[~np.isnan(metric_data)]

        # Check if there's enough data left
        if len(metric_data) > 1:
            mean = np.mean(metric_data)
            se = np.std(metric_data, ddof=1) / np.sqrt(len(metric_data))
            ci95 = 1.96 * se

            row_pos = i // 5
            col_pos = i % 5

            # Plotting the mean with the 95% CI as error bars
            axs[row_pos, col_pos].bar(1, mean, color=colors[i % len(colors)])
            axs[row_pos, col_pos].errorbar(1, mean, yerr=ci95, fmt='o', color='black')
            axs[row_pos, col_pos].set_title(metric)
            axs[row_pos, col_pos].set_xticks([])
            axs[row_pos, col_pos].set_ylim([0, mean + ci95 * 3])  # Adjust y-limit to fit the CI
        else:
            # Handle the case where there's not enough data
            axs[i // 5, i % 5].text(0.5, 0.5, 'Insufficient data for\n' + metric, ha='center')

    plt.tight_layout()
    plt.savefig(os.path.join(evalOutputPath, "evalSimPhyStatsGlobal_meanCI.png"))
    plt.close()  # Close the plot


def main(evalData, campaignName, campaignLabel):

    metrics = ["ratioPathsFound","sysLossRatio", "avgNHops", "minRelayPathCount", "maxRelayPathCount",
            "avgFlowLossRatio", 
            "avgFlowMeanDelay", "avgFlowMeanJitter", "minNonZeroRelayPathCount", 
            "meanNonZeroRelayPathCount" ]

    phy_metrics = ["percentSuccessRxData", "percentCorruptRxData", "percentHdRxData", "percentIgnoredData", "percentUnaccountedData",
                    "percentSuccessRxCtrl", "percentCorruptRxCtrl","percentHdRxCtrl", "percentUnaccountedCtrl"]

    phy_global_metrics = ["nTxData","nRxData","nCorruptRxData","nHdRxData","nIgnoredData"]


    colors = plt.cm.tab10.colors
    all_eval_data = []
    phy_all_eval_data = []
    phy_global_all_eval_data = []
    eval_names = []
    all_stats_by_hops = []

    for evalPath, evalName in evalData.items():
        print(f"Processing {evalName} at {evalPath}")

        # Merge CSV files and gather data for plotting simStats
        merged_rows, header = merge_csv_files_simStats(evalPath)

        # Plot metrics per Run
        plot_metrics_per_run(merged_rows, metrics, colors, evalPath, header)

        # Plot simulation statistics (mean and 95% of the metrics)
        plot_metrics_CI(evalPath, metrics, colors)

        # Gather and calculate flow statistics
        flows_data = gather_flows_stats(evalPath)
        stats_by_hops, _ = calculate_flow_statistics(flows_data, evalName)

        # Plot individual flow statistics
        plot_flow_statistics([stats_by_hops], [evalName], evalPath, "flowStats.png")

        # Store data for combined plots
        eval_data = {metric: [float(row[header.index(metric)]) for row in merged_rows[1:]] for metric in metrics}
        eval_data['RngRun'] = [int(row[header.index('RngRun')]) for row in merged_rows[1:]]
        all_eval_data.append(eval_data)
        eval_names.append(evalName)
        all_stats_by_hops.append(stats_by_hops)

        # Merge CSV files and gather data for plotting simPhyStats
        merge_csv_files_simPhyStats(evalPath)
        phy_merged_rows, phy_header = process_csv_files_simPhyStats (evalPath)

        plot_phy_stats_metrics_per_run(phy_merged_rows, phy_metrics, colors, evalPath, phy_header)
        plot_phy_stats_metrics_CI(evalPath, phy_metrics, colors)
        phy_eval_data = {phy_metric: [float(row[phy_header.index(phy_metric)]) for row in phy_merged_rows[1:]] for phy_metric in phy_metrics}
        phy_eval_data['RngRun'] = [int(row[phy_header.index('RngRun')]) for row in phy_merged_rows[1:]]
        phy_all_eval_data.append(phy_eval_data)

        # Merge CSV files and gather data for plotting simStats globals (total)
        phy_global_merged_rows, phy_global_header = merge_csv_files_simPhyStats_global_data(evalPath)
        plot_phy_stats_global_metrics_CI(evalPath, phy_global_metrics, colors)
        phy_global_eval_data = {phy_global_metric: [float(row[phy_global_header.index(phy_global_metric)]) for row in phy_global_merged_rows[1:]] for phy_global_metric in phy_global_metrics}
        phy_global_eval_data['RngRun'] = [int(row[phy_global_header.index('RngRun')]) for row in phy_global_merged_rows[1:]]
        phy_global_all_eval_data.append(phy_global_eval_data)


    # Generate combined plots
    combined_metrics_filename = f"{campaignName}_evalSimStats"
    combined_statistics_filename = f"{campaignName}_evalSimStats_meanCI"
    combined_flow_stats_filename = f"{campaignName}_flowStats"

    plot_combined_metrics(metrics, all_eval_data, eval_names, combined_metrics_filename, campaignLabel)
    calculate_and_plot_combined_statistics(metrics, all_eval_data, eval_names, combined_statistics_filename, campaignLabel, campaignName)
    plot_flow_statistics(all_stats_by_hops, eval_names, "", combined_flow_stats_filename)

    phy_combined_statistics_filename = f"{campaignName}_evalSimPhyStats_meanCI"
    calculate_and_plot_combined_statistics(phy_metrics, phy_all_eval_data, eval_names, phy_combined_statistics_filename, campaignLabel, campaignName)
    
    phy_global_combined_statistics_filename = f"{campaignName}_evalSimPhyStatsGlobal_meanCI"
    calculate_and_plot_combined_statistics(phy_global_metrics, phy_global_all_eval_data, eval_names, phy_global_combined_statistics_filename, campaignLabel, campaignName)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse evaluation and generate plots.")
    parser.add_argument('--eval', action='append', nargs=2, metavar=('PATH', 'NAME'),
                        help="Add evaluation with PATH and NAME", required=True)
    parser.add_argument('--campaignName', type=str, required=True,
                        help="The name of the campaign to prefix the output files")
    parser.add_argument('--campaignLabel', type=str, required=True,
                        help="The label to show in the X axis of the campaign plots")
    args = parser.parse_args()

    evalData = {path: name for path, name in args.eval}
    main(evalData, args.campaignName, args.campaignLabel)