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
    fig, axs = plt.subplots(2, 6, figsize=(20, 10))

    for i, metric in enumerate(metrics):
        # Calculate position in the grid
        row_pos = i // 6
        col_pos = i % 6

        # Bar plot for each metric with different color
        axs[row_pos, col_pos].bar(rng_run, data[metric], color=colors[i % len(colors)])
        axs[row_pos, col_pos].set_title(metric)
        axs[row_pos, col_pos].set_xlabel('RngRun')
        axs[row_pos, col_pos].set_ylabel(metric)
        axs[row_pos, col_pos].set_xticks(rng_run)

    plt.tight_layout()
    plot_path = os.path.join(output_path, "evalSimStats_allRuns.png")
    plt.savefig(plot_path)
    plt.close()  # Close the plot


def plot_metrics_CI(evalOutputPath, metrics, colors):
    
    # Reading the merged data
    filepath = os.path.join(evalOutputPath, 'evalSimStats.csv')
    data = np.genfromtxt(filepath, delimiter=',', names=True, dtype=None, encoding='utf-8')

    # Prepare the figure for plotting
    fig, axs = plt.subplots(2, 6, figsize=(20, 10))

    for i, metric in enumerate(metrics):
        metric_data = data[metric]  # Get data for each metric using the header name
        
        # Filter out NaN values
        metric_data = metric_data[~np.isnan(metric_data)]

        # Check if there's enough data left
        if len(metric_data) > 1:
            mean = np.mean(metric_data)
            se = np.std(metric_data, ddof=1) / np.sqrt(len(metric_data))
            ci95 = 1.96 * se

            row_pos = i // 6
            col_pos = i % 6

            # Plotting the mean with the 95% CI as error bars
            axs[row_pos, col_pos].bar(1, mean, color=colors[i % len(colors)])
            axs[row_pos, col_pos].errorbar(1, mean, yerr=ci95, fmt='o', color='black')
            axs[row_pos, col_pos].set_title(metric)
            axs[row_pos, col_pos].set_xticks([])
            axs[row_pos, col_pos].set_ylim([0, mean + ci95 * 3])  # Adjust y-limit to fit the CI
        else:
            # Handle the case where there's not enough data
            axs[i // 6, i % 6].text(0.5, 0.5, 'Insufficient data for\n' + metric, ha='center')

    plt.tight_layout()
    plt.savefig(os.path.join(evalOutputPath, "evalSimStats_meanCI.png"))
    plt.close()  # Close the plot


def plot_combined_metrics_perRun(metrics, all_eval_data, eval_names, output_file_name, campaignLabel):
    fig, axs = plt.subplots(2, 6, figsize=(20, 10))
    
    colors = plt.cm.tab10.colors

    for i, metric in enumerate(metrics):
        row_pos = i // 6
        col_pos = i % 6

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
        axs[row_pos, col_pos].set_xlabel('RngRun')
        axs[row_pos, col_pos].set_xticklabels(rng_runs)
        axs[row_pos, col_pos].legend(eval_names)

    plt.tight_layout()
    png_file_name = f"{output_file_name}.png"
    plt.savefig(png_file_name)
    plt.close()  # Close the plot



def plot_combined_metrics_meanCI(metrics, all_eval_data, eval_names, output_file_name, campaignLabel, campaignName):
    fig, axs = plt.subplots(2, 6, figsize=(20, 10))

    colors = plt.cm.tab10.colors

    # Initialize a dictionary to hold the statistics for each campaignLabel
    stats_dict = {metric: [] for metric in metrics}

    for i, metric in enumerate(metrics):
        row_pos = i // 6
        col_pos = i % 6
        max_mean = 0  # Initialize the maximum mean value for the current metric

        for j, eval_data in enumerate(all_eval_data):
            
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



def main(evalData, campaignName, campaignLabel):

    metrics = ["ratioPathsFound","sysLossRatio", "avgMinNHops", "avgMaxNHops",
            "avgFlowLossRatio", "avgFlowMeanDelay", "avgFlowThroughput",
             "sysNRouteChanges", "sysRoutingOhPkt", "sysRoutingOhBytes", "avgNDisruptions" ]

    colors = plt.cm.tab10.colors
    all_eval_data = []
    eval_names = []

    for evalPath, evalName in evalData.items():
        print(f"Processing {evalName} at {evalPath}")

        # Merge CSV files and gather data for plotting simStats
        merged_rows, header = merge_csv_files_simStats(evalPath)

        # Plot metrics per Run
        plot_metrics_per_run(merged_rows, metrics, colors, evalPath, header)

        # Plot simulation statistics (mean and 95% of the metrics)
        plot_metrics_CI(evalPath, metrics, colors)

        # Store data for combined plots
        eval_data = {metric: [float(row[header.index(metric)]) for row in merged_rows[1:]] for metric in metrics}
        eval_data['RngRun'] = [int(row[header.index('RngRun')]) for row in merged_rows[1:]]
        all_eval_data.append(eval_data)
        eval_names.append(evalName)


    # Generate combined plots
    combined_metrics_filename = f"{campaignName}_evalSimStats_perRun"
    combined_statistics_filename = f"{campaignName}_evalSimStats_meanCI"

    plot_combined_metrics_perRun(metrics, all_eval_data, eval_names, combined_metrics_filename, campaignLabel)
    plot_combined_metrics_meanCI(metrics, all_eval_data, eval_names, combined_statistics_filename, campaignLabel, campaignName)





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