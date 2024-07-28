
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import csv
from collections import defaultdict
import argparse
from math import cos, sin, radians
import numpy as np

def plot_nodes_and_paths(filename, simOutputPath):
    
    # Initialize lists to store data
    node_x = []
    node_y = []
    node_labels = []
    node_colors = []
    path_data = {}
    path_colors = ['blue', 'orange', 'green', 'red', 'purple']  # Add more colors if needed
    d = 0

    # Read the data from the CSV file
    csv_file = f"{simOutputPath}/{filename}"
    with open(csv_file, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            if row[0].isdigit():  # Node data rows start with a digit (node ID)
                node_x.append(float(row[1]))
                node_y.append(float(row[2]))
                node_labels.append(row[0])
                color = 'blue' if row[3] == 'relay' else 'black'
                node_colors.append(color)
            elif row[0] == 'Path':  # Path data rows start with the word 'Path'
                path_id = int(row[1])
                node_ids = row[2:]
                path_data[path_id] = node_ids
            elif row[0] == 'D':  # Deployment side
                d = float(row[1])

    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=(8, 8))  # Setting both width and height to the same value

    # Plot nodes
    ax.scatter(node_x, node_y, c=node_colors)
    for i, label in enumerate(node_labels):
        ax.text(node_x[i] + 1, node_y[i] + 1, label, fontsize=9)  # Offset the node labels for clarity

    # Create custom legend entries for node types
    legend_elements = [
        mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=10, label='Non-relay'),
        mlines.Line2D([], [], color='blue', marker='o', linestyle='None', markersize=10, label='Relay')
    ]

    # Define line width range for paths, and set transparency
    min_line_width = 1
    max_line_width = 5
    line_alpha = 0.3  # Transparency level

    # Calculate line width increment if there are multiple paths
    line_width_increment = (max_line_width - min_line_width) / (len(path_data) - 1) if len(path_data) > 1 else 0

    # Plot paths and add them to the legend entries
    for i, (path_id, nodes) in enumerate(sorted(path_data.items())):
        path_x = [node_x[node_labels.index(node_id)] for node_id in nodes if node_id in node_labels]
        path_y = [node_y[node_labels.index(node_id)] for node_id in nodes if node_id in node_labels]
        if path_x and path_y:
            line_width = max_line_width - i * line_width_increment
            line_color = path_colors[path_id % len(path_colors)]  # Cycle through path_colors
            path_label = ' -> '.join(nodes)  # Creates a string representing the actual path
            line, = ax.plot(path_x, path_y, label=f'Path {path_id}: {path_label}', color=line_color, 
                            linewidth=line_width, alpha=line_alpha)
            legend_elements.append(line)

    # Add custom legend including both node types and paths
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')

    # Set the x and y limits to the value of 'D'
    ax.set_xlim(0, d)
    ax.set_ylim(0, d)

    # Set aspect of the plot to be equal
    ax.set_aspect('equal', adjustable='box')

    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')

    plt.tight_layout()  # Adjust the padding between and around subplots

    # Save the plot as an image file
    plot_file_path = f"{simOutputPath}/deployment.png"
    plt.savefig(plot_file_path, bbox_inches='tight', dpi=300)
    plt.close()  # Close the plot


def plot_flow_statistics(filename, simOutputPath):
    # Initialize lists to store data
    flow_ids = []
    loss_ratios = []
    mean_delays = []
    mean_jitters = []
    nTxPkts = []
    nRxPkts = []
    nHops = []
    
    # Read the CSV file
    csv_file = f"{simOutputPath}/{filename}"
    with open(csv_file, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip the header row
        for row in csv_reader:
            flow_ids.append(int(row[0]))
            loss_ratios.append(float(row[8]))
            mean_delays.append(float(row[9]))
            mean_jitters.append(float(row[10]))
            nTxPkts.append(int(row[6]))
            nRxPkts.append(int(row[7]))
            nHops.append(int(row[5]))

    # Calculate Packet Delivery Ratio
    packet_delivery_ratio = [1 - lr for lr in loss_ratios]
    # Define colors for each plot
    plot_colors = ['skyblue', 'tomato', 'deepskyblue', 'mediumseagreen', 'gold', 'orchid']

    # Create a figure for the plots
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(16, 10))  # Adjusted for 2 rows, 3 columns
    fig.suptitle('Path/Flow Statistics')

     # Plot Packet Stats, Packet Delivery Ratio, Packet Loss Ratio, Number of Hops, Mean Delay, Mean Jitter in the specified order
     # Plot Packet Stats, Packet Delivery Ratio, Packet Loss Ratio, Number of Hops, Mean Delay, Mean Jitter in the specified order
    plot_order = [
        (nTxPkts, nRxPkts, 'Packet Stats', axs[0, 0]),
        (packet_delivery_ratio, None, 'Packet Delivery Ratio', axs[0, 1]),
        (loss_ratios, None, 'Packet Loss Ratio', axs[0, 2]),
        (nHops, None, 'Number of Hops', axs[1, 0]),
        (mean_delays, None, 'Mean Delay (ms)', axs[1, 1]),
        (mean_jitters, None, 'Mean Jitter (ms)', axs[1, 2])
    ]

    for i, (data1, data2, title, ax) in enumerate(plot_order):
        if title == 'Packet Stats':
            # Special handling for Packet Stats
            width = 0.35
            tx_positions = [j - width/2 for j in range(len(flow_ids))]
            rx_positions = [j + width/2 for j in range(len(flow_ids))]
            ax.bar(tx_positions, data1, width, label='nTxPkts')
            ax.bar(rx_positions, data2, width, label='nRxPkts')
            # Dynamically adjust y-axis limit
            max_height = max(max(data2), max(data1))
            ax.set_ylim(0, max_height * 1.15)  # Add extra space at the top
            ax.legend(loc='lower right')
            ax.set_xlabel('Flow Id')

            # Annotate each bar in Packet Stats
            for j, (x, y) in enumerate(zip(tx_positions, data1)):
                ax.text(x, y*1.01, f'{data1[j]}', ha='center', va='bottom', rotation=90, fontsize=14)
            for j, (x, y) in enumerate(zip(rx_positions, data2)):
                ax.text(x, y*1.01, f'{data2[j]}', ha='center', va='bottom', rotation=90, fontsize=14)


        else:
            bars = ax.bar(flow_ids, data1, color=plot_colors[i])
            # Annotate each bar for other plots
            for bar in bars:
                yval = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2, yval*1.01, f'{yval:.2f}', ha='center', va='bottom', rotation=90, fontsize=14)
            ax.set_ylim(0, max_height * 1.15)  # Add extra space at the top
            ax.set_title(title, fontsize=14)
            ax.set_xlabel('Flow Id')
            ax.set_ylabel(title, fontsize=14)
            if 'Ratio' in title:  # Set Y-axis limit to 1.0 for ratio plots
                ax.set_ylim(0, 1.1)
            else:
                max_height = max(data1)
                ax.set_ylim(0, max_height * 1.15)  # Add extra space at the top

        ax.set_xticks(range(len(flow_ids)))

    # Adjust the bottom margin and layout, then show plot
    plt.subplots_adjust(bottom=0.15)  # Adjust this value as needed
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plot_file_path = f"{simOutputPath}/flowStats.png"

    plt.savefig(plot_file_path, bbox_inches='tight', dpi=300)
    plt.close()  # Close the plot

def plot_statistics_per_hops(filename, simOutputPath):
    # Initialize dictionaries to store sum and count for each metric per number of hops
    sum_loss_ratios = defaultdict(float)
    sum_mean_delays = defaultdict(float)
    sum_mean_jitters = defaultdict(float)
    sum_nTxPkts = defaultdict(float)
    sum_nRxPkts = defaultdict(float)
    count_per_hops = defaultdict(int)

    # Read the CSV file
    csv_file = f"{simOutputPath}/{filename}"
    with open(csv_file, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip the header row
        for row in csv_reader:
            hops = int(row[5])
            sum_loss_ratios[hops] += float(row[8])
            sum_mean_delays[hops] += float(row[9])
            sum_mean_jitters[hops] += float(row[10])
            sum_nTxPkts[hops] += int(row[6])
            sum_nRxPkts[hops] += int(row[7])
            count_per_hops[hops] += 1

    # Calculate averages
    avg_loss_ratios = {hops: total / count_per_hops[hops] for hops, total in sum_loss_ratios.items()}
    avg_mean_delays = {hops: total / count_per_hops[hops] for hops, total in sum_mean_delays.items()}
    avg_mean_jitters = {hops: total / count_per_hops[hops] for hops, total in sum_mean_jitters.items()}
    avg_nTxPkts = {hops: total / count_per_hops[hops] for hops, total in sum_nTxPkts.items()}
    avg_nRxPkts = {hops: total / count_per_hops[hops] for hops, total in sum_nRxPkts.items()}

    # Now calculate the packet delivery ratio
    avg_packet_delivery_ratio = {hops: (avg_nRxPkts[hops] / avg_nTxPkts[hops]) if avg_nTxPkts[hops] != 0 else 0 for hops in count_per_hops}

  # Sorting keys for consistent plotting
    sorted_hops = sorted(count_per_hops.keys())

    # Create a figure for the plots
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(18, 12))  # Adjusted to 2 rows, 3 columns
    fig.suptitle('Average Path/Flow Statistics per Number of Hops')

    # Define colors for each plot
    plot_colors = ['skyblue', 'tomato', 'deepskyblue', 'mediumseagreen', 'gold', 'orchid']

    # Define the metrics to be plotted
    metrics = [(avg_nTxPkts, avg_nRxPkts), avg_packet_delivery_ratio, avg_loss_ratios, count_per_hops, avg_mean_delays, avg_mean_jitters]
    titles = ['Packet Stats', 'Average Packet Delivery Ratio', 'Average Loss Ratio', 'Flow Count', 'Average Mean Delay (ms)', 'Average Mean Jitter (ms)']

    # Plot each metric
    for i, (metric, title) in enumerate(zip(metrics, titles)):
        ax = axs[i // 3, i % 3]
        if title == 'Packet Stats':
            # Grouped bar chart for nTxPkts and nRxPkts
            width = 0.35
            tx_positions = [x - width/2 for x in sorted_hops]
            rx_positions = [x + width/2 for x in sorted_hops]
            bars_tx = ax.bar(tx_positions, [avg_nTxPkts[h] for h in sorted_hops], width, label='nTxPkts')
            bars_rx = ax.bar(rx_positions, [avg_nRxPkts[h] for h in sorted_hops], width, label='nRxPkts')
            ax.set_title(title)
            ax.legend(loc='lower right')
            
            # Annotate bars for nTxPkts and nRxPkts
            for bar in bars_tx:
                yval = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.2f}', ha='center', va='bottom', rotation=90, fontsize=10)
            for bar in bars_rx:
                yval = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.2f}', ha='center', va='bottom', rotation=90, fontsize=10)
            max_height = max([bar.get_height() for bar in bars_tx])
            ax.set_ylim(0, max_height * 1.15)  # Add extra space at the top
        
        elif title == 'Flow Count':
            # Plot for the count of flows per number of hops
            bars = ax.bar(sorted_hops, [count_per_hops[h] for h in sorted_hops], color=plot_colors[i])
            ax.set_title(title)

            # Annotate bars for flow count
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2, height, f'{int(height)}', ha='center', va='bottom', rotation=90, fontsize=10)
        else:
            bars = ax.bar(sorted_hops, [metric[h] for h in sorted_hops], color=plot_colors[i])
            ax.set_title(title)
            
            # Annotate bars for other metrics
            for bar in bars:
                yval = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.2f}', ha='center', va='bottom', rotation=90, fontsize=10)

            if 'Ratio' in title:  # Set Y-axis limit to 1.0 for ratio plots
                ax.set_ylim(0, 1.1)
            else:
                max_height = max([bar.get_height() for bar in bars])
                ax.set_ylim(0, max_height * 1.15)  # Add extra space at the top
        
        ax.set_xlabel('Number of Hops')
        ax.set_ylabel(title)
        ax.set_xticks(sorted_hops)
    # Adjust layout
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plot_file_path = f"{simOutputPath}/flowStatsPerNHop.png"

    plt.savefig(plot_file_path, bbox_inches='tight', dpi=300)
    plt.close()  # Close the plot


def plot_phy_stats_bars(filename, simOutputPath):


    myFontsize = 14
    colors = plt.get_cmap('tab10').colors

    csv_file = f"{simOutputPath}/{filename}"
    # Read the data from the CSV file
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        data = next(reader)

    # Convert all data to integers
    data = {k: int(v) for k, v in data.items()}

    # Prepare data for Control and Data bar charts
    control_receptions = [data['totalnCorruptRxCtrl'], data['totalnHdRxCtrl']]  # Excluded 'Successful'
    data_receptions = [data['totalnCorruptRxData'], data['totalnHdRxData'], data['totalnIgnoredData']]  # Excluded 'Successful'
    labels_control = ['Corrupted', 'Missed (Half-duplex)']
    labels_data = ['Corrupted', 'Missed (Half-duplex)', 'Ignored (No PSCCH)']
    all_labels_control = ['Successful'] + labels_control  # Include 'Successful' for legend
    all_labels_data = ['Successful'] + labels_data  # Include 'Successful' for legend

    # Plotting
    fig, ax = plt.subplots(1, 2, figsize=(15, 8))

    # Function to calculate percentages
    def calculate_percentages(values, total):
        return [value / total * 100 if total > 0 else 0 for value in values]

    total_control = sum([data['totalnRxCtrl'], data['totalnCorruptRxCtrl'], data['totalnHdRxCtrl']])
    total_data = sum([data['totalnRxData'], data['totalnCorruptRxData'], data['totalnHdRxData'], data['totalnIgnoredData']])
    percentages_control = calculate_percentages(control_receptions, total_control)
    percentages_data = calculate_percentages(data_receptions, total_data)
    all_percentages_control = calculate_percentages([data['totalnRxCtrl']] + control_receptions, total_control)
    all_percentages_data = calculate_percentages([data['totalnRxData']] + data_receptions, total_data)

    # Find the max percentage value and set y-limit accordingly
    max_percentage_control = max(percentages_control)
    max_percentage_data = max(percentages_data)

    # Function to create legend patches
    def create_legend_patches(labels, percentages, colors):
        return [Patch(color=color, label=f'{label} - {percentage:.2f}%') for label, percentage, color in zip(labels, percentages, colors)]

    # Control Bar Chart
    x_control = np.arange(len(labels_control))
    ax[0].bar(x_control, percentages_control, color=colors[1:len(labels_data)+1])  # Start colors from 1 to skip 'Successful'
    ax[0].set_title('PHY Control (PSCCH) Receptions', fontsize=myFontsize)
    ax[0].set_ylabel('Percentage (%)', fontsize=myFontsize)
    ax[0].set_xticks(x_control)
    ax[0].set_xticklabels(labels_control, fontsize=myFontsize)
    control_legend_patches = create_legend_patches(all_labels_control, all_percentages_control, colors[:len(all_labels_control)])
    ax[0].legend(handles=control_legend_patches, fontsize=myFontsize, title="PHY Control (PSCCH) Receptions", title_fontsize=myFontsize)
    ax[0].set_ylim(0, 1.4*max_percentage_control) 

    # Data Bar Chart
    x_data = np.arange(len(labels_data))
    ax[1].bar(x_data, percentages_data, color=colors[1:len(labels_data)+1])  # Start colors from 1 to skip 'Successful'
    ax[1].set_title('PHY Data (PSSCH) Receptions', fontsize=myFontsize)
    ax[1].set_ylabel('Percentage (%)', fontsize=myFontsize)
    ax[1].set_xticks(x_data)
    ax[1].set_xticklabels(labels_data, fontsize=myFontsize)
    data_legend_patches = create_legend_patches(all_labels_data, all_percentages_data, colors[:len(all_labels_data)])
    ax[1].legend(handles=data_legend_patches, fontsize=myFontsize, title="PHY Data (PSSCH) Receptions", title_fontsize=myFontsize)
    ax[1].set_ylim(0, 1.4*max_percentage_data) 

    # Add notes for unaccounted receptions
    note_text_ctrl = f"Unaccounted Control (PSCCH) Rxs: {data['totalnTxCtrl'] - sum([data['totalnRxCtrl']] + control_receptions)}"
    note_text_data = f"Unaccounted Data (PSSCH) Rxs: {data['totalnTxData'] - sum([data['totalnRxData']] + data_receptions)}"
    ax[0].text(0.5, -0.15, note_text_ctrl, ha='center', va='center', fontsize=myFontsize, transform=ax[0].transAxes)
    ax[1].text(0.5, -0.15, note_text_data, ha='center', va='center', fontsize=myFontsize, transform=ax[1].transAxes)

    plt.tight_layout(rect=[0, 0.1, 1, 0.9])
    # Save firgure
    plot_file_path = f"{simOutputPath}/simPhyStats_bars.png"
    plt.savefig(plot_file_path, bbox_inches='tight', dpi=300)
    plt.close()  # Close the plot


#####################################################################################
################################# Main script #######################################
#####################################################################################
    
# Create the parser
parser = argparse.ArgumentParser(description='Plot the stats for a simulation.')

# Add the arguments
parser.add_argument('Path',
                    metavar='path',
                    type=str,
                    help='the path to the run folder')

# Execute the parse_args() method
args = parser.parse_args()

# Now you can use args.Path in your script
simOutputPath = args.Path

plot_nodes_and_paths('nodes_and_paths.csv', simOutputPath)
plot_flow_statistics('flowsStats.csv', simOutputPath)
plot_statistics_per_hops('flowsStats.csv', simOutputPath)
plot_phy_stats_bars ('simPhyStats.csv', simOutputPath)

