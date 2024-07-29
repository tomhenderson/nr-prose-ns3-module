
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import csv
from collections import defaultdict
import argparse
from math import cos, sin, radians
import numpy as np
import pandas as pd
import os 

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

    if os.path.exists(csv_file):
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

    else:
        print (f"File {csv_file} not found!")


def plot_flow_statistics(filename, simOutputPath):
    # Initialize lists to store data
    flow_ids = []
    loss_ratios = []
    mean_delays = []
    mean_jitters = []
    nTxPkts = []
    nRxPkts = []
    nRouteChanges = []
    nMinNHops = []
    nMaxNHops = []
    avgThroughput = []
    # Read the CSV file
    csv_file = f"{simOutputPath}/{filename}"

    if os.path.exists(csv_file):
        with open(csv_file, 'r') as file:
            csv_reader = csv.reader(file)
            next(csv_reader)  # Skip the header row
            for row in csv_reader:
                flow_ids.append(int(row[0]))
                loss_ratios.append(float(row[10]))
                mean_delays.append(float(row[11]))
                mean_jitters.append(float(row[12]))
                avgThroughput.append(float(row[13]))
                nTxPkts.append(int(row[8]))
                nRxPkts.append(int(row[9]))
                nRouteChanges.append(int(row[7]))
                nMinNHops.append(int(row[5]))
                nMaxNHops.append(int(row[6]))

        # Calculate Packet Delivery Ratio
        packet_delivery_ratio = [1 - lr for lr in loss_ratios]
        # Define colors for each plot
        plot_colors = ['skyblue', 'tomato', 'deepskyblue', 'mediumseagreen', 'gold', 'orchid','skyblue', 'tomato','mediumseagreen']

        # Create a figure for the plots
        fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(16, 10))  # Adjusted for 2 rows, 3 columns
        fig.suptitle('Path/Flow Statistics')

        # Plot Packet Stats, Packet Delivery Ratio, Packet Loss Ratio, Number of Hops, Mean Delay, Mean Jitter in the specified order
        plot_order = [
            (nTxPkts, nRxPkts, 'Packet Stats', axs[0, 0]),
            (packet_delivery_ratio, None, 'Packet Delivery Ratio', axs[0, 1]),
            (loss_ratios, None, 'Packet Loss Ratio', axs[0, 2]),
            (mean_delays, None, 'Mean Delay (ms)', axs[0, 3]),
    #        (mean_jitters, None, 'Mean Jitter (ms)', axs[1, 0]),
            (avgThroughput, None, 'Average Throughput (kbps)', axs[1, 0]),
            (nRouteChanges, None, 'Number of Route changes', axs[1, 1]),
            (nMinNHops, None, 'Min number of hops', axs[1, 2]),
            (nMaxNHops, None, 'Max number of hops', axs[1, 3])
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
    else:
        print (f"File {csv_file} not found!")




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

