import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
#from scipy.integrate import simps
from numpy import trapz
import matplotlib.dates as mdates
from datetime import datetime, timedelta

# Define the base path and number of WID folders
base_path_wp = "..\\rescue-output\\WithPowerline_consider_wosvi_OnlyDyn\\"
base_path = "..\\rescue-output\\NoPowerline_consider_wosvi_OnlyDyn\\"
num_wid_folders = 37

# Define the list of strategies
#strategies = ['Equal', 'Fixed1', 'Fixed2', 'Fixed3', 'dyne', 'dyn1', 'dyn2', 'dyn3']
strategies = ['dyne', 'dyn1', 'dyn2', 'dyn3']
strategies1 = ['dyne', 'dyn_1', 'dyn_2', 'dyn_3']

# Define colors for matching strategy pairs
strategy_colors = {
    'Equal': 'blue',
    'dyne': 'blue',
    'Fixed1': 'green',
    'dyn1': 'green',
    'Fixed2': 'red',
    'dyn2': 'red',
    'Fixed3': 'black',
    'dyn3': 'black'
}
start_time = datetime(2021, 9, 2, 1, 0)  # 2021-09-02 01:00
time_interval = timedelta(hours=1)
# Generate a list of datetime labels for the x-axis
time_labels = [start_time + i * time_interval for i in range(26)]
# Generate xticks every 5 hours
xtick_labels = [start_time + timedelta(hours=i) for i in range(0, 26, 5)]
print(xtick_labels)
# Initialize an array to store the reduction rates
reduction_rates = np.zeros((4, num_wid_folders))

# Initialize an array to store the total scores for each strategy and WID case
total_scores = np.zeros((4, num_wid_folders))

# Initialize an array to store the area under the curve for each strategy and WID case
area_under_curves = np.zeros((4, num_wid_folders))
area_under_curves_wp = np.zeros((4, num_wid_folders))

# Compare all strategies for each WID case and plot them together
for wid in range(0, num_wid_folders):
    wid_folder = f"WID{wid}"
    plt_folder = f"plot"
    rescue_load_m_path = os.path.join(base_path, wid_folder, "rescue_load_M.csv")
    rescue_load_m_path_wp = os.path.join(base_path_wp, wid_folder, "rescue_load_M.csv")
    
    # Read the rescue_load_M.csv file
    rescue_load_m_df = pd.read_csv(rescue_load_m_path, header=None)
    rescue_load_m_df.columns = ['Facility_Strategy'] + [str(i) for i in range(26)]  # Set columns from 0 to 25 (representing the time steps)
    rescue_load_m_df_wp = pd.read_csv(rescue_load_m_path_wp, header=None)
    rescue_load_m_df_wp.columns = ['Facility_Strategy'] + [str(i) for i in range(26)]  # Set columns from 0 to 25 (representing the time steps)
    # Split the 'Facility_Strategy' column to extract the strategy name
    rescue_load_m_df[['Facility', 'Strategy']] = rescue_load_m_df['Facility_Strategy'].str.extract(r"\((\d+),\s*'(\w+)'\)")
    rescue_load_m_df_wp[['Facility', 'Strategy']] = rescue_load_m_df_wp['Facility_Strategy'].str.extract(r"\((\d+),\s*'(\w+)'\)")
    
    # Set time step 1 as the average between step 0 and step 2 for all strategies for data fixing, step 1 has issue, treat as missing
    for strategy in strategies:
        strategy_df = rescue_load_m_df[rescue_load_m_df['Strategy'] == strategy]
        total_rescue_load = strategy_df.iloc[:, 1:27].sum(axis=0)  # Sum across all facilities for this strategy at each time step
        total_rescue_load.iloc[1] = (total_rescue_load.iloc[0] + total_rescue_load.iloc[2]) / 2

        strategy_df_wp = rescue_load_m_df_wp[rescue_load_m_df['Strategy'] == strategy]
        total_rescue_load_wp = strategy_df_wp.iloc[:, 1:27].sum(axis=0)  # Sum across all facilities for this strategy at each time step
        total_rescue_load_wp.iloc[1] = (total_rescue_load_wp.iloc[0] + total_rescue_load_wp.iloc[2]) / 2
        
        # Save the adjusted values for later use
        adjusted_rescue_load_path = os.path.join(base_path, wid_folder, f"adjusted_rescue_load_{strategy}.csv")
        total_rescue_load.to_csv(adjusted_rescue_load_path, index=False)
        adjusted_rescue_load_path_wp = os.path.join(base_path_wp, wid_folder, f"adjusted_rescue_load_{strategy}.csv")
        total_rescue_load_wp[0]=total_rescue_load[0]
        total_rescue_load_wp.to_csv(adjusted_rescue_load_path_wp, index=False)
        
        # Calculate the reduction rate as the average reduction from step 13 to step 22
        reduction_rate = ((total_rescue_load.iloc[12:24] - total_rescue_load.iloc[13:25].values) / total_rescue_load.iloc[12]).mean() * 100
        strategy_index = strategies.index(strategy)
        reduction_rates[strategy_index, wid] = reduction_rate
        
        # Plot the total rescue load trend for this strategy
        if wid==0:
            strategy_=strategy.replace("dyn","dyn_")
            print(strategy_)
            strategy_label_np = strategy_.replace("dyn_", "Plan_").replace("e", "equal") + " (NP)"
            strategy_label_wp = strategy_.replace("dyn_", "Plan_").replace("e", "equal") + " (WP)"
            plt.plot(time_labels, total_rescue_load.values, label=strategy_label_np, linestyle='--', color=strategy_colors[strategy])
            plt.plot(time_labels, total_rescue_load_wp.values, label=strategy_label_wp, linestyle='-', color=strategy_colors[strategy])

        # Calculate the area under the curve
        #area_under_curve = simps(total_rescue_load[12:26].values, dx=1)
        #area_under_curves[strategy_index, wid - 1] = area_under_curve
        area_under_curve = trapz(total_rescue_load.values[12:26], dx=1)
        area_under_curve_nor = area_under_curve/(np.max(total_rescue_load.values[12:26])*13)
        area_under_curves[strategy_index, wid] = area_under_curve_nor

        area_under_curve_wp = trapz(total_rescue_load_wp.values[12:26], dx=1)
        area_under_curve_nor_wp = area_under_curve_wp/(np.max(total_rescue_load_wp.values[12:26])*13)
        area_under_curves_wp[strategy_index, wid] = area_under_curve_nor_wp
        
    # Rank the strategies at each time step from 13 to 25
    for t in range(13, 26):
        time_step_values = rescue_load_m_df.iloc[:, t + 1].groupby(rescue_load_m_df['Strategy']).sum()
        sorted_strategies = time_step_values.sort_values().index.tolist()
        
        # Assign scores based on ranking (lowest value gets highest score)
        unique_values = time_step_values.unique()
        unique_values.sort()
        current_rank = 4
        for rank_value in unique_values:
            rank_strategies = time_step_values[time_step_values == rank_value].index.tolist()
            #print(rank_strategies)
            for strategy in rank_strategies:
                strategy_index = strategies.index(strategy)
                total_scores[strategy_index, wid] += current_rank
                
            current_rank -= len(rank_strategies)
        #if t==15:
        #    print(rank_strategies)
        #print(wid, total_scores)
    # Set plot details
    # Add the vertical dashed line
    critical_time = datetime(2021, 9, 2, 13, 0)
    plt.axvline(x=critical_time, color='gray', linestyle='--')
    plt.xlabel("Datetime")
    plt.xticks(xtick_labels, rotation=20)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))
    plt.ylabel('Total Rescue Load (Summed Across All Facilities)')
    plt.title("Sum of Rescue Load Over Time for Different Plans and Scenarios")
    plt.legend()
    plt.grid(True)
    #plt.xticks(range(0, 26))  # Ensure x-axis shows time steps from 0 to 25
    #plt.xlim(0, 25)  # Limit x-axis to show only time steps from 12 to 25
    plt.ylim(0, 55)

    # Save the plot for this WID
    if wid==0:
        output_plot_path = os.path.join("..\\rescue-output\\", f'rescue_load_trend_all_strategies_WID0.png')
        plt.savefig(output_plot_path)
    
    #output_plot_path = os.path.join(base_path, wid_folder, f'rescue_load_trend_all_strategies_WID{wid}.png')
    #output_plot_path2 = os.path.join(base_path, plt_folder, f'rescue_load_trend_all_strategies_WID{wid}.png')
    #plt.savefig(output_plot_path)
    #plt.savefig(output_plot_path2)
    #plt.close()

# Create a heatmap to visualize the reduction rates
plt.figure(figsize=(14, 7))
sns.heatmap(reduction_rates, annot=True, cmap="viridis", yticklabels=strategies, xticklabels=[f"WID{wid}" for wid in range(0, num_wid_folders)], square=True, fmt=".1f", cbar_kws={'shrink': 0.8}, annot_kws={"size": 8})
plt.xlabel('Strategies')
plt.ylabel('WID Cases')
plt.title('Average Reduction Rate of Rescue Load M for Different Strategies and WID Cases (Time Steps 13-22)')

# Save the heatmap plot
heatmap_output_path = os.path.join(base_path, 'reduction_rate_heatmap.png')
plt.savefig(heatmap_output_path)
plt.close()

# Create a heatmap to visualize the total scores for each strategy and WID case
plt.figure(figsize=(14, 7))
sns.heatmap(total_scores, annot=True, cmap="coolwarm", yticklabels=strategies, xticklabels=[f"WID{wid}" for wid in range(0, num_wid_folders)], square=True, fmt=".0f", cbar_kws={'shrink': 0.8}, annot_kws={"size": 8})
plt.xlabel('Strategies')
plt.ylabel('WID Cases')
plt.title('Total Scores for Each Strategy and WID Case (Time Steps 13-25)')

# Save the heatmap plot for total scores
total_score_heatmap_output_path = os.path.join(base_path, 'total_score_heatmap.png')
plt.savefig(total_score_heatmap_output_path)
plt.close()

# Create a heatmap to visualize the area under the curve for each strategy and WID case
plt.figure(figsize=(14, 7))
sns.heatmap(area_under_curves, annot=True, cmap="plasma", yticklabels=strategies, xticklabels=[f"WID{wid}" for wid in range(0, num_wid_folders)], square=True, fmt=".3f", cbar_kws={'shrink': 0.8}, annot_kws={"size": 5})
plt.ylabel('Strategies')
plt.xlabel('WID Cases')
plt.title('Area Under the Curve for Different Strategies and WID Cases')

# Save the heatmap plot for area under the curve
area_under_curve_heatmap_output_path = os.path.join(base_path, 'area_under_curve_heatmap.png')
plt.savefig(area_under_curve_heatmap_output_path, dpi=300)
plt.close()

plt.figure(figsize=(14, 7))
sns.heatmap(area_under_curves_wp, annot=True, cmap="plasma", yticklabels=strategies, xticklabels=[f"WID{wid}" for wid in range(0, num_wid_folders)], square=True, fmt=".3f", cbar_kws={'shrink': 0.8}, annot_kws={"size": 5})
plt.ylabel('Strategies')
plt.xlabel('WID Cases')
plt.title('Area Under the Curve for Different Strategies and WID Cases')

# Save the heatmap plot for area under the curve
area_under_curve_heatmap_output_path_wp = os.path.join(base_path_wp, 'area_under_curve_heatmap.png')
plt.savefig(area_under_curve_heatmap_output_path_wp, dpi=300)
plt.close()

# Export area under curves to an Excel file
area_under_curves_df = pd.DataFrame(area_under_curves, index=strategies, columns=[f"WID{wid}" for wid in range(0, num_wid_folders)])
area_under_curves_excel_output_path = os.path.join(base_path, 'area_under_curves.xlsx')
area_under_curves_df.to_excel(area_under_curves_excel_output_path)