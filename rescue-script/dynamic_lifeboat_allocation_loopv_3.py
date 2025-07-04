"""
Dynamic Lifeboat Allocation System
This script implements a dynamic allocation system for lifeboats across multiple facilities
based on rescue demand metrics and building damage probabilities.
"""

import os
import pandas as pd
from collections import Counter
import numpy as np

# =============================================================================
# CORE CALCULATION FUNCTIONS
# =============================================================================

def calculate_initial_N(metric_df):
    """
    Calculate initial N values (eligibility for rescue) based on damage probability and distance.
    N is 1 if both Pdamage_avg and inv_dist_avg are non-zero, otherwise 0.
    
    Args:
        metric_df (pd.DataFrame): DataFrame containing building metrics
        
    Returns:
        pd.DataFrame: DataFrame with N column added
    """
    metric_df['N'] = np.where((metric_df['Pdamage_avg'] * metric_df['inv_dist_avg']) != 0, 1, 0)
    return metric_df

def calculate_M(metric_df):
    """
    Calculate M values (rescue metrics) by multiplying N values with priority metrics.
    Handles multiple allocation strategies (fixed and dynamic ratios).
    
    Args:
        metric_df (pd.DataFrame): DataFrame containing N values and metrics
        
    Returns:
        pd.DataFrame: DataFrame with M columns added for different strategies
    """
    # Remove duplicate N columns if they exist
    if 'N' in metric_df.columns:
        metric_df = metric_df.loc[:, ~metric_df.columns.duplicated()]

    # Calculate M values for different allocation strategies
    strategy_mappings = {
        'N_fixed3': ('M_fixed3', 'metric_PR'),
        'N_dyne': ('M_dyne_wosvi', 'metric_wosvi_PR'),
        'N_dyn1': ('M_dyn1_wosvi', 'metric_wosvi_PR'),
        'N_dyn2': ('M_dyn2_wosvi', 'metric_wosvi_PR'),
        'N_dyn3': ('M_dyn3', 'metric_PR'),
        'N': ('M', 'metric_PR'),
        'N': ('M_wosvi', 'metric_wosvi_PR')
    }
    
    for n_col, (m_col, metric_col) in strategy_mappings.items():
        if n_col in metric_df.columns and metric_col in metric_df.columns:
            metric_df[m_col] = metric_df[n_col] * metric_df[metric_col]
    
    return metric_df

# =============================================================================
# LIFEBOAT ALLOCATION FUNCTIONS
# =============================================================================

def initial_allocation(strategy):
    """
    Define initial fixed allocation strategies for the three facilities.
    
    Args:
        strategy (str): Strategy name ('equal', 'ratio1', 'ratio2', 'ratio3')
        
    Returns:
        dict: Mapping of facility ID to number of lifeboats
    """
    initial_ratios = {
        'equal': {7413: 2, 7537: 2, 9689: 2},    # Equal distribution
        'ratio1': {7413: 1, 7537: 3, 9689: 2},   # Based on block demand
        'ratio2': {7413: 3, 7537: 2, 9689: 1},   # Based on associated blocks
        'ratio3': {7413: 2, 7537: 3, 9689: 1}    # Based on summed metrics
    }
    return initial_ratios[strategy].copy()

def dynamic_allocation(lifeboats_per_facility, rescue_loads_M):
    """
    Dynamically reallocate lifeboats based on current rescue demand.
    Facilities with zero rescue load release their lifeboats for reallocation.
    
    Args:
        lifeboats_per_facility (dict): Current lifeboat allocation per facility
        rescue_loads_M (dict): Current rescue demand per facility
        
    Returns:
        dict: Updated lifeboat allocation per facility
    """
    released_lifeboats = 0
    facilities_with_load = []

    # Release lifeboats from facilities with zero rescue load
    for fid, load in rescue_loads_M.items():
        if load == 0 and lifeboats_per_facility[fid] > 0:
            released_lifeboats += lifeboats_per_facility[fid]
            lifeboats_per_facility[fid] = 0
        elif load > 0:
            facilities_with_load.append(fid)

    # If no facilities have rescue demand, return current allocation
    if not facilities_with_load:
        return lifeboats_per_facility

    # Calculate proportional allocation based on rescue loads
    remaining_ratios = {fid: rescue_loads_M[fid] for fid in facilities_with_load}
    total_ratio = sum(remaining_ratios.values())
    
    if total_ratio > 0:
        remaining_ratios = {fid: load / total_ratio for fid, load in remaining_ratios.items()}

    # Allocate released lifeboats proportionally
    lifeboats_allocated = {fid: int(released_lifeboats * remaining_ratios[fid]) 
                          for fid in facilities_with_load}

    # Distribute any remaining lifeboats to facility with highest ratio
    total_allocated = sum(lifeboats_allocated.values())
    remaining_lifeboats = released_lifeboats - total_allocated

    while remaining_lifeboats > 0:
        highest_ratio_facility = max(remaining_ratios, key=lambda fid: remaining_ratios[fid])
        lifeboats_allocated[highest_ratio_facility] += 1
        remaining_lifeboats -= 1

    # Update final allocation
    for fid in facilities_with_load:
        lifeboats_per_facility[fid] += lifeboats_allocated[fid]

    return lifeboats_per_facility

def assign_lifeboats(metric_df, plan, ratio=None):
    """
    Assign lifeboats to specific buildings based on allocation plan.
    Uses metric_PR for priority ranking.
    
    Args:
        metric_df (pd.DataFrame): Building metrics DataFrame
        plan (str): Allocation plan ('fixed_ratio3' or 'dyn_ratio3')
        ratio (dict): Number of lifeboats per facility
        
    Returns:
        pd.DataFrame: Updated DataFrame with lifeboat assignments
    """
    for fid in [7413, 7537, 9689]:
        num_lifeboats = ratio.get(fid, 1)
        
        if plan == 'fixed_ratio3':
            metric_col = 'M_fixed3' if 'M_fixed3' in metric_df.columns else 'M'
            n_col = 'N_fixed3'
        elif plan == 'dyn_ratio3':
            metric_col = 'M_dyn3' if 'M_dyn3' in metric_df.columns else 'M'
            n_col = 'N_dyn3'
        else:
            continue
            
        # Select top buildings based on metric values
        top_blocks = metric_df[metric_df['Nearest_FID'] == fid].sort_values(
            by=metric_col, ascending=False).head(num_lifeboats)
        metric_df.loc[metric_df['BID'].isin(top_blocks['BID']), n_col] = -1

    return metric_df

def assign_lifeboats_new(metric_df, plan, ratio=None):
    """
    Assign lifeboats to specific buildings based on allocation plan.
    Uses metric_wosvi_PR for priority ranking (without SVI consideration).
    
    Args:
        metric_df (pd.DataFrame): Building metrics DataFrame
        plan (str): Allocation plan name
        ratio (dict): Number of lifeboats per facility
        
    Returns:
        pd.DataFrame: Updated DataFrame with lifeboat assignments
    """
    plan_mappings = {
        'equal': ('M_equal_wosvi', 'N_equal', 2),
        'fixed_ratio1': ('M_fixed1_wosvi', 'N_fixed1', None),
        'fixed_ratio2': ('M_fixed2_wosvi', 'N_fixed2', None),
        'dyn_ratioe': ('M_dyne_wosvi', 'N_dyne', None),
        'dyn_ratio1': ('M_dyn1_wosvi', 'N_dyn1', None),
        'dyn_ratio2': ('M_dyn2_wosvi', 'N_dyn2', None)
    }
    
    if plan not in plan_mappings:
        return metric_df
        
    metric_col, n_col, fixed_boats = plan_mappings[plan]
    fallback_metric = 'M_wosvi'
    
    for fid in [7413, 7537, 9689]:
        num_lifeboats = fixed_boats if fixed_boats else ratio.get(fid, 1)
        
        # Choose appropriate metric column
        use_metric = metric_col if metric_col in metric_df.columns else fallback_metric
        
        # Select top buildings for lifeboat assignment
        top_blocks = metric_df[metric_df['Nearest_FID'] == fid].sort_values(
            by=use_metric, ascending=False).head(num_lifeboats)
        metric_df.loc[metric_df['BID'].isin(top_blocks['BID']), n_col] = -1

    return metric_df

# =============================================================================
# RESCUE LOAD CALCULATION FUNCTIONS
# =============================================================================

def calculate_rescue_loads_fixblockM(num_files, base_path1):
    """
    Calculate rescue loads for all dynamic allocation strategies across time periods.
    Analyzes building counts and total rescue metrics for each facility and strategy.
    
    Args:
        num_files (int): Number of time periods to analyze
        base_path1 (str): Base path for output files
        
    Returns:
        tuple: (rescue_loads, rescue_loads_M) - building counts and metric sums
    """
    # Initialize data structures for tracking rescue loads
    rescue_loads = {fid: {'dyne': [], 'dyn1': [], 'dyn2': [], 'dyn3': []} 
                   for fid in [7413, 7537, 9689]}
    rescue_loads_M = {fid: {'dyne': [], 'dyn1': [], 'dyn2': [], 'dyn3': []} 
                     for fid in [7413, 7537, 9689]}

    for i in range(1, num_files + 1):
        n_metric_df = pd.read_csv(f'{base_path1}\\N_Metric_PR_{i}.csv')
        
        for fid in [7413, 7537, 9689]:
            if i > 13:  # Post-event period with specific allocation tracking
                # Count buildings requiring rescue for each strategy
                fid_load_dyne = n_metric_df[(n_metric_df['Nearest_FID'] == fid) & 
                                          (n_metric_df['N_Metric_dyne'] > 0)].shape[0]
                fid_load_dyn1 = n_metric_df[(n_metric_df['Nearest_FID'] == fid) & 
                                          (n_metric_df['N_Metric_dyn1'] > 0)].shape[0]
                fid_load_dyn2 = n_metric_df[(n_metric_df['Nearest_FID'] == fid) & 
                                          (n_metric_df['N_Metric_dyn2'] > 0)].shape[0]
                fid_load_dyn3 = n_metric_df[(n_metric_df['Nearest_FID'] == fid) & 
                                          (n_metric_df['N_Metric_dyn3'] > 0)].shape[0]

                # Calculate total rescue metrics for each strategy
                M_load_dyne = n_metric_df.loc[(n_metric_df['Nearest_FID'] == fid) & 
                                            (n_metric_df['N_Metric_dyne'] > 0), 'block_M'].sum()
                M_load_dyn1 = n_metric_df.loc[(n_metric_df['Nearest_FID'] == fid) & 
                                            (n_metric_df['N_Metric_dyn1'] > 0), 'block_M'].sum()
                M_load_dyn2 = n_metric_df.loc[(n_metric_df['Nearest_FID'] == fid) & 
                                            (n_metric_df['N_Metric_dyn2'] > 0), 'block_M'].sum()
                M_load_dyn3 = n_metric_df.loc[(n_metric_df['Nearest_FID'] == fid) & 
                                            (n_metric_df['N_Metric_dyn3'] > 0), 'block_M'].sum()
            else:  # Pre-event period with general metrics
                # Remove duplicate columns and use general M values
                n_metric_df = n_metric_df.loc[:, ~n_metric_df.columns.duplicated(keep='last')]
                
                fid_load_dyne = fid_load_dyn1 = fid_load_dyn2 = fid_load_dyn3 = \
                    n_metric_df[(n_metric_df['Nearest_FID'] == fid) & (n_metric_df['M'] > 0)].shape[0]
                
                M_load_dyne = M_load_dyn1 = M_load_dyn2 = M_load_dyn3 = \
                    n_metric_df.loc[(n_metric_df['Nearest_FID'] == fid) & (n_metric_df['M'] > 0), 'block_M'].sum()

            # Store results
            rescue_loads[fid]['dyne'].append(fid_load_dyne)
            rescue_loads[fid]['dyn1'].append(fid_load_dyn1)
            rescue_loads[fid]['dyn2'].append(fid_load_dyn2)
            rescue_loads[fid]['dyn3'].append(fid_load_dyn3)

            rescue_loads_M[fid]['dyne'].append(M_load_dyne)
            rescue_loads_M[fid]['dyn1'].append(M_load_dyn1)
            rescue_loads_M[fid]['dyn2'].append(M_load_dyn2)
            rescue_loads_M[fid]['dyn3'].append(M_load_dyn3)

    # Save results to CSV files
    rescue_loads_df = pd.DataFrame.from_dict({(i, j): rescue_loads[i][j]
                                              for i in rescue_loads.keys()
                                              for j in rescue_loads[i].keys()},
                                             orient='index')
    rescue_loads_df.to_csv(f'{base_path1}\\rescue_load.csv')

    rescue_loads_M_df = pd.DataFrame.from_dict({(i, j): rescue_loads_M[i][j]
                                              for i in rescue_loads_M.keys()
                                              for j in rescue_loads_M[i].keys()},
                                             orient='index')
    rescue_loads_M_df.to_csv(f'{base_path1}\\rescue_load_M.csv')

    return rescue_loads, rescue_loads_M

# =============================================================================
# MAIN PROCESSING FUNCTION
# =============================================================================

def process_files(num_files=26, wid_range=range(0, 37)):
    """
    Main function to process all scenarios and time periods.
    Implements dynamic lifeboat allocation across multiple weather scenarios.
    
    Args:
        num_files (int): Number of time periods to process (default: 26)
        wid_range (range): Range of weather scenario IDs to process (default: 0-36)
    """
    # Define file paths
    base_path = "..\\rescue-data\\"
    # example path for scenario nopowerline --Dynamic strategies only
    base_path1 = "..\\rescue-output\\NoPowerline_consider_wosvi_OnlyDyn\\"
    base_path2 = "..\\rescue-data\\Metric-nopowerline\\"
    base_path3 = "..\\rescue-data\\metric-nopowerline-wosvi\\"
    
    # Load building-facility association data
    bid_fid_association_path = base_path + "BID_FID_association.csv"
    bid_fid_association_df = pd.read_csv(bid_fid_association_path)
    
    # Define fixed allocation strategies
    fixed_ratio3 = {7413: 2, 7537: 3, 9689: 1}  # Based on summed metrics

    # Process each weather scenario
    for y in wid_range:
        print(f"Processing Weather Scenario {y}")
        
        # Create output directory
        output_folder = os.path.join(base_path1, f"WID{y}")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        
        # Load baseline block metrics
        base_metric_path_blockM = os.path.join(base_path2, f"WID{y}\\", f"metric13.csv")
        block_M = pd.read_csv(base_metric_path_blockM)

        # Initialize tracking variables for dynamic allocation
        updated_N_values_dyne = {}
        updated_N_values_dyn1 = {}
        updated_N_values_dyn2 = {}
        updated_N_values_dyn3 = {}
        N_now_df = pd.DataFrame()

        # Process each time period
        for i in range(1, num_files + 1):
            # Load input data
            base_metric_path = os.path.join(base_path2, f"WID{y}\\", f"metric{i}.csv")
            metric_wosvi_path = os.path.join(base_path3, f"metric_wosvi_TID_{i}.csv")
            
            metric_df = pd.read_csv(base_metric_path)
            metric_wosvi = pd.read_csv(metric_wosvi_path)
            
            # Merge data
            metric_df = pd.merge(metric_df, bid_fid_association_df, on='BID', how='left')
            metric_df['metric_wosvi_PR'] = metric_wosvi['metric_wosvi_PR']
            metric_df['block_M'] = block_M['metric_PR']

            # Process based on time period
            if i == 13:  # Event occurrence time
                # Initialize dynamic allocation strategies
                dyn_ratioe = initial_allocation('equal')
                dyn_ratio1 = initial_allocation('ratio1')
                dyn_ratio2 = initial_allocation('ratio2')
                dyn_ratio3 = initial_allocation('ratio3')
                
                # Calculate initial values
                metric_df_initial = calculate_initial_N(metric_df)
                metric_df = calculate_M(metric_df_initial)
                metric_df['N_now'] = np.where((metric_df['Pdamage_avg'] * metric_df['inv_dist_avg']) != 0, 1, 0)
                N_now_df[f'N_now_{i}'] = metric_df['N_now']
                
                # Assign lifeboats using initial strategies
                metric_df_dyne = assign_lifeboats(metric_df.copy(), 'dyn_ratioe', dyn_ratioe)
                metric_df_dyn1 = assign_lifeboats(metric_df.copy(), 'dyn_ratio1', dyn_ratio1)
                metric_df_dyn2 = assign_lifeboats(metric_df.copy(), 'dyn_ratio2', dyn_ratio2)
                metric_df_dyn3 = assign_lifeboats(metric_df.copy(), 'dyn_ratio3', dyn_ratio3)
                
                # Save output
                output_df = metric_df[['BID', 'Nearest_FID', 'M', 'M', 'block_M']]
                output_df.to_csv(os.path.join(output_folder, f'N_Metric_PR_{i}.csv'), index=False)
                
                # Prepare for next iteration
                if i < num_files:
                    metric_df_dyne['N_dyne'] = metric_df['N']
                    metric_df_dyn1['N_dyn1'] = metric_df['N']
                    metric_df_dyn2['N_dyn2'] = metric_df['N']
                    metric_df_dyn3['N_dyn3'] = metric_df['N']
                    
                    updated_N_values_dyne[i + 1] = metric_df_dyne['N_dyne']
                    updated_N_values_dyn1[i + 1] = metric_df_dyn1['N_dyn1']
                    updated_N_values_dyn2[i + 1] = metric_df_dyn2['N_dyn2']
                    updated_N_values_dyn3[i + 1] = metric_df_dyn3['N_dyn3']
                    
            elif i > 13:  # Post-event dynamic allocation
                # Update N values based on current conditions and previous allocations
                metric_df['N_now'] = np.where((metric_df['Pdamage_avg'] * metric_df['inv_dist_avg']) != 0, 1, 0)
                
                metric_df['N_dyne'] = updated_N_values_dyne.get(i, metric_df_dyne['N_dyne'])
                metric_df['N_dyne'] = np.where(metric_df['N_dyne'] != -1, metric_df['N_now'], -1)
                
                metric_df['N_dyn1'] = updated_N_values_dyn1.get(i, metric_df_dyn1['N_dyn1'])
                metric_df['N_dyn1'] = np.where(metric_df['N_dyn1'] != -1, metric_df['N_now'], -1)
                
                metric_df['N_dyn2'] = updated_N_values_dyn2.get(i, metric_df_dyn2['N_dyn2'])
                metric_df['N_dyn2'] = np.where(metric_df['N_dyn2'] != -1, metric_df['N_now'], -1)
                
                metric_df['N_dyn3'] = updated_N_values_dyn3.get(i, metric_df_dyn3['N_dyn3'])
                metric_df['N_dyn3'] = np.where(metric_df['N_dyn3'] != -1, metric_df['N_now'], -1)
                
                metric_df = calculate_M(metric_df)
                N_now_df[f'N_now_{i}'] = metric_df['N_now']

                # Calculate rescue loads for dynamic reallocation
                rescue_loads_Me = {fid: 0 for fid in initial_allocation('equal').keys()}
                rescue_loads_M1 = {fid: 0 for fid in initial_allocation('ratio1').keys()}
                rescue_loads_M2 = {fid: 0 for fid in initial_allocation('ratio2').keys()}
                rescue_loads_M3 = {fid: 0 for fid in initial_allocation('ratio3').keys()}
                
                for fid in [7413, 7537, 9689]:
                    # Calculate rescue demand for each strategy
                    rescue_loads_Me[fid] = metric_df[(metric_df['Nearest_FID'] == fid) & 
                                                   (metric_df['N_dyne'] > 0)].shape[0]
                    if rescue_loads_Me[fid] != 0:
                        rescue_loads_Me[fid] = 1 * dyn_ratioe[fid]
                        
                    rescue_loads_M1[fid] = metric_df[(metric_df['Nearest_FID'] == fid) & 
                                                   (metric_df['N_dyn1'] > 0)].shape[0]
                    if rescue_loads_M1[fid] != 0:
                        rescue_loads_M1[fid] = 1 * dyn_ratio1[fid]
                        
                    rescue_loads_M2[fid] = metric_df[(metric_df['Nearest_FID'] == fid) & 
                                                   (metric_df['N_dyn2'] > 0)].shape[0]
                    if rescue_loads_M2[fid] != 0:
                        rescue_loads_M2[fid] = 1 * dyn_ratio2[fid]
                        
                    rescue_loads_M3[fid] = metric_df[(metric_df['Nearest_FID'] == fid) & 
                                                   (metric_df['N_dyn3'] > 0)]['M_dyn3'].sum()

                # Update dynamic allocations or fall back to fixed strategies
                dyn_ratioe = (fixed_ratio3 if sum(rescue_loads_Me.values()) == 0 
                             else dynamic_allocation(initial_allocation('equal'), rescue_loads_Me))
                
                dyn_ratio1 = (fixed_ratio3 if sum(rescue_loads_M1.values()) == 0 
                             else dynamic_allocation(initial_allocation('ratio1'), rescue_loads_M1))
                
                dyn_ratio2 = (fixed_ratio3 if sum(rescue_loads_M2.values()) == 0 
                             else dynamic_allocation(initial_allocation('ratio2'), rescue_loads_M2))
                
                dyn_ratio3 = (fixed_ratio3 if sum(rescue_loads_M3.values()) == 0 
                             else dynamic_allocation(initial_allocation('ratio3'), rescue_loads_M3))

                # Assign lifeboats with updated allocations
                metric_df_dyne = assign_lifeboats_new(metric_df.copy(), 'dyn_ratioe', dyn_ratioe)
                metric_df_dyn1 = assign_lifeboats_new(metric_df.copy(), 'dyn_ratio1', dyn_ratio1)
                metric_df_dyn2 = assign_lifeboats_new(metric_df.copy(), 'dyn_ratio2', dyn_ratio2)
                metric_df_dyn3 = assign_lifeboats(metric_df.copy(), 'dyn_ratio3', dyn_ratio3)
                
                # Calculate final metrics
                metric_df['N_Metric_dyne'] = metric_df_dyne['N_dyne'] * metric_df_dyne['metric_PR']
                metric_df['N_Metric_dyn1'] = metric_df_dyn1['N_dyn1'] * metric_df_dyn1['metric_PR']
                metric_df['N_Metric_dyn2'] = metric_df_dyn2['N_dyn2'] * metric_df_dyn2['metric_PR']
                metric_df['N_Metric_dyn3'] = metric_df_dyn3['N_dyn3'] * metric_df_dyn3['metric_PR']
                
                # Save output
                output_df = metric_df[['BID', 'Nearest_FID', 'N_Metric_dyne', 'N_Metric_dyn1', 
                                     'N_Metric_dyn2', 'N_Metric_dyn3', 'block_M']]
                output_df.to_csv(os.path.join(output_folder, f'N_Metric_PR_{i}.csv'), index=False)

                # Prepare for next iteration
                if i < num_files:
                    updated_N_values_dyne[i + 1] = metric_df_dyne['N_dyne']
                    updated_N_values_dyn1[i + 1] = metric_df_dyn1['N_dyn1']
                    updated_N_values_dyn2[i + 1] = metric_df_dyn2['N_dyn2']
                    updated_N_values_dyn3[i + 1] = metric_df_dyn3['N_dyn3']

            else:  # Pre-event period (i < 13)
                metric_df['N'] = np.where((metric_df['Pdamage_avg'] * metric_df['inv_dist_avg']) != 0, 1, 0)
                metric_df = calculate_M(metric_df)
                N_now_df[f'N_now_{i}'] = metric_df['N']
                
                # Save output
                output_df = metric_df[['BID', 'Nearest_FID', 'N', 'M', 'block_M']]
                output_df.to_csv(os.path.join(output_folder, f'N_Metric_PR_{i}.csv'), index=False)

        # Save tracking data and calculate rescue loads
        N_now_df.to_csv(os.path.join(output_folder, 'N_now.csv'), index=False)
        calculate_rescue_loads_fixblockM(num_files, output_folder)
        
        print(f"Completed Weather Scenario {y}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    process_files()