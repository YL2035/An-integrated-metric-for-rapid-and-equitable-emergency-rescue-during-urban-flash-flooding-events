# An-integrated-metric-for-rapid-and-equitable-emergency-rescue-during-urban-flash-flooding-events
This repository provides code for the journal paper "An integrated metric for rapid and equitable emergency rescue during urban flash flooding events." DOI: https://doi.org/10.1016/j.ijdrr.2025.105209
The file contains two main parts, the integrated metric calculation and the rescue boat allocation. For each of these parts, code, data and results are provided.

### Explanation of Repository Structure

```text
./metric-scripts: including all R scripts used for calculating the integrated metric.

    efficiency.R calculates the efficiency of moving from residential buildings to essential facilities. The greater the value, the faster the residents living in these buildings can be rescued to these essential facilities.

    metric_calculation.R is used to derive the integrated metric for both considering and not considering overhead powerlines.

./metric-data: including the data used for metric calculation. The roads in Manville were abstracted as a network composed of nodes and edges.

    ./edge_withpowerline: a folder that includes the accessible road edges, considering the obstruction of overhead powerlines and flooded road sections. The status of road edges at 27 time steps (t1–t26) was provided.

    ./edge_intact: a csv file containing the edge of the road network without disruptions.

    ./node_res_svi: a csv file containing the road network nodes. The nodes are classified as road connections (C), residential buildings (R), and essential facilities (F).

    ./efficiency_withpowerline: a folder includes the efficiency of moving from residential buildings to essential facilities (considering overhead powerline). The status of road edges at 27 time steps (t1–t26) was provided.

    ./efficiency_nopowerline: a csv file that includes the efficiency of moving from residential buildings to essential facilities (without powerline obstruction).

    ./weights: a csv file that provides the varying weights assigned to the three factors of the integrated metric (percentage of building damage, efficiency, and socioeconomic factor).

./rescue-scripts: including all python scripts used for rescue progress simulation based on different lifeboat allocation strategies.

    dynamic_lifeboat_allocation_looppv_3.py is the main script to simulate the rescue progress based on the metrics calculated. The script uses scenarios with “no powerline consideration” / “without considering SVI in metrics” / “Dynamic strategies only” as an example. The aim of this script is to calculate the rescue load variation during the flooding.

    case_compare.py is the script for plots creation.

./rescue-data: including the data used for rescue simulation.

    ./Metric-nopowerline: a folder includes metric calculated for scenario “no powerline consideration”. 37 possible weight combinations of three factors in the integrated metrics are considered, the results are in folder called WID{i} accordingly. In each WID{i} folder, there are metrics calculated for each time step (t1–t26) during flooding, called metric{t}.csv accordingly.

    ./Metric-withpowerline: a folder includes metric calculated for scenario “with powerline consideration”. 37 possible weight combinations of three factors in the integrated metrics are considered, the results are in folder called WID{i} accordingly. In each WID{i} folder, there are metrics calculated for each time step (t1–t26) during flooding, called metric{t}.csv accordingly.

./rescue-output: including plots and results of rescue progress simulation created by scripts in folder rescue-scripts

    ./NoPowerline_consider_wosvi_OnlyDyn: the results for scenario without considering overhead powerlines during rescue planning.

    ./plot: the collections of visualization of rescue load variation.

    ./area_under_curve_heatmap.png: the visualization of the AUC for four allocation strategies among all 37 possible weight combinations of three factors of the integrated metric for specific scenarios

    ./area_under_curves.csv: the calculated AUC for four allocation strategies among all 37 possible weight combinations of three factors of the integrated metric for specific scenarios

    ./reduction_rate_heatmap.png: the average reduction rate for four allocation strategies among all 37 possible weight combinations of three factors of the integrated metric for specific scenarios

    ./total_score_heatmap.png: the score for four allocation strategies among all 37 possible weight combinations of three factors of the integrated metric for specific scenarios

    ./WID{i}: folders include the calculated rescue load for census blocks during flooding, considering 37 possible combinations of three factors of the integrated metric. In each WID{i} folder, the outputs include:

        ./WID{i}/N_Metric_PR_{t}.csv: rescue progress records and rescue load variation for each time step during flooding

        ./WID{i}/adjusted_rescue_load_{strategy name}.csv: the rescue load summary for four allocation strategies

        ./WID{i}/rescue_load.csv: detailed rescue progress record (emergency response facility-strategy)

        ./WID{i}/rescue_load_M.csv: detailed rescue load (emergency response facility-strategy)

        ./WID{i}/N_now: a temp file during rescue simulation

        ./WID{i}/rescue_load_trend_all_strategies_WID{i}.png: the visualization of rescue load variation.

    ./WithPowerline_consider_wosvi_OnlyDyn: the folder includes the results for scenario considering overhead powerlines during rescue planning. Has the same data structure as rescue-output/NoPowerline_consider_wosvi_OnlyDyn

./integrated metrics and components.png: The plot showing the spatial distribution of rescue demands and its three components.

./rescue_load_trend_all_strategies_WID0.png: the variation of rescue load for all scenarios.
```
