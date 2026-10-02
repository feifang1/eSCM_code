# Reproducing the plots and tables

To generate the plots and tables in each section of the [paper](https://arxiv.org/pdf/2508.00223), please run the R scripts in the order listed below.

## Appendix_A.1_code

To run the simulations on a computing cluster by submitting a job:

1. Upload `Appendix_A.1_Simulation_choose_k.R`, `Appendix_A.1_Simulation_choose_k.sh`, and `params_v7.txt` to the same folder.
2. Adjust all input and output paths to match your directory structure.
3. Submit the job using `sbatch Appendix_A.1_Simulation_choose_k.sh`.
4. To run a specific parameter combination `(n, p, shape, lambda, quantile_val)` locally, use the example script `Appendix_A.1_run_example.R`.
5. Run `Appendix_A.1_visualization.R` to generate the plots shown in Figure 8 of the paper.

## Section_4.1_code

To run the simulations on a computing cluster by submitting a job:

1. Upload `Section4.1_simulation_studies.R`, `Section4.1_simulation_studies.sh`, and `param_simulation.txt` to the same folder.
2. Adjust all input and output paths in `Section4.1_simulation_studies.R` and `Section4.1_simulation_studies.sh` to match your directory structure.
3. Submit the job using `sbatch Section4.1_simulation_studies.sh`.
4. To run a specific parameter combination `(n, shape)` locally, use the example script `Section4.1_run_example.R`.
5. Use the simulation results to generate Table 2 of the paper.

## Section_4.2_code

For the left plot in Section 4.2 (ancestral violation rate for the river discharge data), run:

1. `Section4.2_AVR_River_discharge_data.R`
2. `Section4.2_AVR_River_discharge_data_result_visualize.R`

For the right plot in Section 4.2 (pairwise causal direction identification error rate for the river discharge data), run:

1. `Section4.2_pairwise_River_discharge_data.R`
2. `Section4.2_pairwise_River_discharge_data_results_visualize.R`

## Section 4.3

For the plot in Section 4.3, run the following scripts in order:

1. `Section4.3_pairwise_datasets.R`
2. `Section4.3_pairwise_datasets_result_visualize.R`
