

Instruction to get the E-distance results:

1. Run the script "E-dist.R" to calculate the E-distance between clusters for different resoltions. We used sampling to account for unbalanced dataset sizes. Input are the filtered_feature_bc_matrices from the sequencing runs. Output are csv files with the E-distance between cluster i and cluster j for different resolutions and multiple sampling runs.
2. Aggregate the results with the Jupyter notebook "Edist_createCSVPlot.ipynb". Input are the files generated in step1. Output are aggregated csv files, which can be used to plot the results.
3. Execute the script "plot_Edist_YG.R" to plot the results. Input are the files generated in step2. Output are plots in svg format.
