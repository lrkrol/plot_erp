# plot_erp
plot_erp is a MATLAB-based script to plot event-related potentials (ERPs) from any number of given epoched datasets (in EEGLAB format), for a single channel. For each ERP curve, any number of datasets can be given. It can optionally calculate and plot a difference wave, standard errors of the mean, and permutation-based statistics (using an included local version of [permutationTest](https://github.com/lrkrol/permutationTest)). Mean curves and statistics can be calculated either within or between datasets.

Sample screenshots:

Two ERPs plus their difference, standard error of the mean for all curves, and statistics highlighting significant differences in grey.

![Screenshot](./plot_erp-diff.png)

ERPs for eight groups of 19 datasets each with custom colour code, labels, and x scale indicator position.

![Screenshot](./plot_erp-mult.png)
