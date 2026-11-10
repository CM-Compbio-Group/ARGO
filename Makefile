.PHONY: all setup simulate analyze plot clean

# set up computational environment 
setup: .setup_complete

.setup_complete: renv.lock .Rprofile renv/activate.R 
	@echo "setting up the environment"
	Rscript -e "renv::restore(prompt=FALSE)"
	touch .setup_complete

all: simulate analyze plot

simulate: .simulation_complete

.simulation_complete: testing/data/raw/PANCREAS_sce.rds
	@echo "fitting model parameters and simulating counts data"
	Rscript testing/src/scripts/simulate_data.R
	touch .simulation_complete

analyze: .analysis_complete

.analysis_complete: testing/data/temp/experiment_params.csv .simulation_complete
	@echo "now running differential expression analysis between conditions"
	Rscript testing/src/scripts/analyze_data.R
	touch .analysis_complete

plot: .plotting_complete
.plotting_complete: .analysis_complete testing/output/diff_testing_eval_metrics.rds
	@echo "now plotting evaluation metrics"
	Rscript testing/src/scripts/plot_results.R
	touch .plotting_complete

clean:
	rm .setup_complete ._simulation_complete .analyze_complete .plotting_Complete testing/output/*.rds testing/figs/*.pdf testing/data/temp/*.rds
