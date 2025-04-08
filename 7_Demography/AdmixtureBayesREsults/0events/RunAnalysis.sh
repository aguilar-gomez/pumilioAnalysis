 
python3.10 ../AdmixtureBayes-main/admixturebayes/analyzeSamples.py  --mcmc_results  Frog2/output1.csv  --burn_in_fraction 0.3

python3.10 ../AdmixtureBayes-main/admixturebayes/makePlots.py --plot top_trees --posterior thinned_samples.csv  --write_rankings chain3subs.txt
