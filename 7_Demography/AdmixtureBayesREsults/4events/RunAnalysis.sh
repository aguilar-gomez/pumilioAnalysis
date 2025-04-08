 
python3.10 ../../AdmixtureBayes-main/admixturebayes/analyzeSamples.py  --mcmc_results  Frog3/output1.csv  --burn_in_fraction 0.5

python3.10 ../../AdmixtureBayes-main/admixturebayes/makePlots.py --plot estimates  --posterior thinned_samples.csv  --write_rankings chain3.txt

python3.10 ../../AdmixtureBayes-main/admixturebayes/makePlots.py --plot top_minimal_topologies   --posterior thinned_samples.csv  --write_rankings chain3_minimal.txt