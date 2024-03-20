#!/usr/bin/env python3
# coding: utf-8
# Author : Diana Aguilar

import os
# must set these before loading numpy
os.environ["OMP_NUM_THREADS"] = '12' # or however many threads you want
os.environ["OPENBLAS_NUM_THREADS"] = '12'
os.environ["MKL_NUM_THREADS"] = '12'

import argparse
import dadi
import nlopt
import numpy
import matplotlib.pyplot as plt
import pandas as pd
import moments

parser = argparse.ArgumentParser()
parser.add_argument("n_sim", type=int, help="number of simulations")
parser.add_argument("pop1", type=str, help="name of population 1")
parser.add_argument("pop2", type=str, help="name of population 2")

args = parser.parse_args()

dataset = args.pop1 + args.pop2
inputfs = dadi.Spectrum.from_file(dataset+".2dfs")

#Calculate Fst
fst=inputfs.Fst()

print("Analizying",args.pop1, args.pop2, "with Fst:", fst, "\n number of ompimizations:",args.n_sim)

#Retrive the sample sizes from the data
ns = inputfs.sample_sizes

#Grid sizes
# Define the grid points based on the sample size.
# For smaller data (largest sample size is about <100) [ns+20, ns+30, ns+40] is a good starting point.
# For larger data (largest sample size is about >=100) or for heavily down projected data [ns+100, ns+110, ns+120] is a good starting point.
pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

demo_model = dadi.Demographics2D.IM

# Wrap the demographic model in a function that utilizes grid points
#which increases dadi's ability to more accurately generate a model frequency spectrum.
demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

# Define starting parameters
#First two parameters effective pop size, 
#Times are given in units of 2Nref (reference population size) generations
#s,nu1,nu2,T,m12,m21 = params

params = [.5,10, 10, .5, 0,0]

# Define boundaries of optimization.
# It is a good idea to have boundaries to avoid optimization
# from trying parameter sets that are time consuming without
# nessicarily being correct.
# If optimization infers parameters very close to the boundaries, we should increase them.
lower_bounds = [1e-3, 1e-3, 1e-3, 1e-3, 0, 0]
upper_bounds = [1,2000, 2000, 3000, 0,0]

try:
  fid = open('results/'+dataset+'_demoIM_fits.txt','a')
except:
  fid = open('results/'+dataset+'_demoIM_fits.txt','w')

for i in range(args.n_sim):
    print("Optimization number: ",i+1)
    # Perturb parameters
    # Optimizers dadi uses are mostly deterministic
    # so we will want to randomize parameters for each optimization.
    p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,
                              lower_bound=lower_bounds)

    # Run optimization
    # At the end of the optimization we will get the
    # optimal parameters and log-likelihood.
    # Verbose is how many evaluations are done before pàrameter is printed

    popt, ll_model = dadi.Inference.opt(p0, inputfs, demo_model_ex, pts_l,
                                    lower_bound=lower_bounds,
                                    upper_bound=upper_bounds,
                                    algorithm=nlopt.LN_BOBYQA,
                                    maxeval=1000, verbose=0)

    # Calculate the synonymous theta
    model_fs = demo_model_ex(popt, ns, pts_l)
    theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, inputfs)

    # Write results to fid
    res = [ll_model] + list(popt) + [theta0]
    fid.write('\t'.join([str(ele) for ele in res])+'\n')
    
fid.close()

#Read results
results=pd.read_table('results/'+dataset+'_demoIM_fits.txt', names=["ll","s","nu1","nu2","T","m12","m21","theta"])

#Select model with maximum likelihood
model_sel=results.loc[results.ll.idxmax()]

print("Best model", model_sel)

'''
The Nref is calculated by the following equation:
Theta = 4 * Nref * mu * L

Nref=Theta/(4*mu*L)

L is the total length of DNA sequence (in bp) that was analyzed in order to obtain the SNP data

Theta is estimated by dadi and written to the main output files as part of this pipeline

mu refers to the mutation rate (specific to your system or a best guess)

We mapped all sequencing reads for 30 sequence data to the  reference  genome,  
then  used  samtools  to  call  SNPs.Heterozygosity across the assembled sections
of theO. pumiliogenome is H=0.0016. 
Assuming a mutation rate of 10-9, this would yield an estimate of Ne=400,000. 

'''

mu=10e-9
L=3306892825 #sites analyses by angsd first round

Nref=model_sel.theta/(4*mu*L)

#s,nu1,nu2,T,m12,m21 = params
pop1_Ne_aftersplit=model_sel.s*Nref 
pop2_Ne_aftersplit=(1-model_sel.s)*Nref 
pop1_Ne=model_sel.nu1 * Nref 
pop2_Ne=model_sel.nu2 * Nref 

'''
model = moments.ModelPlot.generate_model(demo_model_ex, list(model_sel), ns)
moments.ModelPlot.plot_model(model,
                             pop_labels=[args.pop1, args.pop2],
                             draw_scale=True,
                             reverse_timeline=True)

plt.savefig(args.pop1+args.pop2+"bestIMmodel.png")
'''
