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
import sys
from dadi import Numerics, PhiManip, Integration, Spectrum


parser = argparse.ArgumentParser()
parser.add_argument("--n_sim", type=int, help="number of simulations")
parser.add_argument("--pop1", type=str, help="name of population 1")
parser.add_argument("--pop2", type=str, help="name of population 2")

args = parser.parse_args()

dataset = args.pop1 + args.pop2
inputfs = dadi.Spectrum.from_file(dataset+".2dfs")

#Calculate Fst
fst=inputfs.Fst()

print("Analizying",args.pop1, args.pop2, "with Fst:", fst, "number of optimization:",args.n_sim)

#Retrive the sample sizes from the data
ns = inputfs.sample_sizes

############### Specifying a model ########################
def IM(params, ns, pts):
    s, nu1, nu2, T, m12, m21 = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t : s * (nu1/s) ** (t/T)
    nu2_func = lambda t : (1-s) * (nu2/(1-s)) ** (t/T)

    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12 = m12, m21 = m21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

param_names=["s","n1","n2","T","m12",""m21]

#Grid sizes
# Define the grid points based on the sample size.
# For smaller data (largest sample size is about <100) [ns+20, ns+30, ns+40] is a good starting point.
# For larger data (largest sample size is about >=100) or for heavily down projected data [ns+100, ns+110, ns+120] is a good starting point.
pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]


# Wrap the demographic model in a function that utilizes grid points
func=split# set the function
# Make extrapolation function:
func_ex = dadi.Numerics.make_extrap_log_func(func)


# Define starting parameters
#First two parameters effective pop size, 
#Times are given in units of 2Nref (reference population size) generations
#s,n1,n2,T,m,m = params

params = [.5,1,1,1,1,1]

# Define boundaries of optimization.
# It is a good idea to have boundaries to avoid optimization
# from trying parameter sets that are time consuming without
# nessicarily being correct.
# If optimization infers parameters very close to the boundaries, we should increase them.
lower_bounds = [1e-3,1e-3, 1e-3, 1e-3,1e-5,1e-5]
upper_bounds = [1,200, 200, 300,100,100]
maxiter=100

# Perturb parameters
# Optimizers dadi uses are mostly deterministic
# so we will want to randomize parameters for each optimization.
p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,
                              lower_bound=lower_bounds)

# Run optimization
# At the end of the optimization we will get the
# optimal parameters and log-likelihood.
# Verbose is how many evaluations are done before pàrameter is printed
popt = dadi.Inference.optimize_log(p0, inputfs, func_ex, pts_l,
                                    lower_bound=lower_bounds,
                                    upper_bound=upper_bounds,
                                    maxiter=maxiter, verbose=0)

# Calculate the synonymous theta
model_fs = func_ex(popt, ns, pts_l)
theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, inputfs)

# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model_fs, inputfs)

# also get the LL of the data to itself (best possible ll)
ll_data=dadi.Inference.ll_multinom(inputfs, inputfs)


'''
The Nref is calculated by the following equation:
Theta = 4 * Nref * mu * L
Nref=Theta/(4*mu*L)

L is the total length of DNA sequence (in bp) that was analyzed in order to obtain the SNP data

We mapped all sequencing reads for 30 sequence data to the  reference  genome,  
then  used  samtools  to  call  SNPs.Heterozygosity across the assembled sections
of the O. pumiliogenome is H=0.0016. 
Assuming a mutation rate of 10-9, this would yield an estimate of Ne=400,000. 

'''

mu=10e-9
L=3306892825 #sites analyses by angsd first round

Nref=theta0/(4*mu*L)

#generation time is 2 years
g=2
#Scale parameters to be diploid
N01=popt[0]*Nref
N02=(1-popt[0])*Nref
N1=popt[1]*Nref
N2=popt[2]*Nref
#Time is in units T=2Nref generations 
split_time=popt[3]*2*Nref*g
#Migration rate is obtained by dividing by 2Nref
#mij=Mij/(2*Nref) 
#mij=fraction of individuals in each generation in population i who are new migrants from population j.
m12=popt[4]/(2*Nref)
m21=popt[5]/(2*Nref)


scaled_param_names=("Nref","N01","N02,"N1","N2","split_time_years","m12","m21")
scaled_popt=(Nref,N01,N02,N1,N2,split_time,m12,m21)

#### Write output                          
outputFile=open(dataset+".dadi.inference.run."+str(args.n_sim)+".output","w")
# get all param names:
param_names_str='\t'.join(str(x) for x in param_names)
scaled_param_names_str='\t'.join(str(x) for x in scaled_param_names)
header=param_names_str+"\t"+scaled_param_names_str+"\ttheta\tLL\tLL_data\tmu\tL\tmaxiter\trunNumber\tinitialParameters\tupper_bound\tlower_bound" # add additional parameters theta, log-likelihood, model name, run number and rundate
popt_str='\t'.join(str(x) for x in popt) # get opt'd parameters as a tab-delim string
scaled_popt_str='\t'.join(str(x) for x in scaled_popt)
# joint together all the output fields, tab-separated:
output=[popt_str,scaled_popt_str,theta0,ll_model,ll_data,mu,L,maxiter,args.n_sim,p0,upper_bounds,lower_bounds] # put all the output terms together
output='\t'.join(str(x) for x in output) # write out all the output fields
# this should result in a 2 row table that could be input into R / concatenated with other runs
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()

############### Output SFS ########################
print('Writing out SFS **************************************************')                                   

outputSFS=dataset+".dadi.inference.run."+str(args.n_sim)+".sfs"
model_fs.to_file(outputSFS)

############### Output plot ########################
print('Making plots **************************************************')                                   

fig=plt.figure(1)
outputFigure=dataset+".dadi.inference.run."+str(args.n_sim)+".png"
dadi.Plotting.plot_2d_comp_multinom(model_fs, inputfs)
plt.savefig(outputFigure)

###### exit #######
sys.exit()
