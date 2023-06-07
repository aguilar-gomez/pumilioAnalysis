#!/usr/bin/env python3
# coding: utf-8
# Author : Diana Aguilar
import pandas as pd
import argparse

#Arguments
parser = argparse.ArgumentParser()
parser.add_argument("fixed_af", type=float, help="frequency of allele to be considered fixed, for low coverage data it is 
not recommended to use 1")
parser.add_argument("pop1", type=str, help="name of population 1 (prefix of maf file)")
args = parser.parse_args()

#Read files
outfile = args.pop1 + "fil.mafs"
fixed_af = args.fixed_af
pop1=pd.read_csv(args.pop1 + ".mafs",sep="\t")
polysites=pd.read_csv("pumilio_version11.pos",sep="\t",names=["chromo","position"])

#Merge with confirmed polymorphic sites
pop1poly= pop1.merge(polysites, on=["chromo","position"], how='left',indicator='PolyConfirmed')
pop1poly['PolyBinary'] = pop1poly['PolyConfirmed'].eq('both')

#Keep if confirmed polymorphic or fixed
KeepSites=pop1poly[((pop1poly.knownEM<1-fixed_af) | (pop1poly.knownEM>fixed_af)) | (pop1poly.PolyBinary==True)]
KeepSites.to_csv(outfile,columns=["chromo","position","major","minor","ref","anc","knownEM","nInd"],sep="\t",index=False)
