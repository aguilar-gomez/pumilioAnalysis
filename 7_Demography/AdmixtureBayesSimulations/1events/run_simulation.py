import msprime
from math import exp, log
demography = msprime.Demography()
outgroupsize =  -(7000/2)/log(1-.180)
demography.add_population(name="outgroup", initial_size=outgroupsize)   #edit this.
PopulationsToSizes = [  
    ["n1", -(1000/2)/log(1-0.026631869540764526)],   #done 
    ["n2", -(1000/2)/log(1-0.05961244240609137)],   #done
    ["n3", -(1000/2)/log(1-0.017408711906661142)],  #done   
    ["n4", -(1000/2)/log(1-0.04725148202144411)],   #done
    ["n5", -(2000/2)/log(1-0.04890169022573293)],  #done
    ["n6", -(1000/2)/log(1-0.12976622201367452)],  
    ["n7", -(4000/2)/log(1-0.06603735046700507)],   
    ["n8", -(1000/2)/log(1-0.027059443095721534)],    #done 
    ["n9", -(1000/2)/log(1- 0.06275448341800476)],    #done 
    ["R", -(1000/2)/log(1-.115)],   
    ["a1", -(500/2)/log(1-0.16253357576349323)],    #done
    ["a2", -(500/2)/log(1-0.14737563090666114)],    #done
    ["ancestral", -(2000/2)/log(1-.115)],    
    ["PD", -(2000/2)/log(1-0.03153178318481301)],    #done
    ["TB", -(1000/2)/log(1-0.14119398564083704)],   #done 
    ["HP", -(500/2)/log(1-0.0759845319093546)] ,  #done
    ["CL", -(2000/2)/log(1-0.06521784568445044)], #done
    ["CM", -(1000/2)/log(1-0.12612118026344143)],   #done 
    ["PO", -(5000/2)/log(1-0.0886506846458096)],    #done
    ["RA", -(3000/2)/log(1-0.024614890126696362)] ,    #done
    ["SC", -(2000/2)/log(1-0.23619542012690356)], #done
    ["DB", -(1000/2)/log(1-0.00643706618149798)],   ##done! 
    ["SK", -(1000/2)/log(1-1.2284299388791048/100000)],     #done
]

for i in PopulationsToSizes: print(i)
####################################################################
for i in PopulationsToSizes:
    demography.add_population(name=i[0], initial_size=i[1])

DivPairs = []
AbsoluteTimes = []
AdmixProp = []
InitialList = []

def getalphabetical(val):
    for i in DivPairs:
        if val in i and i[0] <  i[1]:
            return(i[0])
        if val in i and i[1] <  i[0]:
            return(i[1])
    return(val)

def getadmixtureproportions(val):
    for i in AdmixProp:
        if val == i[0]:
            return(i[1])
    return 1.0

def getabsolutetime(val):
    for i in AbsoluteTimes:
        if val == i[0]:
            return(i[1])
    return 0.0

def addrelevantedges(absolutetime, child, ancestral):
    childreturn = getalphabetical(child)
    for i in PopulationsToSizes:
        if child == i[0]:
            timee = 1-exp(-(absolutetime - getabsolutetime(childreturn))/(2*i[1]))
            break
    admixture = getadmixtureproportions(child)
    AbsoluteTimes.append([ancestral, absolutetime])
    return(childreturn + " " + ancestral + " " + str(timee) + " "  + str(admixture))

def addadmixture(absolutetime, child, ancestral1, ancestral2, prop1, prop2):
    childreturn = ancestral1
    if ancestral2 < ancestral1:
        childreturn = ancestral2
    for i in PopulationsToSizes:
        if child == i[0]:
            timee = 1-exp(-(absolutetime - getabsolutetime(child))/(2*i[1]))
            break
    admixture = getadmixtureproportions(child)
    AbsoluteTimes.append([childreturn, absolutetime])
    AdmixProp.append([ancestral1, prop1])
    AdmixProp.append([ancestral2, prop2])
    DivPairs.append([ancestral1, ancestral2])
    return(child + " " + childreturn + " " + str(timee) + " "  + str(admixture))

def causeevent(timee, der1, der2, anc):
    demography.add_population_split(time=timee, derived=[der1, der2], ancestral=anc)
    InitialList.append(addrelevantedges(timee, der1, anc))
    InitialList.append(addrelevantedges(timee, der2, anc))

####################################################################

demography.add_admixture(time=500, derived="HP", ancestral=["a1", "a2"], proportions=[0.5, 0.5])       #PLACE3

InitialList.append(addadmixture(500, "HP",  "a1", "a2", 0.5, 0.5))

causeevent(1000, "SK", "DB", "n1")
causeevent(1000.01, "TB", "a1", "n4")
causeevent(1000.001, "a2", "CM", "n2")
causeevent(2000, "PD", "n4", "n7")
causeevent(2000.001, "SC", "n1", "n3")
causeevent(2000.0001, "CL", "n2", "n5")
causeevent(3000, "n3", "RA", "n6")
causeevent(4000, "n5", "n6", "n8")
causeevent(5000, "n8", "PO", "n9")
causeevent(6000, "n9", "n7", "R")

demography.add_population_split(time=7000, derived=["R", "outgroup"], ancestral="ancestral")

demography.sort_events()
ts = msprime.sim_ancestry( {
    "PD" : 3, 
    "TB" : 35, 
    "HP" : 30, 
    "CL" : 15 , 
    "CM" :  63, 
    "PO" :  34 , 
    "RA" : 7 , 
    "SC" : 26 , 
    "DB" : 104 , 
    "SK" :  30  ,
    "outgroup" :  4 }, 
      demography=demography, sequence_length=50000000, recombination_rate = 1e-8)      #PLACE3 ENDS
######################################################################
mts = msprime.sim_mutations(ts, rate =  2* (10**(-9)) ) # 376292 snps about

with open("TemporaryFiles/Data.vcf", "w") as vcf_file:
    mts.write_vcf(vcf_file)

with open('TemporaryFiles/TrueTree.txt', 'w') as f:
    for line in InitialList:
        f.write(line + "\n")
        print(line)
