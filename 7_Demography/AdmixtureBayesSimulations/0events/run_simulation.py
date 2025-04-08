import msprime
from math import exp, log
demography = msprime.Demography()
outgroupsize =  -(8000/2)/log(1-.180)
demography.add_population(name="outgroup", initial_size=outgroupsize)   #edit this.
PopulationsToSizes = [   #PLACE2
    ["n1", -(500/2)/log(1-.033)],    #FINAL
    ["n2", -(1000/2)/log(1-.026)],    #FINAL
    ["n3", -(500/2)/log(1-.064)],    #FINAL
    ["n4", -(1000/2)/log(1-.0183)],    #FINAL
    ["n5", -(1000/2)/log(1-.139)],    #FINAL
    ["n6", -(1000/2)/log(1-.0182)],    #FINAL
    ["n7", -(1000/2)/log(1-.046)],    #FINAL
    ["n8", -(1000/2)/log(1-.0246)],    #FINAL
    ["R", -(1000/2)/log(1-.115)],    #FINAL
    ["ancestral", -(2000/2)/log(1-.115)],    #FINAL
    ["PD", -(6500/2)/log(1-.044)],    #FINAL
    ["TB", -(6500/2)/log(1-.176)],    #FINAL
    ["HP", -(6000/2)/log(1-.141)] ,  #FINAL
    ["CL", -(4500/2)/log(1-.0767)],
    ["CM", -(4500/2)/log(1-.176)],    #FINAL
    ["PO", -(4000/2)/log(1-.102)],    #FINAL
    ["RA", -(3000/2)/log(1-.024)] ,    #FINAL
    ["SC", -(2000/2)/log(1-.2367)],
    ["DB", -(1000/2)/log(1-.0064)],    #FINAL
    ["SK", -(1000/2)/log(1-.12/100000000)],    #FINAL
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

causeevent(1000, "SK", "DB", "n2")
causeevent(2000, "SC", "n2", "n4")
causeevent(3000, "n4", "RA", "n5")
causeevent(4000, "n5", "PO", "n6")
causeevent(4500, "CL", "CM", "n1")
causeevent(5000, "n1", "n6", "n7")
causeevent(6000, "n7", "HP", "n8")
causeevent(6500, "PD", "TB", "n3")
causeevent(7000, "n3", "n8", "R")

demography.add_population_split(time=8000, derived=["R", "outgroup"], ancestral="ancestral")

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