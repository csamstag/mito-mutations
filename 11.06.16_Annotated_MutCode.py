import numpy

MasterListFile = open("7.19.16_MasterMutListFix.txt", "r")
'''
format:
 [0]refpos [1]refnuc [2] PhyloP [3]gene [4]TMutPred [5]CMutPred
[6]GMutPred [7]AMutPred
For column [3], noncoding regions are represented by "."
For columns [4]-[7], all noncoding regions have only '.'
    all synonymous mutations are given a mutpre==0
    all NS mutations have mutpreds > 0
    Stop codon mutations have a mutpred ==1
    
    NS/S ratio = ratio nonzero:zero mutpreds
    
mutdict (and the nucdicts) will simply have the same format, whereby the key is
the refpos and the value is the remainder of the fields, in list form
refpos : [[0]refnuc [1]phylop [2]gene [3]tmutpred [4]cmutpred [5]gmutpred [6]amutpred
columns [3] - [6] show the mutpred if the base is mutated to the corresponding
    nucleotide in protein coding genes].
    
The following loop will create two dictionary types, both of the same format
1) One dictionary with all mutations, and their associated phylops, mutpreds
2) four dictionaries, one for each nucleotide type (A,C,G,T) according to the ref genome
in either case;
the dictionary key = field[0], aka the refnuc position

'''

mutdict={}
TDict = {}
CDict = {}
ADict={}
GDict={}
for line in MasterListFile:
    x = line.splitlines()
    for mut in x:
        mutline=mut.split()
        mutdict[mutline[0]] = [mutline[1], mutline[2], mutline[3], mutline[4],
        mutline[5], mutline[6], mutline[7]]
        if mutline[1] == "A":
            ADict[mutline[0]] = [mutline[1], mutline[2], mutline[3], mutline[4],
            mutline[5], mutline[6], mutline[7]]
        if mutline[1] == "C":
            CDict[mutline[0]] = [mutline[1], mutline[2], mutline[3], mutline[4],
            mutline[5], mutline[6], mutline[7]]
        if mutline[1] == "G":
            GDict[mutline[0]] = [mutline[1], mutline[2], mutline[3], mutline[4],
            mutline[5], mutline[6], mutline[7]]
        if mutline[1] == "T":
            TDict[mutline[0]] = [mutline[1], mutline[2], mutline[3], mutline[4],
            mutline[5], mutline[6], mutline[7]]
MasterListFile.close()
seqFile = open("7.16.18_NM13_dm6.txt", "r")
SimulOutputFile=open("11.10.16_TestOutput.txt", 'w')

'''
seqFile format:
    [0]refGenome [1]refNuc [2]refPos [3]seqDepth [4]totmuts [5]Tmuts [6]Cmuts
    [7]Gmuts [8]Amuts [9]Ins [10]Dels [11]Und
The following declarations are to create lists, one for each nucleotide type,
for (1) the positions sequenced, (2) the sequencing depth at the corresponding
position, and (3) the fraction of all of sequenced bases of that type that that
position covers

Essentially, we want to control for uneven sequencing depth across the genome
Basically, once we account for the mutation bias of the polymerase, we want to model
the effects of mutation occurring completely randomly, without respect to the position
in the genome or its consequence.  However, the probability of detecting a mutation
through sequencing is based on both the probability of a mutation occurring at this
site (should be completely random according to our simulations) X the probability
of detecting a mutation in this site (aka the seq depth, relative to other sites).
That means that our simulations must perform random draws from the population of
nucleotides, weighted by the proportion of nucleotides that were detected at that
site.  So we create lists that correspond with each nucleotide of that type sequenced,
populate another list with the corresponding seq depth, then iterate through that
list to create a third list that shows the proportion of that nucleotide type
that each individual site represents

TLDR: If two sites are equally mutable, but site A is sequenced to a 1000x depth and 
site B sequenced to 10x, we would be 100 times more likely to observe a mutation
in site A despite the fact that they are equally mutable.  We need to account for
this bias in observation by performing *weighted* random draws for our
mutagenesis simulations
'''

#1: positions of each nucleotide type sequenced
AList=[]
CList=[]
GList=[]
TList=[]
#2L the corresponding sequencing depth at that site
ASeqDepthList=[]
CSeqDepthList=[]
GSeqDepthList=[]
TSeqDepthList=[]
#3: the fraction of all sequenced bases of that nucleotide type
AFracList = []
CFracList = []
GFracList = []
TFracList = []


obsTotMuts = 0
obsAMuts = 0
obsCMuts= 0
obsGMuts = 0
obsTmuts = 0
mutslist = []
mutSpectrumList = []

'''
parses the sequencing file, a .mutpos file
populates list screated above
'''
for line in seqFile:
    v = line.splitlines()
    for x in v:
        mutline = x.split()
        refpos=mutline[2]
        seqdepth = mutline[3]
        if mutline[1] == "A":
            AList.append(int(refpos))
            ASeqDepthList.append(int(seqdepth))
        elif mutline[1] == "C":
            CList.append(int(refpos))
            CSeqDepthList.append(int(seqdepth))
        elif mutline[1] == "G":
            GList.append(int(refpos))
            GSeqDepthList.append(int(seqdepth))
        elif mutline[1] == "T":
            TList.append(int(refpos))
            TSeqDepthList.append(int(seqdepth))
        '''
        if the number of total mutations at a site (mutline[4]) > 0:
            add this site to the list of observed mutations.
            This list takes mutations in the format RefBase+MutatedBase+NucPos
            eg. "CA368".  This is important because the same function will later
            check the consequence of each of these observed mutations, as well
            as the consequece of each of the simulated mutations
        '''
        if int(mutline[4])> 0:
            obsTotMuts += int(mutline[4])
            for w in range(0, int(mutline[5])):
                        mutslist.append([str(mutline[1]), "T", str(mutline[2])])
                        #mutslist.append(str(mutline[2])+"T")
                        mutSpectrumList.append(str(mutline[1])+"T")
            for w in range(0, int(mutline[6])):
                        mutslist.append([str(mutline[1]), "C", str(mutline[2])])
                        #mutslist.append(str(mutline[2])+"C")
                        mutSpectrumList.append(str(mutline[1])+"C")
            for w in range(0, int(mutline[7])):
                        mutslist.append([str(mutline[1]), "G", str(mutline[2])])
                        #mutslist.append(str(mutline[2])+"G")
                        mutSpectrumList.append(str(mutline[1])+"G")
            for w in range(0, int(mutline[8])):
                        mutslist.append([str(mutline[1]), "A", str(mutline[2])])
                        #mutslist.append(str(mutline[2])+"A")
                        mutSpectrumList.append(str(mutline[1])+"A")
seqFile.close()
print len(mutslist)

allMutTypes = ["AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"]
mutSpecDict = {}
for item in allMutTypes:
    mutSpecDict[item] = 0
for item in mutSpectrumList:
    if str(item) not in mutSpecDict:
        print item
    elif str(item) in mutSpecDict:
        mutSpecDict[item] += 1
print mutSpecDict   
for item in allMutTypes:
    print str(item) + '\t' + str(mutSpecDict[item])


'''
this section will go through all the bases sequenced and create lists of all the
nucleotides and their sequencing depth.  Then it will individually sum up all A, 
C, G, and T sites sequenced to get the total number of each nucleotide type sequenced.
Finally, it will make lists of the fraction of total sites of that nucleotide
type at each position in the genome.  This list can be used to correspond with
the nucleotide position list, to be used as a probability for the random choice
function (to account for sequencing depth).  Can these be zipped together???
'''
TotalASeqDepth = sum(ASeqDepthList)
TotalCSeqDepth = sum(CSeqDepthList)
TotalGSeqDepth = sum(GSeqDepthList)
TotalTSeqDepth = sum(TSeqDepthList)


for q in ASeqDepthList:
    AFracList.append(float(q)/TotalASeqDepth)
for q in CSeqDepthList:
    CFracList.append(float(q)/TotalCSeqDepth)
for q in GSeqDepthList:
    GFracList.append(float(q)/TotalGSeqDepth)
for q in TSeqDepthList:
    TFracList.append(float(q)/TotalTSeqDepth)
'''
the parsemuts function takes a list of mutations formatted as outlined above:
    if site 368, a C in the ref genome, is mutated to a A, list has value "CA368"
    this consults a master list, checking
        - the phyloP of the site mutated
        - whether the site is within a protein coding gene
        - if it encodes a protein, is the site synonymous or not (eg, does it
           have a mutpred assigned to 0 or some nonzero value)
        - if mutpred is nonzero, this is a NS site, and it records the mutpred
          of that NS mutation
        - if mutpred == 1, then it also records that this is a stop codon altering
          mutation (the only mutations assigned a mutpred of 1 are those that
          cause premature stop codons or those that affect the normal stop codon)

this function will be used to analzye both the mutations obtained from sequencing
as well as lists of mutations generated from our simulations
'''
def parsemuts(mutsToParse):
    parsedPhyloPList = []
    parsedMutPredList = []
    NCPhyloP = []
    CodingPhyloP = []
    NCmuts =0
    codingMuts = 0
    SynMuts = 0
    NSMuts = 0
    for mut in mutsToParse:
        '''    
         refbase is the reference nucleotide (GCAT)
         mutbase is the nucleotide it is mutated to (GCAT)
         mutpos is the gene that this mutation is located within
            if in protein sequence, it will show the protein
            if it is within a trna, rrna, or noncoding, mutpos will read "."
     
        '''     
        
        refbase=mut[0]
        mutbase=mut[1]
        mutpos = mut[2]
        #print mutpos
        parsedPhyloPList.append(float(mutdict[mutpos][1]))
        #print mutpos, mutdict[mutpos]
        
        '''
        to check the phylops only in rrna sequences
        if int(mut[2]) >12734:
            if int(mut[2]) < 14916:
                parsedPhyloPList.append(float(mutdict[mutpos][1]))
        
        '''
        
        if mutdict[mutpos][2] == '.':
            NCmuts+=1
            NCPhyloP.append(float(mutdict[mutpos][1]))
        elif mutdict[mutpos][2] != '.':
            #parsedPhyloPList.append(float(mutdict[mutpos][1]))
            CodingPhyloP.append(float(mutdict[mutpos][1]))
            codingMuts+=1
            parsedPhyloPList.append(float(mutdict[mutpos][1]))
            if mutbase == "A":
                if mutdict[mutpos][6] == '0':
                    SynMuts+=1
                elif mutdict[mutpos][6] != '0':
                    NSMuts +=1
                    parsedMutPredList.append(float(mutdict[mutpos][6]))
            if mutbase == "C":                 
                if mutdict[mutpos][4] == '0':
                    SynMuts+=1
                elif mutdict[mutpos][4] != '0':
                    NSMuts +=1
                    parsedMutPredList.append(float(mutdict[mutpos][4]))
            if mutbase == "G":
                if mutdict[mutpos][5] == '0':
                    SynMuts+=1
                elif mutdict[mutpos][5] != '0':
                    NSMuts +=1
                    parsedMutPredList.append(float(mutdict[mutpos][5]))
            if mutbase == "T":
                if mutdict[mutpos][3] == '0':
                    SynMuts+=1
                elif mutdict[mutpos][3] != '0':
                    NSMuts +=1
                    parsedMutPredList.append(float(mutdict[mutpos][3]))
    #print len(parsedPhyloPList)
    meanPhyloP = numpy.mean(parsedPhyloPList)
    #print "# coding muts " + str(len(CodingPhyloP))
    #print "# noncoding muts " + str(len(NCPhyloP))
    #print "total muts " + str((len(CodingPhyloP))+(len(NCPhyloP))) 
    meanCodingPhyloP = numpy.mean(CodingPhyloP)
    meanNCPhyloP = numpy.mean(NCPhyloP)
    #print meanPhyloP
    meanMutPred = numpy.mean(parsedMutPredList)
    #print "mutpred :" + str(meanMutPred)
    if float(SynMuts)==0:
        NSSRatio = float(NSMuts)
    elif float(SynMuts) != 0:
        NSSRatio = float(NSMuts)/float(SynMuts)
    '''
    #print parsedPhyloPList
    #print parsedMutPredList
    print "total muts " + str(len(mutsToParse))
    print "NC Muts " + str(NCmuts)
    print "Coding Muts " + str(codingMuts)
    print "Syn Muts " + str(SynMuts)
    print "NS Muts " + str(NSMuts)
    print "NS/S Ratio " + str(NSSRatio)
    '''
    return meanPhyloP, meanCodingPhyloP, meanNCPhyloP, meanMutPred, NSMuts, SynMuts, NSSRatio
    #return meanMutPred
    #return NSMuts
    #return SynMuts
    #return NSSRatio    


#mutdict:
    #refpos = list: [0]refnuc [1] phyloP [2] gene [3] T [4] C [5] G [6]A
#for x in range(1,50):
#    print mutdict[str(x)]
mutpredlist = []
phylolist = []
codinglist = []
nsslist = []
nclist=[]
redundMutSpecList =[]
SynMutList =[]
scount=0

'''
we call a random choice function for each mutation type (eg A->T mutations), equal
to the number of mutations of each type seen through sequencing
the following function generates the list of mutations 
random.choice() takes the following three parameters:
        -the list from which it should draw (in this case, the list of all sites
            sequenced from that particular nucleotide
        -the number of mutations detected of that mutation type
        -the proportion of total sequenced reads each site makes up (eg, the prob
         we would detect a mutation there given unequal sequencing depth--we are
         2x more likely to have observed a mutation in a site that we sequenced
         1000 molecules than a site that we sequenced 500 molecules.  The list of
         proportions used here corresponds to the same indices used in the list
         of sites 

once we generate a proportional number of random mutations to each type observed
for all 12 mutation types, we generate another list for each that adds the mutation
type to each mutation generated (eg, takes the randomly generated nucleotide
position 368 from the simulCAList, and populates a new list with the value "CA368".
Finally, we concatenate a list of ALL simulated mutations in the above form.
This final list should include the same exact number of each type of mutation
observed in the sequencing data. 
'''
def simulmut():

    simulACList= (numpy.random.choice(AList, mutSpecDict["AC"], p=AFracList))
    simulAGList= (numpy.random.choice(AList, mutSpecDict["AG"], p=AFracList))
    simulATList= (numpy.random.choice(AList, mutSpecDict["AT"], p=AFracList))
    simulCAList= (numpy.random.choice(CList, mutSpecDict["CA"], p=CFracList))
    simulCGList= (numpy.random.choice(CList, mutSpecDict["CG"], p=CFracList))
    simulCTList= (numpy.random.choice(CList, mutSpecDict["CT"], p=CFracList))
    simulGAList= (numpy.random.choice(GList, mutSpecDict["GA"], p=GFracList))
    simulGCList= (numpy.random.choice(GList, mutSpecDict["GC"], p=GFracList))
    simulGTList= (numpy.random.choice(GList, mutSpecDict["GT"], p=GFracList))
    simulTAList= (numpy.random.choice(TList, mutSpecDict["TA"], p=TFracList))
    simulTCList= (numpy.random.choice(TList, mutSpecDict["TC"], p=TFracList))
    simulTGList= (numpy.random.choice(TList, mutSpecDict["TG"], p=TFracList))
    
    ACListForm= []
    AGListForm= []
    ATListForm= []
    CAListForm= []
    CGListForm= []
    CTListForm= []
    GAListForm= []
    GCListForm= []
    GTListForm= []
    TAListForm= []
    TCListForm= []
    TGListForm= []
    
    for item in simulACList:
        ACListForm.append(["A", "C", str(item)])
    for item in simulAGList:
        AGListForm.append(["A", "G", str(item)])
    for item in simulATList:
        ATListForm.append(["A", "T", str(item)])
    for item in simulCAList:
        CAListForm.append(["C", "A", str(item)])
    for item in simulCGList:
        CGListForm.append(["C", "G", str(item)])
    for item in simulCTList:
        CTListForm.append(["C", "T", str(item)])
    for item in simulGAList:
        GAListForm.append(["G", "A", str(item)])
    for item in simulGCList:
        GCListForm.append(["G", "C", str(item)])
    for item in simulGTList:
        GTListForm.append(["G", "T", str(item)])
    for item in simulTAList:
        TAListForm.append(["T", "A", str(item)])
    for item in simulTCList:
        TCListForm.append(["T", "C", str(item)])
    for item in simulTGList:
        TGListForm.append(["T", "G", str(item)])
    
    simultotmuts = []

    for item in ACListForm:
        simultotmuts.append(item)
    for item in AGListForm:
        simultotmuts.append(item)
    for item in ATListForm:
        simultotmuts.append(item)
    for item in CAListForm:
        simultotmuts.append(item)
    for item in CGListForm:
        simultotmuts.append(item)
    for item in CTListForm:
        simultotmuts.append(item)
    for item in GAListForm:
        simultotmuts.append(item)
    for item in GCListForm:
        simultotmuts.append(item)
    for item in GTListForm:
        simultotmuts.append(item)
    for item in TAListForm:
        simultotmuts.append(item)
    for item in TCListForm:
        simultotmuts.append(item)
    for item in TGListForm:
        simultotmuts.append(item)    
    #print "# Coding Mutations: " + str(len(simCodPhyloPList))
    #print "# Noncoding Mutations " + str(len(simNCPhyloPList))
    return simultotmuts

simPhyloList = []
simCodPhyloPList=[]
simNCPhyloPList =[]
simMutPredList = []
simNSSRatioList = []
PhyloAboveCount =0
CodPhyloAboveCount=0
NCPhyloAboveCount=0

NSSRatioAboveCount=0
MutPredAboveCount =0


#print mutSpecDict
print "OBSERVED"
obsPhyloP, obsCodingPhyloP, obsNCPhyloP, obsMutPred, obsNS, obsSyn, obsNSSRatio = parsemuts(mutslist)
print "observed PhyloP: " + str(obsPhyloP)
print "observed Coding PhyloP: " + str(obsCodingPhyloP)
print "observed NC PhyloP: " + str(obsNCPhyloP)
print "observed MutPred: " + str(obsMutPred)
print "observed #NS Muts: " + str(obsNS)
print "observed #S Muts: " + str(obsSyn)
print "observed NS:S Ratio: " + str(obsNSSRatio)

print "\n SIMULATED "
#the next line sets the number of simulations you want to run, can be adjusted    
totalsims=100000

'''
this for loop will repeatedly call our simulation function up to the number of
times you want to simulate random mutagenesis.  Each time, it runs a full simulation,
stores the PhyloP,  NS/S ratio, and mutpred obtained from the simulation run and
stores these in their own lists.
It also compares these values with those obtained from sequencing, and gives a quick
print readout as to how many values from simulations were more extreme than those
observed from sequencing.

'''

for t in range(0,totalsims):
    #this line gives a status update to show how many simulations have been run
    if (t % 1000 == 0):
        print (t)
    simulPhyloP, simulCodingPhyloP, simulNCPhyloP, simulMutPred, simulNS, simulSyn, simulNSSRatio = parsemuts(simulmut())
    simPhyloList.append(simulPhyloP)
    simCodPhyloPList.append(simulCodingPhyloP)
    simNCPhyloPList.append(simulNCPhyloP)
    simMutPredList.append(simulMutPred)
    simNSSRatioList.append(simulNSSRatio)
    '''
    the subsequent lines don't affect the simulations being run, but they allow a 
    quick readout of how many simulations are more extreme than the values that you
    observed (necessary for the empirical p-value)

    '''
    if float(simulPhyloP) >= float(obsPhyloP):
        PhyloAboveCount +=1
    if float(simulCodingPhyloP) >= float(obsCodingPhyloP):
        CodPhyloAboveCount +=1
    if float(simulNCPhyloP) >= float(obsNCPhyloP):
        NCPhyloAboveCount +=1
    if float(simulNSSRatio) >= float(obsNSSRatio):
        NSSRatioAboveCount +=1
    if float(simulMutPred) >= float(obsMutPred):
        MutPredAboveCount +=1
    
print "Simulated PhyloPs > Observed: " + str(PhyloAboveCount) + "/ " + str(len(simPhyloList))
print "Simulated Coding PhyloPs > Observed: " + str(CodPhyloAboveCount) + "/ " + str(len(simPhyloList))
print "Simulated NC PhyloPs > Observed: " + str(NCPhyloAboveCount) + "/ " + str(len(simPhyloList))
print "Simulated NSSRatio > Observed: " + str(NSSRatioAboveCount) + "/ " + str(len(simPhyloList))
print "Simulated MutPred > Observed: " + str(MutPredAboveCount) + "/ " + str(len(simPhyloList))


'''
The output file will have the values for each of the parameters outliend above,
first those observed in the sequencing data, and then a column containing
each of the average values from a simulation run.  These can be copied into 
statistical software for further analases and to generate histograms
'''
SimulOutputFile.write("Observed PhyloP \t Observed Coding PhyloP \t Observed NC PhyloP" +
    "\t Observed NS:S Ratio \t Observed MutPred \n" +
    str(obsPhyloP) + "\t" + str(obsCodingPhyloP) + "\t" + str(obsNCPhyloP) + "\t" +
    str(obsNSSRatio) + "\t" + str(obsMutPred) + '\n\n\n'+
    "Simulated PhyloP \t Simulated Coding PhyloP \t Simulated NC PhyloP" +
    "\t Simulated NS:S Ratio \t Simulated MutPred \n"
    )
   
for p in range(0, totalsims):
    
    SimulOutputFile.write(
    str(simPhyloList[p]) + '\t' +
    str(simCodPhyloPList[p]) + '\t' +
    str(simNCPhyloPList[p]) + '\t' +
    str(simNSSRatioList[p]) + '\t' +
    str(simMutPredList[p]) + '\n'
    )    

SimulOutputFile.close()



