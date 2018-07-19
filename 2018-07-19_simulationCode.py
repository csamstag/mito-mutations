'''
takes three inputs: the name of the input file to be run, the name of your output
file, and the number of simulations you want run.

Will eventually have two output files:
    one file with the averages of each of n simultations
    one file with the comparisons with the observed metrics
    
example input 
python 2018-05-21_simulationCode.py InputFile.txt OutputFileNames 10000
'''

import sys, numpy, pprint

MasterListFile = open("7.19.16_MasterMutListFix.txt", "r")
'''
First we read in a file with the consequences for every possible mutation in the
mitochondrial genome--PhyloP for every nucleotide in the genome, and 
NS/S and MutPred score if it is within a protein coding sequence.
format:
 [0]refpos [1]refnuc [2] PhyloP [3]gene [4]TMutPred [5]CMutPred
[6]GMutPred [7]AMutPred
For column [3], noncoding regions are represented by "."
For columns [4]-[7], all noncoding regions have only '.'
    all synonymous mutations are given a mutpred=0 (not factored into calculations
        of average mutpred score)
    all NS mutations have mutpreds > 0
    Stop codon mutations have a mutpred = 1
    
    NS/S ratio = ratio nonzero:zero mutpreds
    
mutdict (and the nucdicts) will simply have the same format, whereby the key is
the refpos
refpos: [0]refnuc [1]phylop [2]gene [3]tmutpred [4]cmutpred [5]gmutpred [6]amutpred
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
#seqFile = open("/Users/student/Desktop/03.05.18_Polg_Redo/polgAnalyzed/SNVs/70Clonal/combined/Comb_NM_SNV_70_NoDLoop.txt", "r")
#SimulOutputFile=open("05.05.18_Comb_NM_NoDLoop_simuls", 'w')
#Next we read in the input and output files based on the command line arguments
seqFile = open(sys.argv[1], "r")
simulOutputFile=open(sys.argv[2]+"simuls.txt", 'w')
summaryOutputFile=open(sys.argv[2]+"Summary.txt", 'w')


'''
For the input of our simulations, we take in a file written in the standard 
formatting of duplex sequecing output
seqFile format:
    [0]refGenome [1]refNuc [2]refPos [3]seqDepth [4]totmuts [5]Tmuts [6]Cmuts
    [7]Gmuts [8]Amuts [9]Ins [10]Dels [11]Und
    
We make lists of every site that is sequenced as well as its associated sequencing
depth. The depth at each site is factored into our simulations, as we make each
site equally likely to "mutate," but weight the probability that we observe
a mutation there by its fractional representation amongst nucleotides of the same
type.  That is, if we want to redistribute a "C" mutation randomly, a C site that was
sequenced 1,000 times will be twice as likely to bear an observed mutation than
a C site that was sequenced 500 times.  That proportional representation is what
we are keeping track of by using these lists
'''

AList=[]
CList=[]
GList=[]
TList=[]

ASeqDepthList=[]
CSeqDepthList=[]
GSeqDepthList=[]
TSeqDepthList=[]

AFracList = []
CFracList = []
GFracList = []
TFracList = []


'''
obsTotMuts = 0
obsAMuts = 0
obsCMuts= 0
obsGMuts = 0
obsTmuts = 0
mutslist = []
mutSpectrumList = []
'''
def readSeqFile(inputFile):
    TotMuts = 0
    AMuts = 0
    CMuts= 0
    GMuts = 0
    TMuts = 0
    mutslist = []
    mutSpectrumList = []

    readfile = inputFile.readlines()
    for line in readfile[1:]:
        mutline = line.split()
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
    
        if int(mutline[4])> 0:
            TotMuts += int(mutline[4])
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
    return(TotMuts, AMuts, CMuts, GMuts, TMuts, mutslist, mutSpectrumList)
    
ObsTotMuts, ObsAMuts, ObsCMuts, ObsGMuts, ObsTMuts, Obsmutslist, ObsmutSpectrumList = readSeqFile(seqFile)
seqFile.close()
#pprint.pprint(Obsmutslist)

allMutTypes = ["AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"]
mutSpecDict = {}
for item in allMutTypes:
    mutSpecDict[item] = 0
for item in ObsmutSpectrumList:
    if str(item) not in mutSpecDict:
        next #CHANGE THIS BACK AT THE END, WAS A PRINT STATEMENT
    elif str(item) in mutSpecDict:
        mutSpecDict[item] += 1
#print mutSpecDict   
#for item in allMutTypes:
#    print str(item) + '\t' + str(mutSpecDict[item])


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
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        print s
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
        refbase=mut[0]
        mutbase=mut[1]
        mutpos = mut[2]
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
print len(AFracList)
print len(CFracList)
print len(GFracList)
print len(TFracList)

The following simulmut function is a somewhat inelegant way of redistributing mutations.
It does a random choice calculation for each muttype (eg, A>C) based on the proportion
each site was sequenced.  It is important that the list of sites is index matched
to the list of proportions. This will soon be replaced with vectorized calculations
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
    #pprint.pprint(simultotmuts)
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
obsPhyloP, obsCodingPhyloP, obsNCPhyloP, obsMutPred, obsNS, obsSyn, obsNSSRatio = parsemuts(Obsmutslist)

#THIS CAN ALL BE WRITTEN TO ONE FILE IF THIS IS PARALLELIZED SINCE IT IS
#ONLY BASED ON THE INPUT FILE
summaryOutputFile.write(str(sys.argv[1]) +'\n')
simulOutputFile.write(str(sys.argv[1]) +'\n')
summaryOutputFile.write("OBSERVED\n" +
    "observed PhyloP:\t" + str(obsPhyloP) + '\n' +
    "observed Coding PhyloP:\t" + str(obsCodingPhyloP) + '\n' +
    "observed NC PhyloP:\t" + str(obsNCPhyloP) + '\n' +
    "observed MutPred:\t" + str(obsMutPred) + '\n' +
    "observed #NS Muts:\t" + str(obsNS) + '\n' +
    "observed #S Muts:\t" + str(obsSyn) + '\n' +
    "observed NS:S Ratio:\t" + str(obsNSSRatio) + '\n'
    )
#summaryOutputFile.write(pprint.pformat(mutSpecDict))

 
totalsims=int(sys.argv[3])
for t in range(0,totalsims):
    if (t % 1000 == 0):
        print (t)
    simulPhyloP, simulCodingPhyloP, simulNCPhyloP, simulMutPred, simulNS, simulSyn, simulNSSRatio = parsemuts(simulmut())
#    print simulPhyloP
    simPhyloList.append(simulPhyloP)
    simCodPhyloPList.append(simulCodingPhyloP)
    simNCPhyloPList.append(simulNCPhyloP)
    simMutPredList.append(simulMutPred)
    simNSSRatioList.append(simulNSSRatio)

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
'''
THE BELOW PRINT STATEMENTS REFER TO THE AGGREGATED DATA, AKA HOW MANY SIMULTATIONS
OUT OF THE TOTAL ARE > THAN THE OBSERVED VALUE.  THIS NEEDS TO BE DONE ON ALL OF THE
SIMULATIONS, SO THAT IF PARALLELIZED, THIS WILL NEED TO CORRESPOND TO *ALL* OF
THE SIMULATIONS RUN ACROSS ALL NODES
'''
summaryOutputFile.write("SIMULATED\n" +
    "Simulated PhyloPs > Observed:\t" + str(PhyloAboveCount) + "/" + str(len(simPhyloList)) + '\n' +
    "Simulated Coding PhyloPs > Observed:\t" + str(CodPhyloAboveCount) + "/" + str(len(simPhyloList)) + '\n' +
    "Simulated NC PhyloPs > Observed:\t" + str(NCPhyloAboveCount) + "/" + str(len(simPhyloList)) + '\n' +
    "Simulated NSSRatio > Observed:\t" + str(NSSRatioAboveCount) + "/" + str(len(simPhyloList)) + '\n' +
    "Simulated MutPred > Observed:\t" + str(MutPredAboveCount) + "/" + str(len(simPhyloList)) + '\n' +
    pprint.pformat(mutSpecDict)
    )

simulOutputFile.write("Observed PhyloP \t Observed Coding PhyloP \t Observed NC PhyloP" +
    "\t Observed NS:S Ratio \t Observed MutPred \n" +
    str(obsPhyloP) + "\t" + str(obsCodingPhyloP) + "\t" + str(obsNCPhyloP) + "\t" +
    str(obsNSSRatio) + "\t" + str(obsMutPred) + '\n\n\n'+
    "Simulated PhyloP \t Simulated Coding PhyloP \t Simulated NC PhyloP" +
    "\t Simulated NS:S Ratio \t Simulated MutPred \n"
    )
   
for p in range(0, totalsims):
    
    simulOutputFile.write(
    str(simPhyloList[p]) + '\t' +
    str(simCodPhyloPList[p]) + '\t' +
    str(simNCPhyloPList[p]) + '\t' +
    str(simNSSRatioList[p]) + '\t' +
    str(simMutPredList[p]) + '\n'
    )    

simulOutputFile.close()
summaryOutputFile.close()