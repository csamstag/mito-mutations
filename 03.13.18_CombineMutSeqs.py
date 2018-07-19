

seqDict = {}
'''

this program simply takes the output files generated from duplex sequencing and
aggregates these files when we wish to combine data
In these files, each nucleotide position has a new row of data with the following
seqFile format:
    [0]refGenome [1]refNuc [2]refPos [3]seqDepth [4]totmuts [5]Tmuts [6]Cmuts
    [7]Gmuts [8]Amuts [9]Ins [10]Dels [11]Und~/Documents/2018-06-18_allFliesNewestAlignment/AllMutPos/SNVs/70Clonal/CombinedNoDNoRNA/CombAllAgesByGenotype
    
when aggregating files, this program just checks checks whether each [2]refPos 
in a sequencing file has been seen in another sequencing data file that we are
combining with.  If it has not, it adds the data obtained from that animal for that
nucleotide.  If it has, it appends the new number of additional reads, and the
new number of each mutation type observed.  Then it generates an output
file of the same format now containing aggregated sequencing depth and mutations
at each site.

To alter, simply add each file name you want to aggregate within the "inputFileList", 
below, then change the output file name in "combinedfile"
'''
inputFileList =["wDah_1_SNV_70.txt",
                "wDah_2_SNV_70.txt",
                "wDah_3_SNV_70.txt"
                #"Fly_BSC+_50d-4.DCS.pileup.mutpos_SNV_0.7_RepeatMasked.txt",
                #"NM_20.DCS.pileup.mutpos_SNV_0.7_RepeatMasked.txt"                
                ]

combinedfile = open("wDah_Combined.txt", "w")


def addfiletodict(seqFile):
    seqFile.readline()
    for line in seqFile:
        v = line.splitlines()
        for x in v:
            mutline = x.split() 
            if mutline[2] not in seqDict:
                seqDict[mutline[2]] = [mutline[0], mutline[1], mutline[2], 
                int(mutline[3]), int(mutline[4]), int(mutline[5]), 
                int(mutline[6]), int(mutline[7]), int(mutline[8])]
            elif mutline[2] in seqDict:
                seqDict[mutline[2]][3] += int(mutline[3])
                seqDict[mutline[2]][4] += int(mutline[4])
                seqDict[mutline[2]][5] += int(mutline[5])
                seqDict[mutline[2]][6] += int(mutline[6])
                seqDict[mutline[2]][7] += int(mutline[7])
                seqDict[mutline[2]][8] += int(mutline[8])
                

 

for inputFile in inputFileList:
    seqFile = open(inputFile, "r")
    addfiletodict(seqFile)
    seqFile.close()


combinedfile.write("CHR\tREF\tPOS\tDCS\tMUTS\tT\tC\tG\tA\n")

for x in range(1, 19524):
    if str(x) in seqDict:
        combinedfile.write(
        seqDict[str(x)][0] + '\t'+
        seqDict[str(x)][1] + '\t'+
        str(seqDict[str(x)][2]) + '\t'+
        str(seqDict[str(x)][3]) + '\t'+
        str(seqDict[str(x)][4]) + '\t'+
        str(seqDict[str(x)][5]) + '\t'+
        str(seqDict[str(x)][6]) + '\t'+
        str(seqDict[str(x)][7]) + '\t'+
        str(seqDict[str(x)][8]) + '\n'
        )


combinedfile.close()