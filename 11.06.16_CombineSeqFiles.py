'''
this program simply takes the output files generated from duplex sequencing and
aggregates these files when we wish to combine data
In these files, each nucleotide position has a new row of data with the following
seqFile format:
    [0]refGenome [1]refNuc [2]refPos [3]seqDepth [4]totmuts [5]Tmuts [6]Cmuts
    [7]Gmuts [8]Amuts [9]Ins [10]Dels [11]Und
    
To aggregate files, this program checks checks whether each [2]refPos 
in a sequencing file has been seen in another sequencing data file that we are
combining with.  If it has not, it adds the data obtained from that animal for that
nucleotide.  If it has, it appends the new number of additional reads, and the
new number of each mutation type observed.  Then it generates an output
file of the same format now containing aggregated sequencing depth and mutations
at each site.
'''
seqDict = {}

def addfiletodict(seqFile):
    for line in seqFile:
        v = line.splitlines()
        for x in v:
            mutline = x.split() 
            if mutline[2] not in seqDict:
                seqDict[mutline[2]] = [mutline[0], mutline[1], mutline[2], 
                int(mutline[3]), int(mutline[4]), int(mutline[5]), 
                int(mutline[6]), int(mutline[7]), int(mutline[8]),
                int(mutline[9]), int(mutline[10]), int(mutline[11])]
            elif mutline[2] in seqDict:
                seqDict[mutline[2]][3] += int(mutline[3])
                seqDict[mutline[2]][4] += int(mutline[4])
                seqDict[mutline[2]][5] += int(mutline[5])
                seqDict[mutline[2]][6] += int(mutline[6])
                seqDict[mutline[2]][7] += int(mutline[7])
                seqDict[mutline[2]][8] += int(mutline[8])
                seqDict[mutline[2]][9] += int(mutline[9])
                seqDict[mutline[2]][10] += int(mutline[10])
                seqDict[mutline[2]][11] += int(mutline[11])
                
seqFile1 = open("file1.mutpos", "r")
seqFile2 = open("file2.mutpos", "r")
seqFile3 = open("file3.mutpos", "r")
seqFile4 = open("file4.mutpos", "r")
seqFile5 = open("file5.mutpos", "r")
combinedfile = open("CombinedFiles.txt", "w")
addfiletodict(seqFile1)
addfiletodict(seqFile2)
addfiletodict(seqFile3)
addfiletodict(seqFile4)
addfiletodict(seqFile5)


for x in range(0, 20000):
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
        str(seqDict[str(x)][8]) + '\t'+
        str(seqDict[str(x)][9]) + '\t'+
        str(seqDict[str(x)][10]) + '\t'+
        str(seqDict[str(x)][11]) + '\n'
        )


seqFile1.close()
seqFile2.close()
seqFile3.close()
seqFile4.close()
seqFile5.close()
combinedfile.close()
