#extracting just protein sequences
'''
240	1263	protein	ND2	1	0	0	0
1474	3009	protein	COX1	1	0	0	0
3083	3767	protein	COX2	1	0	0	0
3907	4068	protein	ATP8	1	0	0	0
4062	4736	protein	ATP6	1	0	0	0
4736	5524	protein	COX3	1	0	0	0
5608	5961	protein	ND3	1	0	0	0
6409	8125	protein	ND5	0	0	0	0
8207	9545	protein	ND4	0	0	0	0
9545	9835	protein	ND4L	0	0	0	0
9971	10495	protein	ND6	1	0	0	0
10499	11635	protein	CYTB	1	0	0	0
11721	12659	protein	ND1	0	0	0	0
'''
geneboundslist=[[240,1263],[1474,3009],[3083,3767],[3907,5524],[5608,5961],
    [6409,8125],[8207,9835],[9971,10495],[10499,11635],[11721,12659]
    ]    
    
seqFile =open("NM_20.DCS.pileup.mutpos_SNV_0.7.txt",'r')
outputFile = open("NM_20_SNV_0.7_proteinOnly.txt", "w")


seqFile.readline()
for line in seqFile:
    v = line.splitlines()
    for x in v:
        mutline = x.split() 
        refPos = int(mutline[2])
        for bound in geneboundslist:
            if refPos >= bound[0] and refPos <= bound[1]:
                outputFile.write(line)
                


seqFile.close()
outputFile.close()