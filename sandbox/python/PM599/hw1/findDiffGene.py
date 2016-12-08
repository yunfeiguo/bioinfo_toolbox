#!/usr/bin/env python

infile = "Piwi_vs_WT.averaged.RPKM.txt"
outfile1 = 'upGeneList.txt'
outfile2 = 'strongUpGeneList.txt'
outfile3 = 'downGeneList.txt'
outfile4 = 'strongDownGeneList.txt'
in_fh = open(infile,'r')
out_fh1 = open(outfile1,'w')
out_fh2 = open(outfile2,'w')
out_fh3 = open(outfile3,'w')
out_fh4 = open(outfile4,'w')

for i in in_fh:
    l = i.split()
    gene = l[0]
    rpkm1 = float(l[1])
    rpkm2 = float(l[2])
    if rpkm2 == 0:
	ratio = float("inf")
    else:
	ratio = rpkm1/rpkm2
    if max(rpkm1,rpkm2) >= 1 and ratio >= 2:
	#upregulated
	out_fh1.write(gene+'\n')
    if max(rpkm1,rpkm2) >= 1 and ratio >= 3:
	#strongly upregulated
	out_fh2.write(gene+'\n')
    if max(rpkm1,rpkm2) >= 1 and ratio <= 1/2:
	#downregulated
	out_fh3.write(gene+'\n')
    if max(rpkm1,rpkm2) >= 1 and ratio <= 1/3:
	#strongly downregulated
	out_fh4.write(gene+'\n')
in_fh.close()
out_fh1.close()
out_fh2.close()
out_fh3.close()
out_fh4.close()
print "All done"
