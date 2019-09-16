#2019-08-27 12:25:12 August Tuesday the 34 week, the 239 day SZ
#Script Owner: Guosong Jia
#Description: This script is the downstream step of "wtf-like.locusUpAndDownSearch_19.02.21.py".
#             Based on the output files: "Total.UP.out.bed" and "Total.DOWN.out.bed"
#Usage: python wtf-like.locusExtraction_19.08.27.py wtf-like.list Total.UP.out.bed Total.DOWN.out.bed

import os
import sys
import commands
#Import "wtf-like.list", "Upstream-unique-gene.bed", "Downstream-unique-gene.bed" 
WtfLikeListFile = open(sys.argv[1])
UpStreamOut = open(sys.argv[2])
DownStreamOut = open(sys.argv[3])
WtfLikeList = [n[:-1] for n in WtfLikeListFile.readlines()]
UpstreamList = [n[:-1] for n in UpStreamOut.readlines()]
DownstreamList = [n[:-1] for n in DownStreamOut.readlines()]
WtfLikeListFile.close()
UpStreamOut.close()
DownStreamOut.close()
#Import "GetSequenceFile" function
def GetSequenceFile(FileName,Sequence):
	'''
	This function is for sequence extraction. Using input information, the function can generate a file for "samtools faidx" sequence extraction.
	'''
	with open(FileName,'w') as Input:
		Input.write(str(os.popen('echo -e ">{0}.locus"'.format(str(FileName))).read()) + str(Sequence))
#Compare Upstream-unique gene and Downstream-unique gene list and generate a list contains passed wtf-like genes.
UpMatchGeneList = list()
UpMatchList = list()
DownMatchGeneList = list()
DownMatchList = list()
for wtf_like_genes in WtfLikeList:
	for items in UpstreamList:
		if str(wtf_like_genes) in items:
			UpMatchGeneList.append(wtf_like_genes)
			UpMatchList.append(items)
	for items in DownstreamList:
		if str(wtf_like_genes) in items:
			DownMatchGeneList.append(wtf_like_genes)
			DownMatchList.append(items)
PassWtf_likeList = [n for n in UpMatchGeneList if n in DownMatchGeneList]
LogFile = open("ExtractionResult.log",'w')
for wtf_like_genes in WtfLikeList:
	if str(wtf_like_genes) in PassWtf_likeList:
		print(str(wtf_like_genes) + " has matched upstream-downstream unique gene records.")
		LogFile.write(str(wtf_like_genes) + " has matched upstream-downstream unique gene records.\n")
		UpMatchIndex = UpMatchList.index("".join([n for n in UpMatchList if str(wtf_like_genes) in n]))
		DownMatchIndex = DownMatchList.index("".join([n for n in DownMatchList if str(wtf_like_genes) in n]))
		chromosome = str(UpMatchList[UpMatchIndex].split("\t")[0])
		StartPos = str(UpMatchList[UpMatchIndex].split("\t")[1])
		EndPos = str(DownMatchList[DownMatchIndex].split("\t")[2])
		Locus = str(chromosome) + ":" + str(StartPos) + "-" + str(EndPos)
		LocusSequence = os.popen('samtools faidx ./TongModified_oct_genome.fa {0} |grep -v ">"'.format(Locus)).read()
		GetSequenceFile(str(wtf_like_genes)+".Locus.fasta",LocusSequence)
	elif str(wtf_like_genes) not in PassWtf_likeList:
		if str(wtf_like_genes) in UpMatchGeneList and str(wtf_like_genes) not in DownMatchGeneList:
			print(str(wtf_like_genes) + " failed. Reason: Lack Down-stream unique gene.")
			LogFile.write(str(wtf_like_genes) + " failed. Reason: Lack Down-stream unique gene.\n")
		elif str(wtf_like_genes) in DownMatchGeneList and str(wtf_like_genes) not in UpMatchGeneList:
			print(str(wtf_like_genes) + " failed. Reason: Lack Up-stream unique gene.")
			LogFile.write(str(wtf_like_genes) + " failed. Reason: Lack Up-stream unique gene.\n")
		elif str(wtf_like_genes) not in UpMatchGeneList and str(wtf_like_genes) not in DownMatchGeneList:
			print(str(wtf_like_genes) + " failed. Reason: Lack both Up-stream and Down-stream unique gene.")
			LogFile.write(str(wtf_like_genes) + " failed. Reason: Lack both Up-stream and Down-stream unique gene.\n")


#print(WtfLikeList)
#print(UpstreamList[1].split("\t")[0][:10])
#print(len(PassWtf_likeList))