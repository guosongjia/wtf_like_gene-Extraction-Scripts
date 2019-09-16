#2019-02-21 10:47:11 February Thursday the 07 week, the 052 day SZ
#Script Owner: Guosong Jia
#Discription: This script is for extracting the up_and_down-stream sequence of wtf-like genes in the reference genome. Besides, this script also check
#             whether the extracted sequence is unique in S.octosporus genome.
#Usage: python wtf-like.locusUpAndDownSearch_19.02.21.py wtf-like_list[1] SOCG_whole_list[2] 


import os
import sys
import commands


with open(sys.argv[1]) as WtfLikeListFile:
	with open(sys.argv[2]) as SOCGWholeListFile:
		WtfLikeList = [n for n in WtfLikeListFile.readlines()]
		SOCGWholeList = [m for m in SOCGWholeListFile.readlines()]

def GetSequenceFile(FileName,Name,Sequence):
	'''
	This function is for sequence extraction. Using input information, the function can generate a file for "samtools faidx" sequence extraction.
	'''
	with open(FileName,'w') as Input:
		Input.write(">" + str(FileName)[:-12]+"-"+str(Name)+"\n"+str(Sequence))
def BLASTResultToBED(FileName):
	'''
	This function can use BLAST search result as input to generate a bed file for downstream locus extraction.
	'''
	OutName = str(FileName) + ".bed"
	ResultList = list()
	with open(FileName) as BLASTResult:
		for lines in BLASTResult.readlines():
			mod_lines = lines.split("\t")
			ResultList.append(str(mod_lines[1]) + "\t" + str(mod_lines[8]) + "\t" + str(mod_lines[9]) + "\t" + str(mod_lines[0]) + "\n")
	OutputFile=open(str(OutName),'w')
	for items in ResultList:
		OutputFile.write(items)
	OutputFile.close()

SuccessLog = open("Success.Up.log",'w')
FailLog = open("Fail.Up.log",'w')
LogFile = open("UpStreamLocus.log",'w')
for wtf_like_name in WtfLikeList:
	LogFile.write(str(wtf_like_name)[:-1] + "\n")
	print(str(wtf_like_name)[:-1])
	MatchHit = "".join([n for n in SOCGWholeList if str(wtf_like_name)[:-1] in n])
	MatchIndex = SOCGWholeList.index(MatchHit)
	i = 1
	while i < 7:
		MatchChromosome = str(SOCGWholeList[int(MatchIndex)].split("\t")[0])
		UpChromosome = str(SOCGWholeList[int(MatchIndex)-int(i)].split("\t")[0])
		if str(MatchChromosome) == str(UpChromosome):
			UpStreamName = str(SOCGWholeList[int(MatchIndex)-int(i)].split("\t")[3])
			UpStreamLocus = str(SOCGWholeList[int(MatchIndex)-int(i)].split("\t")[0]) + ":" + str(SOCGWholeList[int(MatchIndex)-int(i)].split("\t")[1]) + "-" + str(SOCGWholeList[int(MatchIndex)-int(i)].split("\t")[2])
			LogFile.write(str(UpStreamLocus) + "\n") 
			print(UpStreamLocus)
			UpstreamSequence = str(os.popen('samtools faidx ./TongModified_oct_genome.fa {0} |grep -v ">"'.format(UpStreamLocus)).read())[:-1]
			UpName = str(wtf_like_name)[:-1] + ".UP.inter.fasta"
			GetSequenceFile(UpName,UpStreamName,UpstreamSequence)
			os.system('~/Software/ncbi-blast-2.7.1/bin/blastn -query {0} -db ./genomeBLASTdatabase/TongModified.octosporus.database -outfmt 7 -out {0}.UP.out'.format(UpName))
			UpNums = os.popen('cat {0}.UP.out|grep -v "#"|grep "SOCG"|grep -v "^$" |wc -l '.format(UpName)).read()
			print(str(wtf_like_name)[:-1] +" is the " +str(i) + " upstream gene, which has " + str(UpNums)[:-1] + " BLAST results.")
			LogFile.write("Now processing " + str(wtf_like_name)[:-1] +", "+str(SOCGWholeList[int(MatchIndex)-int(i)].split("\t")[3])+" is the " + str(i) + " upstream gene, which has " + str(UpNums)[:-1] + " BLAST results.\n")
			if int(UpNums) == 1:
				print(str(UpStreamName) + " is the upstream unique gene of " + str(wtf_like_name)[:-1])
				LogFile.write(str(UpStreamName) + " is the upstream unique gene of " + str(wtf_like_name)[:-1] + ".\n")
				os.system('cat {0} >> Total.Up.fasta'.format(UpName))
				os.system('''less {0}.UP.out |grep -v "#"  >>Total.UP.out '''.format(UpName))
				os.system('rm {0}.UP.out '.format(UpName))
				os.system('rm {0}'.format(UpName))
				SuccessLog.write(str(wtf_like_name)[:-1] + "\n")
				break
			else:
				i = i + 1
				os.system('rm {0}'.format(UpName))
				os.system('rm {0}.UP.out'.format(UpName))
		elif str(MatchChromosome) != str(UpChromosome):
			print(str(wtf_like_name)[:-1] + " failed in extracting upstream unique gene. Reason: No other unique gene appears till the end of chromosome.")
			FailLog.write(str(wtf_like_name)[:-1] + " failed in extracting upstream unique gene. Reason: No other unique gene appears till the end of chromosome.\n")
			break
	else:
		FailLog.write(str(wtf_like_name)[:-1] + " failed in extracting upstream unique gene. Reason: No unique genes in 7 flanking-up genes.\n")
		continue
SuccessLog.close()
FailLog.close()
LogFile.close()

SuccessLog = open("Success.Down.log",'w')
FailLog = open("Fail.Down.log",'w')
LogFile = open("DownStreamLocus.log",'w')
for wtf_like_name in WtfLikeList:
	LogFile.write(str(wtf_like_name)[:-1] + "\n")
	print(str(wtf_like_name)[:-1])
	MatchHit = "".join([n for n in SOCGWholeList if str(wtf_like_name)[:-1] in n])
	MatchIndex = SOCGWholeList.index(MatchHit)
	i = 1
	while i < 7:
		MatchChromosome = str(SOCGWholeList[int(MatchIndex)].split("\t")[0])
		DownChromosome = str(SOCGWholeList[int(MatchIndex)+int(i)].split("\t")[0])
		if str(MatchChromosome) == str(DownChromosome):
			DownStreamName = str(SOCGWholeList[int(MatchIndex)+int(i)].split("\t")[3])
			DownStreamLocus = str(SOCGWholeList[int(MatchIndex)+int(i)].split("\t")[0]) + ":" + str(SOCGWholeList[int(MatchIndex)+int(i)].split("\t")[1]) + "-" + str(SOCGWholeList[int(MatchIndex)+int(i)].split("\t")[2])
			LogFile.write(str(DownStreamLocus) + "\n")
			print(DownStreamLocus)
			DownstreamSequence = os.popen('samtools faidx ./TongModified_oct_genome.fa {0} |grep -v ">"'.format(DownStreamLocus)).read()
			DownName = str(wtf_like_name)[:-1] + ".DOWN.inter.fasta"
			GetSequenceFile(DownName,DownStreamName,DownstreamSequence)
			os.system('~/Software/ncbi-blast-2.7.1/bin/blastn -query {0} -db ./genomeBLASTdatabase/TongModified.octosporus.database -outfmt 7 -out {0}.DOWN.out'.format(DownName))
			DownNums = os.popen('cat {0}.DOWN.out|grep -v "#"|grep -v "^$" |wc -l'.format(DownName)).read()
			print(str(wtf_like_name)[:-1]+" is the " + str(i) + " downstream gene, which has " + str(DownNums)[:-1] + " BLAST results.")
			LogFile.write("Now processing " + str(wtf_like_name)[:-1] +", "+str(SOCGWholeList[int(MatchIndex)+int(i)].split("\t")[3])+" is the " +str(i) + " downstream gene, which has " + str(DownNums)[:-1] + " BLAST results.\n")
			if int(DownNums) == 1:
				print(str(DownStreamName) + " is the downstream unique gene of " + str(wtf_like_name)[:-1])
				LogFile.write(str(DownStreamName) + " is the downstream unique gene of " + str(wtf_like_name)[:-1] + ".\n")
				os.system('cat {0} >> Total.Down.fasta'.format(DownName))
				os.system('''cat {0}.DOWN.out |grep -v "#"  >>Total.DOWN.out '''.format(DownName))
				os.system('rm {0}.DOWN.out'.format(DownName))
				os.system('rm {0}'.format(DownName))
				SuccessLog.write(str(wtf_like_name)[:-1] + "\n")
				break
			else:
				i = i + 1
				os.system('rm {0}'.format(DownName))
				os.system('rm {0}.DOWN.out'.format(DownName))
		elif str(MatchChromosome) != str(DownChromosome):
			print(str(wtf_like_name)[:-1] + " failed in extracting downstream unique gene. Reason: No other unique gene appears till the end of chromosome.")
			FailLog.write(str(wtf_like_name)[:-1] + " failed in extracting downstream unique gene. Reason: No other unique gene appears till the end of chromosome.\n")
			break
	else:
		FailLog.write(str(wtf_like_name)[:-1] + " failed in extracting downstream unique gene. Reason: No unique genes in 7 flanking-down genes.\n")
		continue
SuccessLog.close()
FailLog.close()
LogFile.close()

BLASTResultToBED("Total.UP.out")
BLASTResultToBED("Total.DOWN.out")