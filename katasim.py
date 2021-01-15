#!/usr/bin/python3

import sys
import random

helpdoc="Created by Xue Lin <xue.lin@njmu.edu.cn>\n\n\
Description：	This is a tool to mutate any DNA sequence to generate sequences with or without kataegis.\n\n\
Usage:	python katasim.py -i <input file> -o <output file> -v <output vcf> -k 1\n\n\
Options:\n\
-i The input file of the fasta format.\n\
-o	The output file of the fasta format.\n\
-v The output of added mutations of VCF format.\n\n\
-k	A boolean value whether to add kataegis, 1 is yes and 0 is no, default is 0.\n\
-kn	An integer to define how many kataegises to add to the sequence. If the customer input kn is larger than the biggest value according to the length of the seqeunce, the program will give the biggist value and return a warning.\n\
-tmb	A float number to define the tumor mutation burden of the total sequence, default is 100 mutations/10^6bp, i.e. 0.0001.\n"

katalen=1000 #the length of the range to define kataegis, usually is 10kb.
katamut=6 #the number of the mutations to define kataegis, usually is 5.

def getopts():
	k=0
	tmb=0.0001
	kn=0
	infile=''
	outfile=''
	vcfout=''
	if len(sys.argv[1:])>0 and len(sys.argv[1:])%2==0:
		#print(sys.argv[1:])
		opts=sys.argv[1:]
		key=opts[::2]
		val=opts[1::2]
		if '-i' in key and '-o' in key and '-v' in key:
			args=dict(zip(key,val))
			#print(args)
			for key in args:
				if key=='-i':
					try:
						infile=open(args[key],'r')
					except:
						print(helpdoc)
						print("WARNINGS: The file you input cannot be open correctly, please check the file name or directory!\n")
						exit(1)
					#infile=args[key]s
				if key=='-o':
					try:
						outfile=open(args[key],'w')
					except:
						print(helpdoc)
						print("WARNINGS: The file you output cannot be created, please check the -o parameter settings!\n")
						exit(1)
					#outfile=args[key]
				if key=='-v':
					try:
						vcfout=open(args[key],'w')
					except:
						print(helpdoc)
						print("WARNINGS: The file you output cannot be created, please check the -v parameter settings!\n")
						exit(1)
					#vcfout=args[key]
				if key=='-k':
					if args[key]=="1" or args[key]=="0":
						k=int(args[key])
					else:
						print(helpdoc)
						print("WARNINGS: The -k should be 1 or 0, please check your -k settings!\n")
						exit(1)
				if key=='-tmb':
					try:
						tmb=float(tmb)
					except(ValueError):
						print(helpdoc)
						print("WARNINGS: The -tmb should be a float, please check your -tmb settings!\n")
						exit(1)
					if tmb>1:
						print(helpdoc)
						print("WARNINGS: The -tmb should be between 0 and 1, please check your -tmb settings!\n")
				if key=='-kn':
					try:
						kn=int(args[key])
					except(ValueError):
						print(helpdoc)
						print("WARNINGS: The -kn should be integer, please check your -kn settings!\n")
						exit(1)
			return infile,outfile,vcfout,k,kn,tmb
		else:
			print(helpdoc)
			print("ERRORS: The input and output files are obliged, please check the -i, -o, -v settings!\n")
			exit(1)	
	else:
		print(helpdoc)
		exit(1)

def mut(infile,outfile,vcfout,k,kn,tmb):
	#check kn
	size=0
	value=''
	genome={}
	for line in infile.readlines():
		line=line.strip("\n")
		if line[0]==">":
			key=line[1:]
		else:
			value+=line.upper()
	size=len(value)
	value=list(value)
	maxseg=int(size/katalen)
	maxmut=int(size*tmb)
	print(str(size)+"\n"+str(maxmut)+"\n"+str(maxseg))
	if k==1:
		if maxmut>=kn*katamut:
			#parameters are correctly setted.
			if kn==0:
				#give a default kn setting or the user give a wrong parameter,
				#at this senario, the kn will be calculated as the largest according to he seuqncine size.
				kn=maxseg
			elif kn>0 and kn<=maxseg:
				#the kn is correctly setted.
				kn=kn
			else:
				print("ERROR！The kn is too large according to the size of the genome!\n")
				print("The program will set the kn to "+maxseg+" automaticallys.\n")
				kn=maxseg
		else:
			print("ERROR! The tumor mutation burden is not sufficient for calling the number of kataegis you set!\n")
			print("Please set larger tmb and smaller kn.\n")
			exit(1)
		#start the mutation process
		#to mutate the sequence with kataegis
		#the mutation is happened in randomly selected regions of the kataegis
		vcfout.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
		loc=[]
		div=[]
		tmp=[]
		kataregion=[]
		for i in range(0,maxseg):
			div.append(i*katalen) #left end
		kataregion=random.sample(div,kn)
		kataregion.sort()
		i=0
		for sep in kataregion:
			if i==0:
				for i in range(0,katamut):
					tmp.append(random.randrange(sep,sep+katalen))
				tmp.sort()
				loc+=tmp
			else:
				while tmp[0]-loc[-1]<katalen:
					tmp=[]
					for i in range(0,katamut):
						tmp.append(random.randrange(sep,sep+katalen))
					tmp.sort
				loc+=tmp
			i+=1
		#print(len(loc))
		div=[]
		i=0
		for sep in kataregion:
			if sep==0:
				div.append(seq+katalen+1) #left end
			else:
				if i==0:
					div.append(1)
					div.append(sep)
					div.append(sep+katalen+1)
				else:
					div.append(sep)
					div.append(sep+katalen+1)
				
			i+=1
		if len(div)%2==1:
			div.append(size)
		#print(len(div)%2)
		thres=len(loc)
		while len(loc)<maxmut:
			for i in range(0,len(div),2):
				times=random.randrange(1,katamut)
				if i==0:
					tmp=[]
					for j in range(0,times):
						tmp.append(random.randrange(div[i],div[i+1]))
					tmp.sort()
				else:
					while div[i+1]-tmp[-1]<katalen or tmp[0]-div[i]<katalen:
						tmp=[]
						for j in range(0,times):
							tmp.append(random.randrange(div[i],div[i+1]))
						tmp.sort()
					loc+=tmp
		print(kataregion)
		sub=['A','T','C','G']
		x=0
		for i in loc:
			x+=1
			origin=value[i]
			value[i]=sub[random.randrange(0,len(sub))]#the substitutional nucleotide is generated randomly
			mutated=value[i]
			if x <=thres:
				vcfout.write(str(key)+"\t"+str(i)+"\t"+'.'+"\t"+origin+"\t"+mutated+"\t"+'.'+"\t"+"PASS"+"\t"+"KATA\n")
			else:
				vcfout.write(str(key)+"\t"+str(i)+"\t"+'.'+"\t"+origin+"\t"+mutated+"\t"+'.'+"\t"+"PASS"+"\t"+".\n")
	else:
		#to mutate the sequnce without kataegis
		#the mutation is randomly happened on genome.
		vcfout.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n")
		loc=[]
		while len(loc)<maxmut:
			loc.append(random.randrange(1,size+1))
			if len(loc)>=katamut:
				loc.sort()
				j=0
				for pos in loc:
					if (j+katamut)<=len(loc):
						if (loc[j+katamut-1]-pos+1)<=katalen:
							abe=loc.pop(j+katamut-1)
							#print("The element "+str(abe)+" is removed, and the indice is "+str(j+katamut-1)+".\n")
					j+=1
		sub=['A','T','C','G']
		for i in loc:
			origin=value[i]
			value[i]=sub[random.randrange(0,len(sub))]#the substitutional nucleotide is generated randomly
			mutated=value[i]
			vcfout.write(str(key)+"\t"+str(i)+"\t"+'.'+"\t"+origin+"\t"+mutated+"\t"+'.'+"\t"+"PASS"+"\n")
	
	#return loc
			
	
pars=getopts()#infile,outfile,vcfout,k,kn,tmb

mut(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5])
#print(loc)
