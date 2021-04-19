##############################################################################################################################
#																															 #
#      sci-L3-target-seq Primer Designer																					 #
#																															 #
#																															 #
# Author: Nidhi Koundinya																									 #
# Date Created: March 31, 2021																								 #
##############################################################################################################################
# Written for Python 3.6
#
# Usage:
# python3 run.py -i <path to input BED or fasta file> -g <path to reference genome> [options]
#
# Arguments:
# 
#		-h, --help                  show this help message and exit
#		-i, --input                 Enter the path to the BED or FASTA file
#		-g, --genome                Enter the path to the reference genome
#		--P5, -P5                   P5 Adapter sequence(Default=AATGATACGGCGACCACCGAGA)
#		--P7, -P7                   P7 Adapter sequence(Default=CAAGCAGAAGACGGCATACGAGAT)
#		--read1, -r1                Read1 sequence(Default=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG)
#		--read2, -r2                Read2 sequence(Default=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG)
#		--SSSPR, -ssspr             SSS Primer Region(Default=GGGATGCAGCTCGCTCCTG)
#		--indices, -indices         Enter the path to the CSV file containing sample specific indices
#		--barcodes, -b              Enter the path to the CSV file containing Barcodes
#
#
#
# Inputs:
# 1. A 6 column BED file containing coordinates of regions to be targeted
# 2. A CSV file with a list of sample specific indices.
# 3. A CSV file with a list of barcodes.
#
# Outputs:
# 1. A directory "primer_out" at current location with directory "final" inside.
# 2. A csv file "primers.csv" inside the final directory with a list of primers.
# 3. A csv file "forwardPrimers.csv" inside the final directory with a list of the forward primers.
#
#
############################################################################################################################


import os
import re
import math
import time
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO


new_path = os.path.join(os.getcwd() + "/primer_out")
temp_path = os.path.join(new_path, "temp")
final_path = os.path.join(new_path, "final")
seqs_path = os.path.join(temp_path, "sequences.fasta")

#Deletes any folder with primers named "primer_out"
if os.path.isdir(new_path):
	rm_cmd = ["rm","-rf",new_path]
	subprocess.run(rm_cmd)

#Deletes temp file from "primer_out"
def delete_temp():
	if os.path.isdir(temp_path):
		rm_cmd = ["rm","-rf",temp_path]
		subprocess.run(rm_cmd)

os.mkdir(new_path)
os.mkdir(temp_path)
os.mkdir(final_path)
fragment_size_max=200


#parsing command line arguments
def command_Parse():
	parser = argparse.ArgumentParser(prog='python3 run.py', usage='%(prog)s -i <path to input BED or fasta file> -g <path to reference genome> [options]',
    		description='Generates primers for sci-L3-target-seq', epilog="***** sci-L3-target-seq - PRIMER DESIGNER *****\n")
	def check_input(file):
		base, ext = os.path.splitext(file)
		if not os.path.exists(file):
			raise argparse.ArgumentTypeError('{} not found. Check path.'.format(file))
		elif ext.lower() not in ('.bed', '.bed'):
			raise argparse.ArgumentTypeError('\nInput must be a BED file or a FASTA file')
		return file
	
	def check_ref(file):
		base, ext = os.path.splitext(file)
		if not os.path.exists(file):
			raise argparse.ArgumentTypeError('{} not found. Check path.'.format(file))
		elif ext.lower() not in ('.fasta'):
			raise argparse.ArgumentTypeError('\nInput must be a FASTA file.')
		return file
	
	def check_indices(file):
		base, ext = os.path.splitext(file)
		if not os.path.exists(file):
			raise argparse.ArgumentTypeError('{} not found. Check path.'.format(file))
		elif ext.lower() not in ('.csv', '.csv'):
			raise argparse.ArgumentTypeError('\nIndices must be in a CSV file.')
		return file
	
	def check_barcodes(file):
		base, ext = os.path.splitext(file)
		if not os.path.exists(file):
			raise argparse.ArgumentTypeError('{} not found. Check path.'.format(file))
		elif ext.lower() not in ('.csv', '.csv'):
			raise argparse.ArgumentTypeError('\nBarcodes must be in a CSV file.')
		return file
	
	##Arguments
	parser.add_argument("-i", "--input", help='Enter the path to the BED or FASTA file', type=check_input , required=True)
	parser.add_argument("-g", "--genome", help='Enter the path to the reference genome', type=check_ref , required=True)
	
	parser.add_argument("--P5", "-P5", type=str, default="AATGATACGGCGACCACCGAGA", help="P5 Adapter sequence(Default=AATGATACGGCGACCACCGAGA)")
	parser.add_argument("--P7", "-P7", type=str, default="CAAGCAGAAGACGGCATACGAGAT", help="P7 Adapter sequence(Default=CAAGCAGAAGACGGCATACGAGAT)")
	
	parser.add_argument("--read1", "-r1", type=str, default="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", help="Read1 sequence(Default=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG)")
	parser.add_argument("--read2", "-r2", type=str, default="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG", help="Read2 sequence(Default=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG)")
	parser.add_argument("--SSSPR", "-ssspr", type=str, default="GGGATGCAGCTCGCTCCTG", help="SSS Primer Region(Default=GGGATGCAGCTCGCTCCTG)")
	
	parser.add_argument("--indices", "-indices", help='Enter the path to the CSV file containing sample specific indices', type=check_indices, required=True)
	parser.add_argument("--barcodes", "-b", help='Enter the path to the CSV file containing Barcodes', type=check_barcodes, required=True)
		
	return parser


class DesignPrimer:
	def __init__(self, args):
		self.input=args.input
		self.genome=args.genome
		self.P5=args.P5 #fwd
		self.P7=args.P7	#rev
		self.read1=args.read1	#fwd
		self.read2=args.read2	#rev
		self.SSSPR=args.SSSPR	#fwd
		self.indices=args.indices	#rev
		self.barcodes=args.barcodes	#fwd
		
		self.bed_input=False #Default
	
	
	##Get sequences from reference human genome
	def getSequences(self):
		base, ext = os.path.splitext(self.input)
		self.ext=ext
		if ext=='.bed':
			self.bed_input=True
			bed_cmd = ["bedtools", "getfasta", "-fi", self.genome, "-bed", self.input, "-fo", seqs_path, "-name", "-s"]
			subprocess.run(bed_cmd)
			for record in SeqIO.parse(seqs_path,"fasta"):
				record_id = record.id.split('::')[0]
				cmd="sed -i 's/>"+record.id+"/>"+record_id+"/' " + seqs_path
				os.system(cmd)
		elif ext=='.fasta':
			fasta_cmd=["cp", self.input, seqs_path]
			subprocess.run(fasta_cmd)
		return self.bed_input
	
	
	## Fragment each target region into fragments less than 200bp
	def fragment_sequences(self,record_sequence):
		fragment_size=201
		fragments=[]
		i=1
		pos=0
		seq_length = len(record_sequence)
		while (fragment_size > fragment_size_max):
			fragment_size = math.ceil(seq_length/i)
			i+=1
			if (fragment_size < fragment_size_max):
				break
		while pos < len(record_sequence):
			fragments.append(record_sequence[pos:pos+fragment_size])
			pos += fragment_size
		return fragments
	
	
	##Creates input file for primer3
	def make_p3Input(self,header,seq,primer_values):
		with open(f'{temp_path}/inputs.txt','a') as txt_file:
			txt_file.write('SEQUENCE_ID={}\nSEQUENCE_TEMPLATE={}\n'.format(header.upper(),seq.upper()))

			for i,items in enumerate(primer_values):
				txt_file.write(items)
				if i==len(primer_values)-1:
					break
				txt_file.write('\n')
			txt_file.write('\n')
	
	
	##Makes target length and product size settings for each region
	def make_p3Settings(self):
		for record in SeqIO.parse(seqs_path,"fasta"):
			record_id = record.id
			record_sequence = str(record.seq)
			fragments = DesignPrimer.fragment_sequences(self,record_sequence)
	
			for index,seq in enumerate(fragments):
				header = '{}_{}'.format(record_id,index+1)
				primer_values = ["PRIMER_TASK=generic","PRIMER_PICK_LEFT_PRIMER=1","PRIMER_PICK_INTERNAL_OLIGO=0","PRIMER_PICK_RIGHT_PRIMER=1",
									"PRIMER_NUM_RETURN=1","PRIMER_OPT_SIZE=20","PRIMER_MIN_SIZE=16","PRIMER_MAX_SIZE=27","PRIMER_EXPLAIN_FLAG=1"]
				if (len(seq)<27):
					print(f'{header} is too short.\n')
				#for longer fragments
				if (len(seq)>61):
					targ=30
					target=f'SEQUENCE_TARGET=33,{targ}'            
					length=len(seq)
					prod_max_size=210
					a=int(.5*(length))
					prod_min_size=30
					product_range=f'PRIMER_PRODUCT_SIZE_RANGE={a}-{prod_max_size} {prod_min_size}-{prod_max_size}'
					primer_values.append(product_range)
					primer_values.append(target)
					primer_values.append("=")
				#for short fragments
				else:
					target=f'SEQUENCE_TARGET=33,5'
					product_range=f'PRIMER_PRODUCT_SIZE_RANGE=27-{prod_max_size}'
					primer_values.append(product_range)
					primer_values.append(target)
					primer_values.append("=")
				self.make_p3Input(header,seq,primer_values)


##Parsing primer3 output -> csv
def parseP3_Output(args):
	def generate_p3log(line,counter):
		with open(f'{final_path}/primer3_log.txt','a') as log_file:
			if re.search('(?<=PRIMER_RIGHT_EXPLAIN=).*$',line):
				right=re.findall('(?<=PRIMER_RIGHT_EXPLAIN=).*$',line)[0]
				log_file.write(f'Right Primers: {right}\n')
			if re.search(f'(?<=PRIMER_ERROR=).*$',line):
				log_file.write(f'{seq_ids[counter]}\n')
				err=re.findall('(?<=PRIMER_ERROR=).*$',line)[0]
				log_file.write(f'Error: {err}\n\n')
				counter+=1
	def parseIdAndTemplate():
		with open(f'{temp_path}/primerOutput','r') as p_out:
			for line in p_out:
				if re.search('SEQUENCE_ID=',line):
					m=re.findall('(?<=SEQUENCE_ID=).*$',line)[0]
					
					indices_df = pd.read_csv(f'{args.indices}', header=None)
					indices = indices_df.loc[:,0]
					for i,x in enumerate(indices):
						index=str(i+1)
						name=m+"_index"+index
						seq_ids.append(name)
				if re.search('SEQUENCE_TEMPLATE=',line):
					m=re.findall('(?<=SEQUENCE_TEMPLATE=).*$',line)
					templates.append(m[0])
	def parsePrimers():
		with open(f'{temp_path}/primerOutput','r') as p_out:
			for line in p_out:
				if re.search('SEQUENCE_ID=',line):
					next
				generate_p3log(line,counter)
				if re.search('(?<=PRIMER_RIGHT_0_SEQUENCE=).*$',line):
					seq=re.findall('(?<=PRIMER_RIGHT_0_SEQUENCE=).*$',line)[0]
					
					indices_df = pd.read_csv(f'{args.indices}', header=None)
					indices = indices_df.loc[:,0]
					for x in indices:
						primer=args.P7+x+args.read2+seq
						right_primers.append(primer)

	def createOutputFile():
		primers_out = {}
		primers_out['Sequence_ID'] = seq_ids
		primers_out['Right_Primers'] = right_primers
		data = list(zip(primers_out['Sequence_ID'],primers_out['Right_Primers']))
		df = pd.DataFrame(data=data)
		df.to_csv(f'{final_path}/primers.csv', index=False, header=fieldnames)
		print(f'primers.csv created in FINAL directory')
	seq_ids=[]
	templates=[]
	right_primers=[]
	fieldnames=["Sequence_ID","Primers"]
	parseIdAndTemplate()
	with open(f'{final_path}/primer3_log.txt','a') as log_file:
		log_file.write(f'\n')
		counter=0
	
	parsePrimers()
	createOutputFile()

##Running primer3
def runPrimer3():
	input_file = temp_path + "/inputs" + ".txt"
	
	#removing extra newline at the end of input file
	truncInput_cmd = ["truncate", "-s", "-1", input_file]
	subprocess.run(truncInput_cmd)
	
	#Running primer3
	output_file = "--output=" + temp_path + "/primerOutput"
	primer3_cmd = ["primer3_core", input_file , output_file]
	subprocess.run(primer3_cmd)


##Creating Forward primers based on barcodes given
def makeForwardPrimer(args):
	barcodes_df = pd.read_csv(f'{args.barcodes}', header=None)
	barcodes = barcodes_df.loc[:,0]
	seq_ids = []
	primers = []
	for i,x in enumerate(barcodes):
		index=str(i+1)
		seq_id="universalPrimer_"+index
		
		seq=args.P5+args.read1+x+args.SSSPR
		primers.append(seq)
		seq_ids.append(seq_id)
	primers_out = {}
	primers_out['Sequence_ID'] = seq_ids
	primers_out['Primers'] = primers
	data = list(zip(primers_out['Sequence_ID'],primers_out['Primers']))
	df = pd.DataFrame(data=data)
	
	fieldnames=["Sequence_ID","Primers"]
	df.to_csv(f'{final_path}/forwardPrimers.csv', index=False, header=fieldnames)
	print(f'forward primers.csv created in FINAL directory')
	
	

def main():
	parser = command_Parse()
	args = parser.parse_args()
	print()
	print("Input:\t",args.input)
	print("Genome:\t",args.genome,"\n")
	print("P5 Adapter sequence:\t",args.P5)
	print("P7 Adapter sequence:\t",args.P7,"\n")
	print("Read1 sequence:\t",args.read1)
	print("Read2 sequence:\t",args.read2,"\n")
	print("SSS Primer Region sequence:\t",args.SSSPR,"\n")
	print("Indices:\t",args.indices)
	print("Barcodes:\t",args.barcodes,"\n")
	
	dp = DesignPrimer(args)
	
	bed_input=dp.getSequences()
	print("Creating Input File for Primer3")
	dp.make_p3Settings()
	print("Input File Created")
	
	runPrimer3()
	parseP3_Output(args)
	makeForwardPrimer(args)
	
	delete_temp()


if __name__ == "__main__":
	main()
