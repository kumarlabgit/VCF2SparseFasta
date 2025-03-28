import argparse
import os
import subprocess
import math

bcft_path = "~/VCFconverter/bcftools/bcftools"
hg38_chrom_sizes = {"chr1": 249000000, "chr2": 243000000, "chr3": 199000000, "chr4": 191000000, "chr5": 182000000, "chr6": 171000000, "chr7": 160000000,
					"chr8": 146000000, "chr9": 139000000, "chr10": 134000000, "chr11": 136000000, "chr12": 134000000, "chr13": 115000000, "chr14": 108000000,
					"chr15": 102000000, "chr16": 91000000, "chr17": 84000000, "chr18": 81000000, "chr19": 59000000, "chr20": 65000000, "chr21": 47000000,
					"chr22": 51000000, "chrX": 157000000, "chrY": 58000000,
					"1": 249000000, "2": 243000000, "3": 199000000, "4": 191000000, "5": 182000000, "6": 171000000, "7": 160000000,
					"8": 146000000, "9": 139000000, "10": 134000000, "11": 136000000, "12": 134000000, "13": 115000000, "14": 108000000,
					"15": 102000000, "16": 91000000, "17": 84000000, "18": 81000000, "19": 59000000, "20": 65000000, "21": 47000000,
					"22": 51000000, "X": 157000000, "Y": 58000000}


def main(args):
	if os.path.splitext(args.vcf_path)[1] == ".vcf":
		if os.path.exists(args.vcf_path + ".gz"):
			print("Compressed version already exists, switching to {}.".format(args.vcf_path + ".gz"))
			args.vcf_path = args.vcf_path + ".gz"
		else:
			print("Compressing to .gz format for indexing")
			bcft_cmd = "{}*view*{}*-Oz*-o*{}".format(bcft_path, args.vcf_path, args.vcf_path + ".gz")
			subprocess.run(bcft_cmd.replace("*", " "), shell=True)
			args.vcf_path = args.vcf_path + ".gz"
	if os.path.splitext(args.vcf_path)[1] == ".gz" and os.path.splitext(os.path.splitext(args.vcf_path)[0])[1] == ".vcf":
		if os.path.exists(args.vcf_path + ".tbi"):
			print("Index file {} detected.".format(args.vcf_path + ".tbi"))
		else:
			print("Generating index file.")
			bcft_cmd = "{}*index*-t*{}".format(bcft_path, args.vcf_path, args.vcf_path)
			subprocess.run(bcft_cmd.replace("*", " "), shell=True)
	else:
		raise Exception("Please specify a *.vcf or *.vcf.gz file.")
	if args.output is None:
		args.output = os.path.splitext(os.path.basename(args.vcf_path))[0].replace(".vcf", "")
	if not os.path.exists(args.output):
		os.mkdir(args.output)
	if args.end is None:
		args.end = hg38_chrom_sizes[args.chrom]
	[region_name, region_start, region_end] = [args.chrom, args.start, args.end]
	temp_base = os.path.splitext(os.path.basename(args.vcf_path))[0]
	last_pos = max(region_start, 0)
#	temp_path = "tempfile.txt"
	bad_pos_list = []
	gt_dict = {}
	header = None
	while last_pos < region_end:
		temp_path = temp_base + "_temp.txt"
		region = "{}:{}-{}".format(region_name, last_pos, last_pos + args.query_size)
#		last_pos = args.query_size
		bcft_cmd = "{}*query*-i*'strlen(REF)=1 && strlen(ALT)=1'*--print-header*-f*'%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'*{}".format(bcft_path, args.vcf_path)
		bcft_cmd = "{}*-r*{}".format(bcft_cmd, region)
		if header is None and args.verbosity > 0:
			print(bcft_cmd.replace("*", " "))
		with open(temp_path, "w") as file:
			if args.verbosity > 0:
				print("Querying region {}.".format(region))
			subprocess.run(bcft_cmd.replace("*", " "), stderr=subprocess.STDOUT, stdout=file, shell=True)
		with open(temp_path, 'r') as file:
			if header is None:
				header = [x[x.find("]")+1:] for x in file.readline().strip().split('\t')]
				sample_ids = header[4:]
			else:
				file.readline()
			for line in file.readlines():
				data = line.strip().split('\t')
				if int(data[1]) in gt_dict.keys():
					[ref, alt] = data[2:4]
					if gt_dict[int(data[1])][0] != ref:
						if args.verbosity > 0:
							print("Returned multiple reference alleles for position {}, dropping position.".format(data[1]))
						bad_pos_list += [int(data[1])]
					else:
						gt = [ref] + [ref if x[0] == "0|0" and x[1] == ref else alt if x[1] == ref else x[1] if x[0] == "0|0" else "-" for x in zip(data[4:], gt_dict[int(data[1])][1:])]
						if "-" in gt:
							#bad_pos_list += [int(data[1])]
							if args.verbosity > 1:
								print("Found multiple non-ref alleles for same sample at position {}, setting allele to -.".format(data[1]))
						gt_dict[int(data[1])] = gt
				else:
					[ref, alt] = data[2:4]
					gt_dict[int(data[1])] = [ref] + [ref if x == "0|0" else alt for x in data[4:]]
		last_pos = min(last_pos + args.query_size, region_end)
		os.remove(temp_path)
		for bad_pos in bad_pos_list:
			del gt_dict[bad_pos]
		bad_pos_list = []
		while len(gt_dict) > args.chunk_size or (last_pos >= region_end and len(gt_dict) > 0):
			fasta_dict = {sample_id: "" for sample_id in sample_ids}
			pos_list = sorted(gt_dict.keys())[0:args.chunk_size]
			pos_map_fname = os.path.join(args.output, "{}_{}-{}".format(region_name, pos_list[0], pos_list[-1]) + "_positions.txt")
			fasta_fname = os.path.join(args.output, "{}_{}-{}".format(region_name, pos_list[0], pos_list[-1]) + ".fasta")
			with open(pos_map_fname, 'w') as file:
				fasta_pos = 0
				for pos in pos_list:
					for [sample_id, sample_idx] in zip(sample_ids, range(0, len(sample_ids))):
						fasta_dict[sample_id] += gt_dict[pos][sample_idx + 1]
					var_coord = "{}:{}".format(region_name, pos)
					file.write("{}_{}\t{}\n".format(os.path.splitext(os.path.basename(fasta_fname))[0], fasta_pos, var_coord))
					del gt_dict[pos]
					fasta_pos += 1
			with open(fasta_fname, 'w') as file:
				for sample_id in sample_ids:
					file.write(">{}\n".format(sample_id))
					file.write("{}\n".format(fasta_dict[sample_id]))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Gene contribution visualizer for ESL pipeline.")
	parser.add_argument("vcf_path", help="VCF file to convert.", type=str)
	parser.add_argument("chrom", help="Chromosome to query.", type=str)
	parser.add_argument("--start", help="First chromosome position to query.", type=int, default=0)
	parser.add_argument("--end", help="Last chromosome position to query.", type=int, default=None)
	parser.add_argument("--output", help="Output fasta file.", type=str, default=None)
	parser.add_argument("--chunk_size", help="Maximum number of positions per FASTA file.", type=int, default=5000)
	parser.add_argument("--query_size", help="Coordinate span when querying VCF.", type=int, default=1000000)
	parser.add_argument("--verbosity", help="Level of verbosity.", type=int, default=1)
	args = parser.parse_args()
	main(args)
