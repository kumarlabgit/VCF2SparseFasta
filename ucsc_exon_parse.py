import sys
import os

# Expected header:
# #bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames

make_bed = True
print_exons = False

wdir = "{}_exons".format(os.path.splitext(sys.argv[1])[0])

if not os.path.exists(wdir):
	os.mkdir(wdir)

with open("exon_parse_log.txt", 'w') as log:
	with open(sys.argv[1], 'r') as file:
		header = file.readline().strip().split("\t")
		name_set = set()
		name2_set = set()
		name_counts = {}
		name2_counts = {}
		coords_counts = {}
		coords_set = set()
		count = 0
		for line in file:
			data = line.strip().split("\t")
			count += 1
			if len(data) != len(header):
				raise Exception("Header/data field count mismatch on line {}.".format(count))
			data_dict = {val[0]:val[1] for val in zip(header, data)}
			if len(data_dict["chrom"]) > 6:
				continue
			for key in ["cdsStart", "cdsEnd", "txStart", "txEnd"]:
				data_dict[key] = int(data_dict[key])
#			if data_dict["name"] in name_set:
#				log.write("Refseq identifier {} observed multiple times.\n".format(data_dict["name"]))
#			if data_dict["name2"] in name2_set:
#				log.write("Gene identifier {} observed multiple times.\n".format(data_dict["name2"]))
			if data_dict["cdsEnd"] - data_dict["cdsStart"] < 1 or data_dict["txEnd"] - data_dict["txStart"] < 1:
				log.write("Skipping {} because it has a coding span of 0.\n".format(data_dict["name"]))
				continue
			name_set.add(data_dict["name"])
			name2_set.add(data_dict["name2"])
			name_counts[data_dict["name"]] = name_counts.get(data_dict["name"], 0) + 1
			name2_counts[data_dict["name2"]] = name2_counts.get(data_dict["name2"], 0) + 1
			exon_starts = [int(val) for val in data_dict["exonStarts"].split(",") if val != '']
			exon_ends = [int(val) for val in data_dict["exonEnds"].split(",") if val != '']
			tx_start = data_dict["txStart"]
			tx_end = data_dict["txEnd"]
			cds_start = data_dict["cdsStart"]
			cds_end = data_dict["cdsEnd"]
			if len(exon_starts) != len(exon_ends):
				raise Exception("Exon start/end count mismatch on line {}.".format(count))
			exon_count = 0
			with open(os.path.join(wdir, "{}.bed".format(data_dict["name2"])), 'w') as exon_file:
				for [exon_start, exon_end] in zip(exon_starts, exon_ends):
					if exon_end < cds_start or exon_end < cds_start or exon_start > cds_end or exon_start > tx_end:
						continue
					exon_count += 1
					exon_start = max(exon_start, cds_start, tx_start)
					exon_end = min(exon_end, cds_end, tx_end)
					coord_string = "{}:{}-{}".format(data_dict["chrom"], exon_start, exon_end)
	#				if coord_string in coords_set:
	#					log.write("Coordinate span {} observed multiple times.\n".format(coord_string))
					coords_counts[coord_string] = coords_counts.get(coord_string, 0) + 1
					if print_exons:
						print("{}-{}\t{}-{}\t{}\t{}\t{}".format(data_dict["name2"], exon_count, data_dict["name"], exon_count, data_dict["chrom"], exon_start, exon_end))
					if make_bed:
						exon_file.write("{}\t{}\t{}\n".format(data_dict["chrom"], exon_start, exon_end))
			if not make_bed:
				os.remove(os.path.join(wdir, "{}.bed".format(data_dict["name2"])))
		for key in name_counts.keys():
			if name_counts[key] > 1:
				log.write("Refseq ID {} observed {} times\n".format(key, name_counts[key]))
		for key in name2_counts.keys():
			if name2_counts[key] > 1:
				log.write("Gene symbol {} observed {} times\n".format(key, name2_counts[key]))
		for key in coords_counts.keys():
			if coords_counts[key] > 1:
				log.write("Coordinate range {} observed {} times\n".format(key, coords_counts[key]))
