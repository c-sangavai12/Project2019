import vcf
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd

vcf_variants = defaultdict(set)
vcf_positions = defaultdict(set)
#vcf_features = defaultdict(set)
#feature = {}

def parseVCF(vcffile):
	vcfreader = vcf.Reader(open(vcffile))
	samples = vcfreader.samples
	sampleName = None
	if len(samples) == 0:
		print("No sample found, consolidating it as known_mutations")
		sampleName = vcffile.split("/")[-1].replace(".vcf", "") + "_known_mutations"

	for record in vcfreader:
		alts = record.ALT
		chrm = record.CHROM.replace('chr', '')
		chrm = chrm.replace('Chr', '')
		info = record.INFO
		pos = chrm + "_" + str(record.POS)
		if sampleName and "_known_mutations" in sampleName:
				vcf_positions[sampleName].add(pos)
				alts = record.ALT
				#feature[sampleName]
				for alt in alts:
					vcf_variants[sampleName].add("_".join([pos, record.REF, str(alt)]))

		else:
			for sample in samples:
				sampleName = vcffile.split("/")[-1].replace(".vcf", "") + "_" + sample
				vcf_positions[sampleName].add(pos)
				try:
				    alts = record.genotype(sample).gt_bases.split("|")
				except:
					continue
				if len(alts) == 1:
					alts = record.genotype(sample).gt_bases.split("/")
				for alt in alts:
					vcf_variants[sampleName].add("_".join([pos, record.REF, str(alt)]))
	print(sampleName, "=" , len(vcf_variants[sampleName]))

def parseCommon(allcommon):
	vcfreader = vcf.Reader(open(allcommon))
	vcfreader._parse_info
	commonMutations = set()
	for record in vcfreader: 
		for alt in record.ALT:
			chrm = record.CHROM.replace('chr', '')
			chrm = chrm.replace('Chr', '')
			identifier = "_".join([chrm, str(record.POS), record.REF, str(alt)])
			commonMutations.add(identifier)
				
	return commonMutations

def non_common_mutations(vcf_positions, vcf_variants, commonMutations):
	samples = vcf_positions.keys()
	print("SAMPLES being compared:", samples)

	non_common = {}
	for sample1 in samples:
		non_common[sample1] = (vcf_variants[sample1] - (commonMutations))
		
	for sample1 in samples:
		fout = open(sample1 + "_filteredVariants.txt", 'w')
		for variant in non_common[sample1]:
			fout.writelines(variant + "\n")
		fout.close()
		print(len(non_common[sample1]))
	return non_common
	


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--vcfs", action="append", help="list vcfs")
    parser.add_argument("-ac", "--all_common", help="all common vcf")
    parser.add_argument("-o", "--output_path", help="output path")

    args = parser.parse_args()
    for eachVcf in args.vcfs:
    	print("PARSING", eachVcf)
    	parseVCF(eachVcf)

    commonMutations = parseCommon(args.all_common)
    non_common = non_common_mutations(vcf_positions, vcf_variants, commonMutations)
    #vcf_exclusives = unique_mutations(vcf_variants)

if __name__ == "__main__":
    main()

	## Common variants between every 2 samples and with cosmic:
	


