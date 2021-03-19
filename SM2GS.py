#!/usr/bin/python

import argparse as ap
import sys
import re
from scipy.stats import chisquare
from itertools import compress
gens_d = {'0/0':0, '0/1':1, '1/1':2, '1/2':3, '2/2':4, '0/0/0/0':0, '0/0/0/1':1, '0/0/1/1':2, '0/1/1/1':3, '1/1/1/1':4, '1/1/1/2':5, '1/1/2/2':6, '1/2/2/2':7, '2/2/2/2':8}
gens_tassel = {'0/0':0, '0/1':1, '1/1':2, '0/0/0/0':0, '0/0/0/1':1, '0/0/1/1':2, '0/1/1/1':3, '1/1/1/1':4, '1/1/1/2':5, '1/1/2/2':6, '1/2/2/2':7, '2/2/2/2':8}
gens_pa = {'0/0':0, '0/1':1, '1/1':1, '0/0/0/0':0, '0/0/0/1':1, '0/0/1/1':1, '0/1/1/1':1, '1/1/1/1':1}

par = ap.ArgumentParser(description = 'Convert from VCF format to the SNP Loc format expected by TetraploidSNPMap. Variant sites with more than two alleles, along with variant sites for which one or both parents has an unknown genotype, will be removed because TetraploidSNPMap can\'t handle such sites.')
par.add_argument('-v', '--vcf', metavar = 'input.vcf', required = True, help = 'VCF file to process.', dest = 'v')
par.add_argument('-o', '--output', metavar = 'output.snploc', required = True, help = 'SnpLoc file to output.', dest = 'o')
par.add_argument('-p1', metavar = 'parent1', required = False, help = 'Sample name of first parent (must match a column in the VCF header line).')
par.add_argument('-p2', metavar = 'parent2', required = False, help = 'Sample name of second parent (must match a column in the VCF header line).')
par.add_argument('-f', '--offset', metavar = 'offset', required = False, help = 'Column of first genotype, numbered from 0. Default 9.', dest = 'f')
par.add_argument('-d', '--fix-double-reductions', required = False, help = 'Look for double-reducion genotypes and mark as unknown instead of the actual genotype. TSM can\'t analyze loci with double-reducitons and will refuse any locus that has any.', dest='dr', action = "store_true")
par.add_argument('-z', '--remove-homozygotic-parents', required = False, help = 'Remove loci where the parents are identical homozygotes (AAAAxAAAA or BBBBxBBBB). These loci are not informative for linkage mapping.', dest='rh', action = "store_true")
par.add_argument('-gq', '--min-gen-qscore', metavar = 'phred-qscore', required = False, help = 'Mark as unknown any genotypes with a genotype phred q-score (GQ) less than phred-score', dest = 'q')
par.add_argument('-vq', '--min-var-qscore', metavar = 'phred-qscore', required = False, help = 'Remove any variant sites with a locus phred q-score (GQ) less than phred-score', dest = 'vq')
par.add_argument('-c', '--chromosome', metavar = 'chromosome', required = False, help = 'Output SNPs on the specified chromosome', dest = 'c')
par.add_argument('-r', '--min-genotype-reads', metavar = 'M', required = False, help = 'Drop (mark as 9, unknown) genotypes supported by fewer than N total reads.', dest='r')
par.add_argument('-u', '--max-unknown-genotypes', metavar='f', required = False, help = 'Discard loci for which the fraction of unknown genotypes is more than f. The value f should be a decimal number between 0 and 1. Default is not to do any such filtering (equivalent to f = 1.0). Parents do not count toward this total (but if either is unknown then the locus will be discarded regardless of this setting).', dest='u')
par.add_argument('-a', '--presence-absence', required = False, help = 'Produce binary, presence/absence output only, with 1 = alternative allele is present in the genotype, 0 = it is not present, 9 = unknown', dest='a', action = "store_true")
par.add_argument('-p', '--min-chisq-p-val', metavar='cutoff', required = False, help = 'Cutoff for chi-square significance (p-value). Smaller p means less chance offspring ratios are good. Loci with p < cutoff will be rejected. Default is 0.001, the same as TetraploidSNPMap. Calculated after any adjustments that reject individual genotypes such as -d, -r, and -gq. Setting cutoff = 0.0 disables this filtering.', dest = 'p')
par.add_argument('-m', '--minor-allele-freq', metavar = 'freq', required = False, help = 'Minimum minor allele frequency accross all genotypes at a locus. Loci wherein either allele is less frequent than freq will be dropped. Must be a number between 0.0 and 1.0. Each allele at a locus counts separately, so e.g. three genotypes, AAAA ABBB AABB, would have frequency 5/12 = .42 for allele B.', dest = 'm')
par.add_argument('-t', '--tassel-output', required = False, help = 'Alternative output format for TASSEL 5. Header line with sample names is included instead of matrix dimensions, genotype signifiers up to 8 for 2/2/2/2 are included (with allelic series no longer filtered out) and unknowns are represented by -1 rather than 9. Cannot be used with chi-square filtering (-p)', dest = 't', action = "store_true")
par.add_argument('-us', '--max-unknown-sample-genotpyes', metavar = 'f', required = False, help = 'Discard samples for which the fraction of unknown genotypes is greater than f. The value f should be a decimal number between 0 and 1. If parents are specified, this filtering will not be applied to them.', dest = 'us')
par.add_argument('--debug', action = 'store_true', help = 'Show debug messages', dest = 'debug')
newline = '\n'
# Expected offspring ratios for each possible cross. Used for chi-square
# filtering (-p). 3-dimensional array with first index = parent 1 dosage,
# second index = parent 2 dosage, 3rd index = offspring dosage, value =
# fraction of offspring of those parents expected to have that dosage.
expected =    [[[1.0, 0.0, 0.0, 0.0, 0.0],				# AAAAxAAAA
				[1.0/2.0,1.0/2.0,0.0,0.0,0.0],			# AAAAxAAAB
				[1.0/6.0,4.0/6.0,1.0/6.0,0.0,0.0],		# AAAAxAABB
				[0.0,1.0/2.0,1.0/2.0,0.0,0.0],			# AAAAxABBB
				[0.0,0.0,1.0,0.0,0.0]],					# AAAAxBBBB

			   [[1.0/2.0,1.0/2.0,0.0,0.0,0.0],				# AAABxAAAA
			    [1.0/4.0,2.0/4.0,1.0/4.0,0.0,0.0],          # AAABxAAAB
				[1.0/12.0,5.0/12.0,5.0/12.0,1.0/12.0,0.0],  # AAABxAABB
				[0.0,1.0/4.0,2.0/4.0,1.0/4.0,0.0],          # AAABxABBB
				[0.0,0.0,1.0/2.0,1.0/2.0,0.0]],             # AAABxBBBB

			   [[1.0/6.0,4.0/6.0,1.0/6.0,0.0,0.0],				# AABBxAAAA
				[1.0/12.0,5.0/12.0,5.0/12.0,1.0/12.0,0.0],      # AABBxAAAB
			    [1.0/36.0,8.0/36.0,18.0/36.0,8.0/36.0,1.0/36.0],# AABBxAABB
				[0.0,1.0/12.0,5.0/12.0,5.0/12.0,1.0/12.0],      # AABBxABBB
				[0.0,0.0,1.0/6.0,4.0/6.0,1.0/6.0]],             # AABBxBBBB

			   [[0.0,1.0/2.0,1.0/2.0,0.0,0.0],				# ABBBxAAAA
				[0.0,1.0/4.0,2.0/4.0,1.0/4.0,0.0],          # ABBBxAAAB
				[0.0,1.0/12.0,5.0/12.0,5.0/12.0,1.0/12.0],  # ABBBxAABB
				[0.0,0.0,1.0/4.0,2.0/4.0,1.0/4.0],          # ABBBxABBB
				[0.0,0.0,0.0,1.0/2.0,1.0/2.0]],             # ABBBxBBBB

			   [[0.0,0.0,1.0,0.0,0.0],					# BBBBxAAAA
				[0.0,0.0,1.0/2.0,1.0/2.0,0.0],          # BBBBxAAAB
				[0.0,0.0,1.0/6.0,4.0/6.0,1.0/6.0],      # BBBBxAABB
				[0.0,0.0,0.0,1.0/2.0,1.0/2.0],          # BBBBxABBB
				[0.0,0.0,0.0,0.0,1.0]]]                 # BBBBxBBBB

args = par.parse_args()
def dbg_msg(message):
	if(args.debug):
		sys.stderr.write(message+'\n')
have_parents = args.p1 and args.p2
if not have_parents:
	if not args.t:
		sys.stderr.write('Parents must be specified for TetraploidSNPMap output.')
		exit(1)
	if args.dr:
		sys.stderr.write('Fixing double reductions (-d) requires both parents be specified (-p1 and -p2)')
		exit(1)
	if args.rh:
		sys.stderr.write('Removing heterozygous parents (-z) requires both parents be specified (-p1 and -p2)')
		exit(1)
	if args.p:
		sys.stderr.write('Filtering by p-value (-p) requires both parents be specified (-p1 and -p2)')
		exit(1)
if(args.f):
	offs = int(args.f)
else:
	offs = 9
if(args.u):
	u = float(args.u)
else:
	u = 0.0
if(args.a):
	gens = gens_pa
else:
	gens = gens_d
if(args.p):
	p_cutoff = float(args.p)
else:
	p_cutoff = 0.001
if(args.m):
	m = float(args.m)
	neg_m = 1.0 - m
else:
	m = 0.0
if(args.t):
	unknown = 'NA'
	separator = '\t'
	if(args.p):
		sys.stderr.write("Error: Cannot use p-value filtering (-p) with tassel output format (-t)")
		exit(1)
	if(args.a):
		sys.stderr.write("Error: Cannot use presence-absence output (-a) with tassel output format (-t)")
		exit(1)
else:
	separator = ' '
	unknown = 9
try:
	infile = open(args.v)
except IOError as e:
	sys.stderr.write("Could not open input file '%s': %s.\n"%(args.v,e.strerror))
	exit(1)

try:
	outfile = open(args.o,'w')
except IOError as e:
	sys.stderr.write("Could not open output file '%s' for writing: %s.\n"%(args.o,e.strerror))
	exit(1)

# In this section, we read the VCF header and find what columns the two parents
# are in, stored as p1 and p2 (these are stored as their ordering among the
# offspring, so minus the fixed, non-genotype columns as specified by --offset)
while True:
	headerline = infile.readline()
	if (not headerline):	# Ran out of file and haven't found a header line. Exit.
		sys.stderr.write("Error: VCF file %s has no header line.\n"%args.v)
		exit(1)
	if (headerline[0:6] == '#CHROM'): break	# Found the header line, continue on
	if (headerline[0] != '#'):	# Encountered data and haven't found a header line. Exit.
		sys.stderr.write("Error: VCF file %s has no header line.\n"%args.v)
		exit(1)
p1 = -1
p2 = -1
fields = headerline.split()
n_fields = len(fields)
num_gt = n_fields - offs
if(have_parents):
	for i in range(offs,n_fields):
		if(fields[i] == args.p1):
			p1 = int(i-offs)
		if(fields[i] == args.p2):
			p2 = int(i-offs)
	if(p1 < 0):
		sys.stderr.write("Error: Could not find parent %s in header line.\n"%args.p1)
		exit(1)
	if(p2 < 0):
		sys.stderr.write("Error: Could not find parent %s in header line.\n"%args.p2)
		exit(1)
# pmin and pmax are here so we can make sure we remove the parents in the
# reverse order they appear in the list of samples, regardless of whether p1 or
# p2 is first. If we removed the first one first then the list would shift over
# and we'd end up removing the wrong one for the second parent.
	pmin = min(p1,p2)
	pmax = max(p1,p2)
	p1_name = fields[p1+offs]
	p2_name = fields[p2+offs]
	fields.pop(pmax)
	fields.pop(pmin)
	samples = [p1_name,p2_name]+fields[offs:]
	#tassel_header = separator.join(["Chr_loc", p1_name, p2_name] + fields[offs:] + ["Alleles"])	# Output header for tassel format (-t) if parents specified
else:
	samples = fields[offs:]
	#tassel_header = separator.join(["Chr_loc"] + fields[offs:] + ["Alleles"])	# Output header for tassel format (-t) if parents unspecified

all_dosages = []	# 2D array that contains all genotype dosages in the file. First index is locus, second is sample. Created progressively by appending the array of locus genotypes to it for each locus in the input file (i.e. for each loop iteration)
alleles_list = []	# List of strings describing alleles at each locus
names_list = []		# List of locus names
unknown_by_sample = [0] * (num_gt+1)	# Number of unknown genotypes for each sample. Used in -us filtering.
snp_n = 0
havent_warned_no_qscore = True
havent_warned_no_reads = True
qscore_pos = -1
reads_pos = -1
if(args.q): q = float(args.q)
if(args.r): r = int(args.r)
if(args.vq): vq = float(args.vq)

# outstr will hold the entire contents of the snploc file we're going to write.
# We cache everything in this string and write it out at the end so that we can
# put an accurate count of SNPs in the snploc header line, since we won't know
# until the end how many SNPs will end up making it through filtering.
outstr = ''

# Counts of loci and genotypes rejected for various reasons, used to output the
# statistics at the end
rejected_chr = 0
rejected_alleles = 0
rejected_monomorphic = 0
rejected_unk_parent = 0
rejected_homozygous = 0
rejected_vqscore = 0
rejected_unknown = 0
rejected_chisq = 0
rejected_maf = 0

rejected_qscore = 0
rejected_reads = 0
rejected_dr = 0

retained = 0


# The main loop, one iteration per line in the VCF file
while True:
	line = infile.readline()
	if (not line): break	# Reached EOF
	fields = line.split()
	chromosome = fields[0]
	if (args.c and chromosome != args.c):
		rejected_chr  += 1
		continue
	unknown_gt = 0	# Number of unknown genotypes at this locus. Used for filtering by fraction of unknown genotypes (-u option)
	loc = fields[1]
	# Name the locus with chromosome and physical map position, so we can
	# easily compare the physical vs linkage map afterward
	name = '%s_%s'%(chromosome,loc)
	try:
		if args.vq and float(fields[5]) < vq:
			rejected_vqscore += 1
			continue
	except ValueError:
		pass
	dosage = [None]*(n_fields-offs)

	# Find which field within the genotypes has the q-score (GQ)
	if (args.q and qscore_pos < 0 and havent_warned_no_qscore):
		gen_field_names = fields[offs-1].split(':')
		for i in range(0,len(gen_field_names)):
			if (gen_field_names[i] == 'GQ'):
				qscore_pos = i
				break
		if (qscore_pos < 0 and havent_warned_no_qscore):
			sys.stderr.write("Warning: Filter by genotype q-score requested, but one or more loci do not have genotype q-scores.\n")
			havent_warned_no_qscore = False
	# Find which field within the genotypes has the reads (DP)
	if (args.r and reads_pos < 0 and havent_warned_no_reads):
		gen_field_names = fields[offs-1].split(':')
		for i in range(0,len(gen_field_names)):
			if (gen_field_names[i] == 'DP'):
				reads_pos = i
				break
		if (reads_pos < 0 and havent_warned_no_reads):
			sys.stderr.write("Warning: Filter by genotype read depth requested, but one or more loci do not have genotype read depth.\n")
			havent_warned_no_reads = False
	try:
		# Iterate over the genotypes at this locus
		ploidy = 0
		for g in range(offs,n_fields):
			gen_fields = fields[g].split(':')
			ploidy = max(ploidy, len(gen_fields[0].split('/')))
			if (gen_fields[0] == '.'):
				dosage[g-offs] = unknown
				unknown_gt += 1
				continue
			# Filter by qscore if requested (a nonnegative qscore_pos must mean
			# it was requested, so no need to check again explicitly)
			if qscore_pos >= 0 and float(gen_fields[qscore_pos]) < q:
				rejected_qscore  += 1
				dosage[g-offs] = unknown
				unknown_gt += 1
				continue
			# Filter by read depth if requested (a nonnegative reads_pos must
			# mean it was requested, so no need to check again explicitly)
			if reads_pos >= 0 and int(gen_fields[reads_pos]) < r:
				rejected_reads  += 1
				dosage[g-offs] = unknown
				unknown_gt += 1
				continue
			# Look up the genotype number (will throw an exception if it isn't
			# found, caught below)
			dosage[g-offs] = gens[gen_fields[0]]
		# Rescue non-reference genotypes. If no genotypes have the reference
		# (allele 0), then we may still have a useable locus even if it at
		# first appears to be an allelic series. Try to fix it by moving all
		# alleles "down" by 1 (i.e. 1 -> 0, 2 -> 1)
		ref_absent = True
		allelic_series = False
		# Go through and determine whether there are no reference alleles and
		# whether this is an allelic series
		for i in range(0,len(dosage)):
			if(dosage[i] != unknown):
				if(dosage[i] <= ploidy-1):
					ref_absent = False
				if(dosage[i] > ploidy):
					allelic_series = True
		if (ref_absent):
			# No reference alleles, go through and apply the correction
			for i in range(0,len(dosage)):
				if(dosage[i] != unknown):
					dosage[i] = dosage[i]-ploidy
		elif (allelic_series):
			# There was a reference allele 0 somewhere (ref_absent is false)
			# and there's a 2 somewhere (allelic_series is true) and if there's
			# a 2 then there's also a 1. It's therefore a real allelic series,
			# and TPM can't handle those. Reject.
			dbg_msg('Rejected '+name+' due to series logic')
			rejected_alleles += 1
			continue
	except KeyError as e:
		# Found a locus with >3 alleles. Can't rescue these at all. Reject.
		dbg_msg('Rejected '+name+' due to KeyError')
		rejected_alleles += 1
		continue
	# If either parent is unknown, TPM won't be able to do anything with this
	# snp. Reject.
	if have_parents and (dosage[p1] == unknown or dosage[p2] == unknown):
		rejected_unk_parent  += 1
		continue

	# Fix (some) double reductions if requested by the user; if parents are
	# AAAB (1) x AAAA (0) then mark any dosages in the offspring other than
	# AAAB and AAAA as unknown.  Also do the equivalent for the reverse
	# (ABBB (3) x BBBB (4)).
	if (args.dr):
		if (ploidy != 4):
			sys.stderr.write('Error: fixing double reductions only works with ploidy 4, but this file has genotypes of ploidy %i. Genotype is %s\n'%(ploidy,gen_fields[0]))
			exit(1)
		if ((dosage[p1] == 0 and dosage [p2] == 1) or (dosage[p1] == 1 and dosage[p2] == 0)):
			for i in range(0,len(dosage)):
				if(dosage[i] != 1 and dosage[i] != 0):
					rejected_dr += 1
					dosage[i] = unknown
					unknown_gt += 1
		if ((dosage[p1] == 3 and dosage [p2] == 4) or (dosage[p1] == 4 and dosage[p2] == 3)):
			for i in range(0,len(dosage)):
				if(dosage[i] != 4 and dosage[i] != 3):
					rejected_dr += 1
					dosage[i] = unknown
					unknown_gt += 1
	if(args.u):
		if(float(unknown_gt) / num_gt > u):
			rejected_unknown += 1
			continue
	for i in range(num_gt):
		if(dosage[i] == unknown):
			unknown_by_sample[i+1] += 1
	# Remove AAAA x AAAA, BBBB x BBBB, and AAAA x BBBB loci if requested by the user
	if (args.rh):
		if((dosage[p1] == 0 or dosage[p1] == 4) and (dosage[p2] == 4 or dosage[p2] == 0)):
			rejected_homozygous += 1
			continue
	if(have_parents):
		p1_str = dosage[p1]
		p2_str = dosage[p2]
		dosage.pop(pmax)
		dosage.pop(pmin)

	# Now we will calculate the p-value for difference between the observed
	# offspring ratios at this locus and what we'd expect based on the parents'
	# genotypes. A low p-value means less chance that observed deviation from
	# expected ratios could be due to chance rather than an actual problem with
	# the calls. TetraploidSNPMap rejects loci where p < 0.001, so that is our
	# cutoff unless overridden by the user. In these calculations, we ignore
	# offspring genotypes for which the expected number is 0 (which would cause
	# a divide by zero error) because this is what TPM does and we're trying to
	# replicate their p-value calculation.
	if(args.p):
		if (ploidy !=4 and ploidy != 1):
			sys.stderr.write('Error: p-value filtering is only supported for ploidy 4 at this time; ploidy is %i.\n'%ploidy)
			exit(1)
		# Iterate over offspring and get totals of each genotype
		observed_offspr = [0,0,0,0,0]
		expected_offspr = expected[p1_str][p2_str]
		known_offspr = 0
		for offspr in dosage:
			offspr = int(offspr)
			if(offspr <= 4):
				observed_offspr[offspr] += 1
				known_offspr += 1

		# Iterate over the five possible genotypes and keep only those for
		# which the expected number is nonzero, since chi-squared will divide
		# by it
		expected_nz = []
		observed_nz = []
		for i in range(0,len(expected_offspr)):
			if(expected_offspr[i] != 0):
				expected_nz  += [expected_offspr[i]]
				observed_nz  += [observed_offspr[i]]
		# Calculate the chi-squared and p value and reject this locus if the
		# latter is below the cutoff
		(chisq,p) = chisquare(observed_nz,[e*known_offspr for e in expected_nz])
		if(p < p_cutoff):
			rejected_chisq += 1
			continue
	# Filter by minor allele frequency > m
	n_alt_allele = 0
	known_alleles = 0
	for offspr in dosage:
		if(offspr != unknown):
			n_alt_allele += offspr
			known_alleles += 4
	maf = float(n_alt_allele) / known_alleles
	if args.m is not None and (maf <= m or maf >= neg_m):
		rejected_maf += 1
		continue
	snp_n += 1
	all_dosages += [[name]+dosage]

	if args.t:
		alleles_list += [separator+fields[3] + ',' + fields[4]]
	else:
		alleles_list += [""]

samples = ['<Marker>']+samples
rejected_sample_unknowns = 0
if(args.us is not None):	# Filtering samples by fraction unknown genotypes
	if have_parents:	# We don't want to filter out either parent even if they fail this test because the whole file would be rendered useless if we did. If there are parents, then we want to remove their entries from the list of unknown counts and replace them with 0's tacked onto the beginning, to match the configuration of the all_dosages arrays. The 0's garauntee that they'll be kept no mater what the cutoff is.
		unknown_by_sample.pop(pmax)
		unknown_by_sample.pop(pmin)
		unknown_by_sample = [0,0]+unknown_by_sample
	# Here we're going to do the filtering. First calculate a threshold telling us how many unknowns are allowable for the given number of snps (faster then dividing each time), then generate a list of booleans telling us which samples are under this limit and should therefore be kept (samples_to_keep), then apply that to each entry in all_dosages using the compress() funciton
	keep_sample_threshold = float(args.us) * snp_n
	dbg_msg('Unknowns by sample: %s'%unknown_by_sample)
	dbg_msg('Keep threshold: %s'%keep_sample_threshold)
	samples_to_keep = [sample_unknowns <= keep_sample_threshold for sample_unknowns in unknown_by_sample]
	rejected_sample_unknowns = len(samples)-sum(samples_to_keep)
	dbg_msg('Keeping samples: %s'%samples_to_keep)
	all_dosages = [list(compress(dosages_for_locus, samples_to_keep)) for dosages_for_locus in all_dosages]
	samples = list(compress(samples, samples_to_keep))
if(args.t):
	all_dosages = [samples] + all_dosages	# Add the header for the TASSEL format
	all_dosages = map(list, zip(*all_dosages))	# Transpose the whole things since TASSEL wants samples in rows and markers in columns, the opposite arrangement from that of VCF and TPM
else:	# If we're outputting for TPM, we need to truncate to 300 samples, since TPM won't accept a file with more than this
	all_dosages = [dosages_for_locus[0:min(len(dosages_for_locus),301)] for dosages_for_locus in all_dosages]
	outfile.write(('%i %i'+newline)%(n_fields-11,snp_n))	# The header for the TPM format

for l in all_dosages:
	outfile.write(separator.join([str(d) for d in l])+newline)
#outfile.write(outstr)
infile.close()
outfile.close()
sys.stderr.write("Retained Loci: %i\nRejected Loci:\n\tWrong Chromosome: %i\n\tMore than Two Alleles: %i\n\tUnknown Parent(s): %i\n\tHomozygous Parents: %i\n\tLow Q-score: %i\n\tP < %f: %i\n\tMinor Allele Frequency <= %f: %i\n\tMore than %i%% Unknown Genotypes: %i\nRejected Genotypes:\n\tLow Q-Score: %i\n\tDouble Reduction: %i\n\tRead Count: %i\nRejected Samples:\n\tUnknown Genotypes: %i\n"%(snp_n, rejected_chr, rejected_alleles, rejected_unk_parent, rejected_homozygous, rejected_vqscore, p_cutoff, rejected_chisq, m, rejected_maf, int(100*u), rejected_unknown, rejected_qscore, rejected_dr, rejected_reads, rejected_sample_unknowns))
if(snp_n > 8000) and not args.t:
	sys.stderr.write("Warning: Number of retained loci is %i. As of March 2017, the current version of TetraploidSNPMap cannot handle more than 8,000 loci. This output file may be rejected by the software."%snp_n)
