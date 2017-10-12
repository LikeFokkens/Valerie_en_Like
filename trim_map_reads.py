import os, sys, glob, argparse

def init():
	parser = argparse.ArgumentParser(description='Trim and map sequencing reads, using fastq-mcf for trimming (and quality filtering -q 20) and bowtie2 for mapping. Dependencies: fastq-mcf, bowtie2, samtools')
	parser.add_argument('task', type=str, choices = ['all','trim', 'map2transposons', 'extractMappedReads','map2genome'], default = 'all', help='which step has to be performed')
	parser.add_argument('-v', dest='verbose', action = 'store_true', default = True, help='print messages')
	parser.add_argument('-readDir', dest='readDir', type=str, help='name of the directory, it is assumed this fdirectory contains two subdirs -R1 and R2- that contain the .fastq or fastq.gz files with the reads.\
	 Filenames may end with .fastq or .fq')
	parser.add_argument('-outDir', dest='outDir', type=str, help='name of the output directory, output will be put in outDir/01.Trimming/, \
		outDir/02.Mapping/refName (see -refName option), outDir/03a.UnmappedReads/, outDir/03b.MappedReadsWithUnmappedMates/')
	parser.add_argument('-adapterFasta', dest='adapterFasta', type=str, help='fastafile with adapter sequences')
	parser.add_argument('-transposonFasta', dest='transposonFasta', type=str, help='fastafile with the transposon sequences to which we want to map the reads')
	parser.add_argument('-genomeFasta', dest='refName', type=str, help='fastafile with genome to map to')
	parser.add_argument('-nCPU', dest='nCPU', type=str, default = '10', help='number of CPUs to use during mapping')
	parser.add_argument('-I', dest='minInsertSize', type=str, default = '200', help='lower bound of insert size of read library (-I in bowtie2)')
	parser.add_argument('-X', dest='maxInsertSize', type=str, default = '700', help='upper bound of insert size of read library (-X in bowtie2)')
	parser.add_argument('-mp', dest='isMatePair', action = 'store_true', help='if this flag is set, it is assumed that data is mate pair, hence bowtie2 will run with -rf (default it runs with -fr)')
	parser.add_argument('-path2fastq-mcf', dest='path2fastqmcf', type=str, default = '/Applications/ExpressionAnalysis-ea-utils-27a4809/clipper/', help='path to fastq-mcf')
	args = parser.parse_args()

	print '\n\n##################################'
	print '#'
	print '#   SETTINGS'
	print '#'
	if args.verbose:
		argsdict = vars(args)
		for var in argsdict.keys():
			print '# ', var, argsdict[var]
	print '#'
	print '##################################\n\n'
	return args


def trim(trimmer, adapterfasta, dirname_read_data, outdirname, verbose):

	#trimmer = '/Applications/ExpressionAnalysis-ea-utils-27a4809/clipper/fastq-mcf -q 20 '
	cmnd = ''
	failInTotal = 0
	cmnds       = []
	for fastq_R1 in glob.glob(dirname_read_data+'/R1/*'):

		fastq_R2   = fastq_R1.replace('/R1/', '/R2/').replace('_R1','_R2')
		trimmed_R1 = outdirname+fastq_R1.split('/')[-1].replace('.fastq', '.trimmed_fastq-mcf_q20.fastq').replace('.fq', '.trimmed_fastq-mcf_q20.fq')
		trimmed_R2 = outdirname+fastq_R2.split('/')[-1].replace('.fastq', '.trimmed_fastq-mcf_q20.fastq').replace('.fq', '.trimmed_fastq-mcf_q20.fq')

		if verbose:
			print fastq_R1
			print fastq_R2
			print 'paired data'
			print 'trimmed reads in', outdirname+':'
			print trimmed_R1
			print trimmed_R2

		if not os.path.exists(trimmed_R1) or not os.path.exists(trimmed_R2):
			cmnd = trimmer+'fastq-mcf -q 20 -o '+trimmed_R1+' -o '+trimmed_R2+' '+adapterfasta+' '+fastq_R1+' '+fastq_R2 +' >& '+trimmed_R1.replace('_R1', '_')+'.TRIM.LOG'
			cmnds.append(cmnd)
			print cmnd
			fail = os.system(cmnd)
			print fail
			if fail != 0:
				failInTotal += 1
		else:
			print
			print trimmed_R1, 'and', trimmed_R2, 'exist, skipping...'

	print len(cmnds), 'commands, of which', failInTotal,'failed'
	return failInTotal


def mapReads(fasta, dirname_read_data, outdirname, I, X, matePair, nCPU, verbose):

	base_cmnd = 'bowtie2 -x '+fasta+' --end-to-end'
	if matePair:
		base_cmnd += ' --rf'
	else: base_cmnd += ' --fr'

	base_cmnd += ' -p '+nCPU
	base_cmnd += ' --reorder -I '+I+' -X '+X

	failInTotal = 0
	cmnds       = []
	for fastq_R1 in glob.glob(dirname_read_data+'/R1/*'):
		fastq_R2 = fastq_R1.replace('/R1/', '/R2/').replace('_R1','_R2')

		if verbose:
			print "Will map reads from", fastq_R1, 'and', fastq_R2, 'to', fasta

		outfile_base = outdirname+fastq_R1.split('/')[-1].replace('_R1', '_')+'__mapped2__'+fasta.split('/')[-1].split('.fa')[0]+'.bowtie2_I'+I+'X'+X
		if matePair: outfile_base+='-rf'

		if not os.path.exists(outfile_base+'.bam'):
			cmnd = base_cmnd +' -1 '+fastq_R1 +' -2 '+fastq_R2+' -S '+outfile_base+'.sam >& '+outfile_base+'.MAPPING.LOG'
			cmnds.append(cmnd)
			print cmnd
			fail = os.system(cmnd)
			print fail
			if fail == 0:
				cmnd = 'samtools view -S -b '+outfile_base+'.sam > '+outfile_base+'.bam'
				if os.system(cmnd) == 0:
					os.system('rm '+outfile_base+'.sam')
				else:
					failInTotal += 1
			else:
				failInTotal += 1
			
		else:
			print
			print outfile_base+'.bam exists, skipping...'

	print len(cmnds), 'commands, of which', failInTotal,'failed'
	return failInTotal


def extractMappedReads(dirname_read_data, outdirname, verbose):
	failInTotal = 0
	for bam_fname in glob.glob(dirname_read_data+'/*.bam'):
		fail   = 0
		mapped_bams = []
		mapped_bam = outdirname+bam_fname.split('/')[-1].replace('.bam', '.mapped_mateunmapped.bam')
		cmnd   = 'samtools view -b -F 4 -f 8 '+bam_fname+' > '+ mapped_bam 

		if verbose:
			print '\nExtract mapped reads for '+bam_fname+':'
			print cmnd
		fail = os.system(cmnd)
		if fail != 0:
			print 'Command failed, skipping....'
		else: 
			mapped_bams.append(mapped_bam)
			mapped_bam = outdirname+bam_fname.split('/')[-1].replace('.bam', '.unmapped_matemapped.bam')
			cmnd   = 'samtools view -b -F 8 -f 4 '+bam_fname+' > '+ mapped_bam 
			if verbose:
				print '\nExtract mapped reads for '+bam_fname+':'
				print cmnd
			fail = os.system(cmnd)
			if fail != 0:
				print 'Command failed, skipping....'
			else:
				mapped_bams.append(mapped_bam)
				mapped_bam = outdirname+bam_fname.split('/')[-1].replace('.bam', '.both_mapped.bam')
				cmnd   = 'samtools view -b -F 12 '+bam_fname+' > '+ mapped_bam 
				if verbose:
					print '\nExtract mapped reads for '+bam_fname+':'
					print cmnd
				fail = os.system(cmnd)
				if fail != 0:
					print 'Command failed, skipping....'
				else:
					mapped_bams.append(mapped_bam)
		if fail != 0:
			failInTotal += 1
		else:
			print '\n\nSort 3 classes of mapped reads and merge into a single file...'
			for mapped_bam in mapped_bams:
				cmnd = 'samtools sort -n '+mapped_bam+' '+mapped_bam.replace('.bam', '.sorted')
				fail += os.system(cmnd)
			if fail != 0:
				failInTotal += 1
				print cmnd, ' Failed'
			else:
				merged_mapped_bam_fname = outdirname+bam_fname.split('/')[-1].replace('.bam', '.MAPPED.bam')
				cmnd = 'samtools merge -f -n '+ merged_mapped_bam_fname
				for mapped_bam in mapped_bams:
					cmnd += ' '+mapped_bam.replace('.bam', '.sorted.bam')
				cmnd += ' >& '+merged_mapped_bam_fname.replace('.MAPPED.bam', '.MERGE.LOG')

				if verbose: 
					print cmnd
				fail = os.system(cmnd)
				if fail != 0:
					failInTotal += 1
					print cmnd, ' Failed'

				else:
					cmnd = "bedtools bamtofastq -i "+merged_mapped_bam_fname
					cmnd +=' -fq '+merged_mapped_bam_fname.replace('.bam', '_R1.fastq')
					cmnd +=' -fq2 '+merged_mapped_bam_fname.replace('.bam', '_R2.fastq')
					cmnd +=' >& '+merged_mapped_bam_fname.replace('.bam', '.BAM2FASTQ.LOG')
					if verbose: 
						print cmnd
					fail = os.system(cmnd)
					if fail != 0:
						failInTotal += 1
						print cmnd, ' Failed'

	print 'Processed ', len(glob.glob(dirname_read_data+'/*.bam')), 'files, of which', failInTotal, 'failed'
	return failInTotal





if __name__ == "__main__":

	args   = init()
	readDir = args.readDir 
	can_savely_continue = True
	trimOutDir = args.outDir+'/01.Trimming/'

	if args.task == 'trim' or args.task == 'all':
		if args.verbose:
			print 'Will trim reads in', readDir
		
		os.system('mkdir -p '+args.outDir)
		os.system('mkdir -p '+trimOutDir)

		if trim(args.path2fastqmcf, args.adapterFasta, readDir, trimOutDir, verbose = args.verbose) == 0:
			readDir = trimOutDir
		else: can_savely_continue = False

	if can_savely_continue and (args.task == 'map2transposons' or args.task == 'all'):
		if args.verbose:
			print '\n\nWill map reads in', readDir, 'to', args.transposonFasta

		mapOutDir = args.outDir+'/02.MappingToTransposons/'
		os.system('mkdir -p '+args.outDir)
		os.system('mkdir -p '+mapOutDir)
		if len(glob.glob(args.transposonFasta+'*.bt2')) == 0:
			os.system('bowtie2-build '+args.transposonFasta+' '+args.transposonFasta)

		if mapReads(args.transposonFasta, readDir, mapOutDir, args.minInsertSize, args.maxInsertSize, args.isMatePair, args.nCPU, args.verbose) == 0:
			readDir = mapOutDir
		else: can_savely_continue = False

	if can_savely_continue and (args.task == 'extractMappedReads' or args.task == 'all'):
		if args.verbose:
			print '\n\nWill extract mapped reads from bamfiles in ', readDir

		mappedOutDir = args.outDir+'/03.MappedReads/'
		os.system('mkdir -p '+args.outDir)
		os.system('mkdir -p '+mappedOutDir)

		if extractMappedReads(readDir, mappedOutDir, args.verbose) == 0:
			readDir = mappedOutDir
		else: can_savely_continue = False

	if can_savely_continue and (args.task == 'map2genome' or args.task == 'all'):
		if args.verbose:
			print '\n\nWill map reads in', readDir, 'to', args.refName
		mapOutDir = args.outDir+'/04.MappingToGenome/'
		os.system('mkdir -p '+args.outDir)
		os.system('mkdir -p '+mapOutDir)

		if len(glob.glob(args.refName+'*.bt2')) == 0:
			os.system('bowtie2-build '+args.refName+' '+args.refName)

		print mapReads(args.refName, readDir, mapOutDir, args.minInsertSize, args.maxInsertSize, args.isMatePair, args.nCPU, args.verbose)
			