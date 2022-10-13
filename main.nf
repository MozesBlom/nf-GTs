#!/usr/bin/env nextflow

log.info """\
		From VCF to GTs - Calculate GT distributions  
		===================================
		P1_file    		: ${params.indivs_P1_file}
		P2_file    		: ${params.indivs_P2_file}
		P3_file    		: ${params.indivs_P3_file}
		O_file			: ${params.indivs_outgroup_file}
		chromos      	: ${params.chromos_file}
		sex_chromos  	: ${params.sex_chromos}
		chromo_coords  	: ${params.chromo_coords}
		reference    	: ${params.ref_file}
		vcfdir       	: ${params.inputdir}
		outdir       	: ${params.outputdir}

		Window sizes 	: ${params.window_sizes}
		"""
		.stripIndent()

/*
 * Create five value channels :
 * - P1
 * - P2
 * - P3
 * - O
 *
 * - An aggregate list with all individuals (for which consensus sequences need to be called)
 */

indivs_P1 = channel.value(file(params.indivs_P1_file).readLines())
indivs_P2 = channel.value(file(params.indivs_P2_file).readLines())
indivs_P3 = channel.value(file(params.indivs_P3_file).readLines())
indivs_O = channel.value(file(params.indivs_outgroup_file).readLines())

indivs_all = indivs_P1.concat(indivs_P2, indivs_P3, indivs_O).flatten()

/*
 * Gene-trees for which chromosomes to include?
 */
chromos = file(params.chromos_file).readLines()


/*
 * Value channel for sex chromosomes. For cross compatibility with downstream plotting scripts
 */
sex_chromos = channel.value(params.sex_chromos)

process call_consensus {
    
/*
 * Call consensus sequences for each individual and each chromosome.
 *
 * Directed into four concensus channels that are used to create cartesian products of each combination
 *
 */

	tag "Create consensus sequence for $chromo of $indiv"

	input:
	val(indiv) from indivs_all
	each chromo from chromos

	output:
	tuple val(indiv), \
	val(chromo), \
	path("${indiv}_${chromo}_cons.fa") into consensus_ch1, consensus_ch2, consensus_ch3, consensus_ch4

	when:
	params.call_consensus == true

	script:
	def chromo_var_fn = file("${params.inputdir}/${indiv}/${chromo}_vars_filt_indels.vcf.gz")
	def chromo_mask_fn = file("${params.inputdir}/${indiv}/${chromo}_mask_cov_het.vcf")
	"""

	samtools faidx ${params.ref_file} ${chromo} | bcftools consensus ${chromo_var_fn} -m ${chromo_mask_fn} -o ${indiv}_${chromo}_cons.fa

	sed -i 's/${chromo}/${indiv}/g' ${indiv}_${chromo}_cons.fa

	"""
}

/*
 * The aim is to create cartesian product combinations for each P1 x P2 x P3 x O
 * A multiple sequence alignment of four individuals is created for each set of 4.
 * This is done for each combination and each chromosome.
 *
 * First create a channel where only consensus sequences for each P1/P2/P3/O is stored
 * Will be of format: indiv, chromo, consensus
 *
 */

indivs_P1
	.flatten()
	.combine(consensus_ch1, by:0)
	.set { consensus_ch_p1 }

indivs_P2
	.flatten()
	.combine(consensus_ch2, by:0)
	.set { consensus_ch_p2 }

indivs_P3
	.flatten()
	.combine(consensus_ch3, by:0)
	.set { consensus_ch_p3 }

indivs_O
	.flatten()
	.combine(consensus_ch4, by:0)
	.set { consensus_ch_O }

/*
 * Now make cartesian combinations P1 x P2 x P3 x O
 * Will be of format: chromo, indiv_p1, consensus_P1, indiv_p2, consensus_P2, indiv_p3, consensus_P3, indiv_O, consensus_O
 */
consensus_ch_p1
    .combine(consensus_ch_p2, by: 1)
	.map{ chromo, p1, p1_cons, p2, p2_cons -> tuple(p1, chromo, p1_cons, p2, p2_cons)  }
	.combine(consensus_ch_p3, by: 1)
	.map{ chromo, p1, p1_cons, p2, p2_cons, p3, p3_cons -> tuple(p1, chromo, p1_cons, p2, p2_cons, p3, p3_cons)  }
	.combine(consensus_ch_O, by: 1)
    .set { consensus_combos_ch }



/*
 * For each combination P1 x P2 x P3 x O, and each chromosome: 
 * Create a multiple sequence alignment
 */

process chromo_raw_msa {

	/*
	 * Create multiple sequence alignment for each chromosome
	 * Directed into one channel
	 * 1) For process (subset_MSA_to_windows) -- Subsetting each chromosome into windows
	 */

	tag "Raw multiple sequence alignment per chromosome -- raw"

	input:
	tuple val(chromo), \
	val(p1), \
	path(p1_cons_fn), \
	val(p2), \
	path(p2_cons_fn), \
	val(p3), \
	path(p3_cons_fn), \
	val(out), \
	path(out_cons_fn) from consensus_combos_ch

	output:
	tuple val(chromo), \
	val(p1), \
	val(p2), \
	val(p3), \
	val(out), \
	path("${p1}_${p2}_${p3}_${out}_${chromo}_raw_msa.fa") into combo_msa_by_chromo_raw_ch1

	script:
	"""

	cat ${p1_cons_fn} ${p2_cons_fn} ${p3_cons_fn} ${out_cons_fn} > ${p1}_${p2}_${p3}_${out}_${chromo}_raw_msa.fa

	"""
}


/*
 * Consensus sequences have now been called and evaluated, and raw sequence alignments have been created for each combination and chromosome.
 *
 * We will now subset each genome, if required for a range of subset sizes, and filter the resulting window alignments for missing data.
 *
 */


process subset_MSA_to_windows {

	/*
	 * Subset MSA for each chromosome in windows of size X and minimally spaced Y bp. from eachother
	 * Each candidate window is subsequently inspected for missing data: a) Per column and b) per individual
	 * a) Per column: params.missing_data_per_col_threshold, params.missing_cols_per_window_threshold
	 * b) Per individual: params.missing_data_per_row_threshold, params.missing_indivs_per_window_threshold
	 *
	 * If a given window surpasses the filtering criteria specified above. It will generate two alignments for the same window:
	 * "*all_indivs_aln.fa" = This includes all individuals, even if some individuals may not fulfill (missing_data_per_row_threshold)
	 * "*filt_indivs_aln.fa" = This includes only those individuals that fulfill the (missing_data_per_row_threshold) threshold
	 *
	 * For now only the filtered dataset is used because IQtree runs into issues with too much missing data
	 * 1) For process (window_msa_to_phy_SINGLE) -- Infer each tree as a separate SLURM job
	 * 2) For process (window_msa_to_phy_COMBINED) -- Infer all trees for a subset/chromosomes as a single SLURM job
	 */

	label 'SC_SM_HT'
    
	tag "Create window alignments for each chromo"

	input:
	tuple val(chromo), \
	val(p1), \
	val(p2), \
	val(p3), \
	val(out), \
	path(chromo_msa) from combo_msa_by_chromo_raw_ch1
	each subset from params.window_sizes

	output:
	tuple val(chromo), \
	val(p1), \
	val(p2), \
	val(p3), \
	val(out), \
	val(subset), \
	file("*filt_indivs_aln.fa") optional true into combo_indiv_filt_window_msa_by_subset_ch1, combo_indiv_filt_window_msa_by_subset_ch2
	tuple val(chromo), \
	val(p1), \
	val(p2), \
	val(p3), \
	val(out), \
	val(subset), \
	file("*all_indivs_aln.fa") optional true into combo_indiv_all_window_msa_by_subset_ch1

	script:
	"""

	01_window_subset_MSA.py \
	-m ${chromo_msa} \
	-s ${subset} \
	-j ${params.window_jump} \
	-c ${params.missing_data_per_col_threshold} \
	-w ${params.missing_cols_per_window_threshold} \
	-r ${params.missing_data_per_row_threshold} \
	-i ${params.missing_indivs_per_window_threshold} \
	-p ${chromo}

	"""
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 		TIME FOR TREE BUILDING!! 												   *
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 																				   *
 * 		- Maximum likelihood phylogeny for each window in each subset
 *
 * 																				   *
 ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 

 /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		MAXIMUM-LIKELIHOOD PHYLOGENY FOR EACH WINDOW 	    *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  *
 *
 * 		NOTE, the 01_window_subset_MSA.py script emits a 
 *		single list of MSA per chromo/subset and to 
 *		parallelize the process, a single item for EACH 
 *		MSA is needed. In example, with the following info:
 * 		[chromo, combo, subset, aln_path, basename of alignment]
 * 
 */


combo_indiv_filt_window_msa_by_subset_ch1
	.transpose()
	.map{ it ->

		def chromo = it[0]
		def p1 = it[1]
		def p2 = it[2]
		def p3 = it[3]
		def out = it[4]
		def combo = p1 + '_' + p2 + '_' + p3 + '_' + out
		def subset = it[5]
		def aln_path = it[6]
		def basename = file(aln_path).getSimpleName()

		[chromo, combo, subset, aln_path, basename]
	}
	.set { combo_indiv_filt_window_msa_by_subset_annotated_ch1 }


combo_indiv_filt_window_msa_by_subset_ch2
	.map{ it ->

		def chromo = it[0]
		def p1 = it[1]
		def p2 = it[2]
		def p3 = it[3]
		def out = it[4]
		def combo = p1 + '_' + p2 + '_' + p3 + '_' + out
		def subset = it[5]
		def aln_list = it[6]

		[chromo, combo, subset, aln_list]
	}
	.set { combo_indiv_filt_window_msa_by_subset_annotated_ch2 }


process window_msa_to_phy_SINGLE {

	/*
	 * Take each window MSA and infer phylogeny using IQtree2
	 *
	 * Directed into one channel:
	 * - ch1: process (calc_topo_freqs) -- Calculate the distribution of gene tree freqs
	 */
    

	tag "Create window phylogenies for all windows (across all subsets)"
	publishDir "${params.outputdir}/02.phylogenies/00.windows/${subset}/${combo}/${chromo}/", mode:'copy'

	input:
	tuple val(chromo), \
	val(combo), \
	val(subset), \
	path(msa_fn), \
	val(basename) from combo_indiv_filt_window_msa_by_subset_annotated_ch1

	output:
	tuple val(chromo), \
	val(combo), \
	val(subset), \
	file("${basename}.treefile") into combo_window_phy_by_subset_ch1

	when:
	params.window_msa_to_phy_in_parallel == true

	script:
	"""

	iqtree \
	-s ${msa_fn}  \
	-m MFP \
	-nt AUTO \
	-ntmax ${task.cpus} \
	--redo \
	--prefix ${basename}

	"""
}


process window_msa_to_phy_COMBINED {

	/*
	 * Take each window MSA and infer phylogeny using IQtree2
	 *
	 * Directed into one channel:
	 * - ch1: process (calc_topo_freqs) -- Calculate the distribution of gene tree freqs
	 */
    

	tag "Create window phylogenies for all windows (across all subsets)"
	publishDir "${params.outputdir}/02.phylogenies/00.windows/${subset}/${combo}/${chromo}/", mode:'copy'

	input:
	tuple val(chromo), \
	val(combo), \
	val(subset), \
	file(aln_list) from combo_indiv_filt_window_msa_by_subset_annotated_ch2

	output:
	tuple val(chromo), \
	val(combo), \
	val(subset), \
	file("*.treefile") into combo_window_phy_by_subset_ch2

	when:
	params.window_msa_to_phy_in_parallel == false

	script:
	"""
	for aln_fn in ${aln_list}
	do
		BASENAME=`basename \$aln_fn .fa`

		iqtree \
		-s \$aln_fn \
		-m MFP \
		-nt AUTO \
		-ntmax ${task.cpus} \
		--redo \
		--prefix \$BASENAME
	done

	"""
}



 /*
 *
 * 		NOTE, each tree will be individually or by chromosome added to the channel.
 *		To process them individually would be computationally inefficient.
 *		Moreover, they should be independently calculated by comparison.
 *		Therefore we need to process the channel accordingly:
 *		
 * 		[subset, combo, [list of trees]]
 *
 *		NOTE: Eventually we want to treat windows from sex chromosomes differently
 *		but we will assign chromosome origin in the tsv summary sheet based on filenames 
 */


/*
 * 		Now create an input channel where all trees are stored by subset and combo
 */

if (params.window_msa_to_phy_in_parallel == true) {
	combo_window_phy_by_subset_ch1
		.groupTuple(by: [2,1])
		.map{ it ->

			def chromo = it[0]
			def combo = it[1]
			def subset = it[2]
			def trees = it[3].flatten()

			[subset, combo, trees]
		}
		.set { phy_by_subset_combo_ch1 }
} else {
	combo_window_phy_by_subset_ch2
		.groupTuple(by: [2,1])
		.map{ it ->

			def chromo = it[0]
			def combo = it[1]
			def subset = it[2]
			def trees = it[3].flatten()

			[subset, combo, trees]
    		}
		.set { phy_by_subset_combo_ch1 }
}

/*
 * 		Summarise Gene Tree frequency distributions per combination using custom Python script
 */

process calc_topo_freqs {

	/*
	 * Take a list of gene trees and caluc
	 *
	 * Directed into one channel:
	 * - ch1: process (summary_across_comps) -- Summarise gene tree frequency distributions across all comparisons
	 */

	label 'ete3'

	tag "Summarise gene tree frequencies per combination"
	publishDir "${params.outputdir}/03.summaries/00.windows/${subset}/${combo}/", mode:'copy'

	input:
	tuple val(subset), \
	val(combo), \
	path(window_phy_fn_list) from phy_by_subset_combo_ch1
	val(sex_chromo_list) from sex_chromos

	output:
	tuple val(subset), \
	val(combo), \
	file("${combo}_gene_tree_freqs_WIDE.tsv") into gt_freqs_by_subset_combo_ch1
	file('*.png')

	script:
	"""

	02_calc_topo_freqs.py \
	-c ${combo} \
	-p ${window_phy_fn_list} \
	-i ${params.chromo_coords} \
	-b ${params.bin_size} \
	-s ${sex_chromo_list}

	"""
}

/*
 * 		Now create an input channel where all summary tsv files are stored by subset size
 */

gt_freqs_by_subset_combo_ch1
	.groupTuple(by: [0])
	.map{ it ->
	
		def subset = it[0]
		def combo = it[1]
		def tsv = it[2].flatten()

		[subset, tsv]
	}
	.set { gt_freqs_by_subset_ch1 }


/*
 * 		Summarise topology frequencies ACROSS all comparisons
 */

process summary_across_comps {

	/*
	 * Take a list of tsv summary files and summarise frequency distributions across all comparisons
	 */

	label 'ete3'
	label 'SC_SM_IT'

	tag "Summarise gene tree frequencies across all combinations"
	publishDir "${params.outputdir}/03.summaries/00.windows/${subset}/", mode:'copy'

	input:
	tuple val(subset), \
	path(summary_tsv_fn_list) from gt_freqs_by_subset_ch1
	val(sex_chromo_list) from sex_chromos

	output:
	file('*.png')

	script:
	"""

	03_topo_comp_summary.py \
	-t ${summary_tsv_fn_list} \
	-i ${params.chromo_coords} \
	-s ${sex_chromo_list}

	"""
}

