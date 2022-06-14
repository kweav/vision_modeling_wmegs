#!/bin/bash

#usage: ./run_analysis.sh [nolo | lo] [nearest | standard] [mm10 | hg38]
#note: [mm10 | hg38] command line argument is only needed if $1 is lo

# ran these base analyses without leaving out genome/celltype combos
if [ $1 = "nolo" ]; then
	if [ $2 = "nearest" ]; then
		#nearest gene
		##shuffle tss
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_shuffle_tss \
		-v 3 --nearest-gene --shuffle tss

		#nearest gene
		##shuffle tss
		###tss only
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_shuffle_tss_tssonly \
		-v 3 --nearest-gene --shuffle tss --cre-dist 0

		#nearest gene
		##shuffle cre
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_shuffle_cre \
		-v 3 --nearest-gene --shuffle cre

		#nearest gene
		##shuffle cre
		##cre only
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_shuffle_cre_creonly \
		-v 3 --nearest-gene --shuffle cre --promoter-dist 0

		#nearest gene
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng \
		-v 3 --nearest-gene

		#nearest gene
		##tss only
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_tssonly \
		-v 3 --nearest-gene --cre-dist 0

		#nearest gene
		##cre only
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_creonly \
		-v 3 --nearest-gene --promoter-dist 0
	fi

	if [ $2 = "standard" ]; then
		#standard
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard \
		-v 3

		#standard
		##tss only
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_tssonly \
		-v 3 --cre-dist 0

		#standard
		##cre only
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_creonly \
		-v 3 --promoter-dist 0

		#standard
		##shuffle tss
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_shuffle_tss \
		-v 3 --shuffle tss

		#standard
		##shuffle tss
		##tss only
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_shuffle_tss_tssonly \
		-v 3 --shuffle tss --cre-dist 0

		#standard
		##shuffle cre
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_shuffle_cre \
		-v 3 --shuffle cre

		#standard
		##shuffle cre
		###cre only
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_shuffle_cre_creonly \
		-v 3 --shuffle cre --promoter-dist 0
	fi
fi

# ran these extended analyses leaving out genome/cell type combos for each base analysis
function call_it {
	if [ $3 = "nearest" ]; then
		#nearest gene
		##shuffle tss
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_shuffle_tss_lo_${1}_${2} \
		-v 3 --nearest-gene --shuffle tss -l ${1},${2}

		#nearest gene
		##shuffle tss
		###only use tss
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_shuffle_tss_tssonly_lo_${1}_${2} \
		-v 3 --nearest-gene --shuffle tss --cre-dist 0 -l ${1},${2}

		#nearest gene
		##shuffle cre
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_shuffle_cre_lo_${1}_${2} \
		-v 3 --nearest-gene --shuffle cre -l ${1},${2}

		#nearest gene
		##shuffle cre
		##only use cre
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_shuffle_cre_creonly_lo_${1}_${2} \
		-v 3 --nearest-gene --shuffle cre --promoter-dist 0 -l ${1},${2}

		#nearest gene
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_lo_${1}_${2} \
		-v 3 --nearest-gene -l ${1},${2}

		#nearest gene
		##only use tss
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_tssonly_lo_${1}_${2} \
		-v 3 --nearest-gene --cre-dist 0 -l ${1},${2}

		#nearest gene
		##only use cre
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/ng_creonly_lo_${1}_${2} \
		-v 3 --nearest-gene --promoter-dist 0 -l ${1},${2}
	fi
	if [ $3 = "standard" ]; then
		#standard refinement
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_lo_${1}_${2} \
		-v 3 -l ${1},${2}

		#standard refinement
		##only use tss
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_tssonly_lo_${1}_${2} \
		-v 3 --cre-dist 0 -l ${1},${2}

		#standard refinement
		##only use cre
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_creonly_lo_${1}_${2} \
		-v 3 --promoter-dist 0 -l ${1},${2}

		#standard refinement
		##shuffle tss
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_shuffle_tss_lo_${1}_${2} \
		-v 3 --shuffle tss -l ${1},${2}

		#standard refinement
		##shuffle tss
		###only use tss
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_shuffle_tss_tssonly_lo_${1}_${2} \
		-v 3 --shuffle tss --cre-dist 0 -l ${1},${2}

		#standard refinement
		##shuffle cre
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_shuffle_cre_lo_${1}_${2} \
		-v 3 --shuffle cre -l ${1},${2}

		#standard refinement
		##shuffle cre
		##only use cre
		date; time python JointCRE_reformatted.py -r ../data/mm10_rna_all.npy ../data/hg38_rna_all.npy \
		-s ../data/mm10_state.npy ../data/hg38_state.npy \
		-c ../data/mm10_cre.npy ../data/hg38_cre.npy \
		-g mm10 hg38 -o ../results/nearest_gene_run_06022022/standard_shuffle_cre_creonly_lo_${1}_${2} \
		-v 3 --shuffle cre --promoter-dist 0 -l ${1},${2}
	fi
}

if [ $1 = "lo" ]; then
	genome=$3
	if [ $genome = "mm10" ]; then
		for celltype in CFUE ERY iMK MON CMP G1E LSK NEU ER4 GMP MEP
		do
			call_it $genome $celltype $2
		done
	fi
	if [ $genome = "hg38" ]; then
		for celltype in B EOS HUDEP1 MK MPP CLP ERY HUDEP2 MONc NEU CMP GMP MEP MONp NK
		do
			call_it $genome $celltype $2
		done
	fi
fi
