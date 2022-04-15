## Setup

In order to prepare data, first build the cython library

```
cd scripts
python JointTadsSetup.py build_ext --inplace
cd ..
```

Next, get the data

```
mkdir -p data
wget http://usevision.org/data/mm10/rnaHtseqCountsall_withcoordinates.0.txt -o data/mm10_rna.txt
wget http://usevision.org/data/mm10/ideasJointMay2021/S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state -o data/mm10_state.txt
wget http://usevision.org/data/mm10/ideasJointMay2021/S3V2_IDEAS_mm10_ccre2.cCRE.M.notall0.bed -o data/mm10_cre.txt
wget http://usevision.org/data/hg38/RNA/Oct2021/cntsFeb21v3.tab -o data/hg38_rna.txt
wget http://usevision.org/data/hg38/IDEASstates/ideasJointMay2021/S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state -o data/hg38_rna.txt
wget http://usevision.org/data/hg38/IDEASstates/ideasJointMay2021/S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.bed -o hg38_cre.txt

```

Finally, prepare the data

```
python scripts/mm10_rna_to_npy.py -r data/mm10_rna.txt -o data/mm10_rna_all.npy
python scripts/states_to_npy.py -s data/mm10_state.txt -o data/mm10_state.npy
python scripts/CREs_to_npy.py -s data/mm10_cre.txt -o data/mm10_cre.npy
python scripts/hg38_rna_to_npy.py -r data/hg38_rna.txt -o data/hg38_rna_all.npy
python scripts/states_to_npy.py -s data/hg38_state.txt -o data/hg38_state.npy
python scripts/CREs_to_npy.py -s data/hg38_cre.txt -o data/hg38_cre.npy
```

## Running

In order to run, most of the settings are fine with default values

```
python bin/JointCRE_reformatted.py -r data/mm10_rna_all.npy data/hg38_rna_all.npy \
  -s data/mm10_state.npy data/hg38_state.npy -c data/mm10_cre.npy data/hg38_cre.npy \
  -g mm10 hg38 -o results -v 3
```

In order to leave out on celltype, use the `-l mm10,CLP` for example. To run controls, use `--shuffle cre` or `shuffle tss`.
