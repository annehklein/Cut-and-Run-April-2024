![image](https://github.com/user-attachments/assets/1b5695aa-61df-483e-bec2-29325a2b62ac)![image](https://github.com/user-attachments/assets/8e5c8943-f760-4f98-a28f-d1efa8974e63)# Cut and run analysis

Analysis of Mawj's cut and run data.

These data came from patient hepatocytes, engrafted into FRG mice.  Data were collected at P0 (primary hepatocytes), P1 (after one round of engraftment), or P2 (after two rounds of engraftment).

Mawj also did a spike-in of E.coli and [SNAP-CUTANA](https://www.epicypher.com/content/documents/SNAP-CUTANA_K-MetStat_user_guide.pdf) for normalization. Epicypher provides a [shell script](https://www.epicypher.com/content/documents/SNAP-CUTANA_K-MetStat_Panel_ShellScript.sh) on the [SNAP-CUTANA product page](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) that quantifies the SNAP-CUTANA spike-in using `grep`.

There are several protocols describing this analysis or parts thereof:

 - [CUT and TAG analysis protocol](https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1)
 - [Another Cut and TAG analysis protocol](https://yezhengstat.github.io/CUTTag_tutorial/#I_Introduction)
 - [A description of how to normalize with a scaling factor in deepTools](https://github.com/deeptools/deepTools/issues/1073)
 - [peak calling with MACS2](https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/06_peak_calling_macs.html)

## Workflow

I will write a snakemake workflow to do this analysis.  The main steps are:
s
- QC (fastqc)
- Align to human reference, remove duplicates
- Align to E. Coli reference for normalization, remove duplicates
- Count SNAP-CUTANA spike-in reads
- Quantify human reference and spike-in reads and calculate normalization factors
- Convert to bedgraph 
- Peak calling with SEACR or MACS3


### References

For the human reference, I use [this version](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) [recommended by Heng Li](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use).  

For the E. coli reference, I use the version  (NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome).

I ignored any blacklisted regions, the blacklist was obtained [here](https://github.com/Boyle-Lab/Blacklist/tree/master/lists).

##QC 
reads can be found [here](https://csiroau-my.sharepoint.com/:f:/r/personal/kle065_csiro_au/Documents/CMRI%20projects/Mawj/Cut%20and%20Run%20April%202024/reads?csf=1&web=1&e=0Ju2F9)
fastqc ran just after sequencing can be found [here](https://csiroau-my.sharepoint.com/:f:/r/personal/kle065_csiro_au/Documents/CMRI%20projects/Mawj/Cut%20and%20Run%20April%202024/fastqc?csf=1&web=1&e=MHnzCU)

Filtering reads
```
OTC 3 P0 H3K27Ac
Input:                  	134252768 reads 		20137915200 bases.
Low quality discards:   	5974104 reads (4.45%) 	896115600 bases (4.45%)
Total Removed:          	5974104 reads (4.45%) 	896115600 bases (4.45%)
Result:                 	128278664 reads (95.55%) 	19241799600 bases (95.55%)
```
```
OTC 3 P0 H3K27me3
Input:                  	139364122 reads 		20904618300 bases.
Low quality discards:   	5051566 reads (3.62%) 	757734900 bases (3.62%)
Total Removed:          	5051566 reads (3.62%) 	757734900 bases (3.62%)
Result:                 	134312556 reads (96.38%) 	20146883400 bases (96.38%)
```
```
OTC #3 P0 H3K4me3
Input:                  	130662074 reads 		19599311100 bases.
Low quality discards:   	5326812 reads (4.08%) 	799021800 bases (4.08%)
Total Removed:          	5326812 reads (4.08%) 	799021800 bases (4.08%)
Result:                 	125335262 reads (95.92%) 	18800289300 bases (95.92%)
```
```
OTC #3 P0 IgG
```
```
OTC #3 P1 H3K27Ac
Input:                  	154017766 reads 		23102664900 bases.
Low quality discards:   	9005470 reads (5.85%) 	1350820500 bases (5.85%)
Total Removed:          	9005470 reads (5.85%) 	1350820500 bases (5.85%)
Result:                 	145012296 reads (94.15%) 	21751844400 bases (94.15%)
```
```
OTC #3 P1 H3K27me3
Input:                  	154017766 reads 		23102664900 bases.
Low quality discards:   	9005470 reads (5.85%) 	1350820500 bases (5.85%)
Total Removed:          	9005470 reads (5.85%) 	1350820500 bases (5.85%)
Result:                 	145012296 reads (94.15%) 	21751844400 bases (94.15%)
```
```
OTC #3 P1 H3K4me3
Input:                  	130662074 reads 		19599311100 bases.
Low quality discards:   	5326812 reads (4.08%) 	799021800 bases (4.08%)
Total Removed:          	5326812 reads (4.08%) 	799021800 bases (4.08%)
Result:                 	125335262 reads (95.92%) 	18800289300 bases (95.92%)
```
```
OTC #3 P1 IgG
Input:                  	125958962 reads 		18893844300 bases.
Low quality discards:   	6152010 reads (4.88%) 	922801500 bases (4.88%)
Total Removed:          	6152010 reads (4.88%) 	922801500 bases (4.88%)
Result:                 	119806952 reads (95.12%) 	17971042800 bases (95.12%)
```
```
OTC #3 P2 H3K27Ac
Input:                  	136895994 reads 		20534399100 bases.
Low quality discards:   	7823934 reads (5.72%) 	1173590100 bases (5.72%)
Total Removed:          	7823934 reads (5.72%) 	1173590100 bases (5.72%)
Result:                 	129072060 reads (94.28%) 	19360809000 bases (94.28%)
```
```
OTC #3 P2 H3K27me3
Input:                  	134223916 reads 		20133587400 bases.
Low quality discards:   	6221424 reads (4.64%) 	933213600 bases (4.64%)
Total Removed:          	6221424 reads (4.64%) 	933213600 bases (4.64%)
Result:                 	128002492 reads (95.36%) 	19200373800 bases (95.36%)
```
```
OTC #3 P2 H3K4me3
Input:                  	141659896 reads 		21248984400 bases.
Low quality discards:   	6363856 reads (4.49%) 	954578400 bases (4.49%)
Total Removed:          	6363856 reads (4.49%) 	954578400 bases (4.49%)
Result:                 	135296040 reads (95.51%) 	20294406000 bases (95.51%)
```
```
OTC #3 P2 IgG
Input:                  	24818876 reads 		3722831400 bases.
Low quality discards:   	1369012 reads (5.52%) 	205351800 bases (5.52%)
Total Removed:          	1369012 reads (5.52%) 	205351800 bases (5.52%)
Result:                 	23449864 reads (94.48%) 	3517479600 bases (94.48%)
```
```
FRG #458+FRG 166 H3K27Ac
Input:                  	139611468 reads 		20941720200 bases.
Low quality discards:   	5875070 reads (4.21%) 	881260500 bases (4.21%)
Total Removed:          	5875070 reads (4.21%) 	881260500 bases (4.21%)
Result:                 	133736398 reads (95.79%) 	20060459700 bases (95.79%)
```
```
FRG #458+FRG 166 H3K27me3
Input:                  	144875966 reads 		21731394900 bases.
Low quality discards:   	7367672 reads (5.09%) 	1105150800 bases (5.09%)
Total Removed:          	7367672 reads (5.09%) 	1105150800 bases (5.09%)
Result:                 	137508294 reads (94.91%) 	20626244100 bases (94.91%)
```
```
FRG #458+FRG 166 H3K4me3
```
```
FRG #458+FRG 166 IgG
Input:                  	133547278 reads 		20032091700 bases.
Low quality discards:   	5095822 reads (3.82%) 	764373300 bases (3.82%)
Total Removed:          	5095822 reads (3.82%) 	764373300 bases (3.82%)
Result:                 	128451456 reads (96.18%) 	19267718400 bases (96.18%)
```
```
FRG #474+FRG#182 H3K27Ac
Input:                  	127259888 reads 		19088983200 bases.
Low quality discards:   	8006450 reads (6.29%) 	1200967500 bases (6.29%)
Total Removed:          	8006450 reads (6.29%) 	1200967500 bases (6.29%)
Result:                 	119253438 reads (93.71%) 	17888015700 bases (93.71%)
```
```
FRG #474+FRG#182 H3K27me3
Input:                  	143823766 reads 		21573564900 bases.
Low quality discards:   	8469372 reads (5.89%) 	1270405800 bases (5.89%)
Total Removed:          	8469372 reads (5.89%) 	1270405800 bases (5.89%)
Result:                 	135354394 reads (94.11%) 	20303159100 bases (94.11%)
```
```
FRG #474+FRG#182 H3K4me3
Input:                  	123309158 reads 		18496373700 bases.
Low quality discards:   	4793962 reads (3.89%) 	719094300 bases (3.89%)
Total Removed:          	4793962 reads (3.89%) 	719094300 bases (3.89%)
Result:                 	118515196 reads (96.11%) 	17777279400 bases (96.11%)
```
```
FRG #474+FRG#182 
Input:                  	130432406 reads 		19564860900 bases.
Low quality discards:   	7973768 reads (6.11%) 	1196065200 bases (6.11%)
Total Removed:          	7973768 reads (6.11%) 	1196065200 bases (6.11%)
Result:                 	122458638 reads (93
```

SNAP-CUTANA outputs can be found [here](https://csiroau-my.sharepoint.com/:f:/r/personal/kle065_csiro_au/Documents/CMRI%20projects/Mawj/Cut%20and%20Run%20April%202024?csf=1&web=1&e=hyp0XR)

Macs2 takes into account the different number of reads between sample and control

Normalisation is made for bigwig files 

Normalsiation factor is calcualted this way
average_ecoli_count = sum(ecoli_read_counts.values()) / len(ecoli_read_counts)
norm_factors = average_ecoli_count /count in ecoli_read_counts

# Align against GRCh38 and E.coli

|**Sample ID**	|**Index (i7) barcode**	|**Overall alignement (%)**	|**Number of reads aligned to E.coli**	|**Normalisation factor** 	|** Number of reads aligning agaisnt X chromosome** |
|:-----------|:-----------|:-----------|:-----------|:-----------|:-----------
|OTC #3 P0 H3K27Ac	|ATCACG|	99.6|	145800	|18.6388|2699364|
|OTC #3 P0 H3K27me3|	CGATGT|	99.59|	133275|	20.3904|3913455|
|OTC #3 P0 H3K4me3|	TTAGGC	|99.32	|37710	|72.064|504265|
|OTC #3 P0 IgG	|TGACCA|	99.55|	191892	|14.1618|2937105|
|OTC #3 P1 H3K27Ac	|ACAGTG	|96.75	|1624459|	1.6729|1741714|
|OTC #3 P1 H3K27me3	| GCCAAT	|97.8	|297255	|9.1421|4100476|
|OTC #3 P1 H3K4me3	|CAGATC	|91.44|	5529812|	0.4914|2126629|
|OTC #3 P1 IgG|	ACTTGA|	94.94|	7980688|	0.3405|2007048|
|OTC #3 P2 H3K27Ac|	GATCAG	|94.68|	661148	|4.1103|1830181|
|OTC #3 P2 H3K27me3|	TAGCTT|	95.7|	132118	|20.569|3895655|
|OTC #3 P2 H3K4me3	|GGCTAC	|91.23|	4087586	|0.6648|2180552|
|OTC #3 P2 IgG|	CTTGTA|	93.1|	2559729|	1.0616|2692173|
|FRG #458+FRG 166 H3K27Ac|	AGTCAA	|91.21|	419196|	6.4827|1537297|
|FRG #458+FRG 166 H3K27me3|	AGTTCC|	89.4|	166818|	16.2904|4824740|
|FRG #458+FRG 166 H3K4me3	|ATGTCA|	83.33|	4048658	|0.6712|2321623|
|FRG #458+FRG 166 IgG|	CCGTCC|	84.17|	3905827	|0.6958|2183303|
|FRG #474+FRG#182 H3K27Ac|	GTCCGC|	92.01|	1701344	|1.5973|1353419|
|FRG #474+FRG#182 H3K27me3|	GTGAAA	|91.45|	692115|	3.9264|4159333|
|FRG #474+FRG#182 H3K4me3|	GTGGCC|	83.65|	12723918|	0.2136|2271391|
|FRG #474+FRG#182 IgG|	GTTTCG|	84.88|	7311358|	0.3717|2051329|

# Results
[MACS2](https://csiroau-my.sharepoint.com/:f:/r/personal/kle065_csiro_au/Documents/CMRI%20projects/Mawj/Cut%20and%20Run%20April%202024/results?csf=1&web=1&e=dsjYC5)

[ChromHMM](https://csiroau-my.sharepoint.com/:f:/r/personal/kle065_csiro_au/Documents/CMRI%20projects/Mawj/Cut%20and%20Run%20April%202024/ChromHMM?csf=1&web=1&e=Go7TbX)

[bigWig files](https://csiroau-my.sharepoint.com/:f:/r/personal/kle065_csiro_au/Documents/CMRI%20projects/Mawj/Cut%20and%20Run%20April%202024/results?csf=1&web=1&e=Pc6v8E)

[bigWig - IgG](https://csiroau-my.sharepoint.com/:f:/r/personal/kle065_csiro_au/Documents/CMRI%20projects/Mawj/Cut%20and%20Run%20April%202024/results?csf=1&web=1&e=Pc6v8E)




