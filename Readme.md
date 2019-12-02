# CUBES
<!--Version Nov. 15th, 2018.-->

**Codon Usage Bias Evolutionary Scripts** is a software package designed to study the evolutionary traits of codon bias. To that end, a set of progressively ancestral core-genomes must first be obtained using for example [Edgar](https://edgar.computational.bio.uni-giessen.de) or [Get-Homologues](https://github.com/eead-csic-compbio/get_homologues) software. The core-genome sets (e.g. C1 -> Cn) are constructed by the successive incorporation of more distantly related genomes -following the phylogeny of the family/genus- to the analysis. Next, these Cn gene sets are used to calculate the modal codon use frequencies (for each set), the relative synonymous codon use (RSCU, for each gene), and to perform a Correspondence analysis. Additionally, tools to analyse the evolutionary traits of codon bias are provided.

CUBES is **pipeline** constituted by a set of bash and perl scripts which can be used to:

1. calculate modal codon usage frequencies for a set of genes (e.g. from a core genome obtained from [Edgar software](https://edgar.computational.bio.uni-giessen.de)) using [G. Olsen Software](http://www.life.illinois.edu/gary/programs/codon_usage.html)
2. generate representative DNA sequences for the calcualted modal frequencies (Since the codon usage analysis will be done using codonW which requires coding sequences)
3. calculate the correspondence analysis of RSCU using [CodonW software](http://codonw.sourceforge.net/). At this point, the pipeline generates plots of the first 2 components of COA, both for genes and codons, including the modal sequences for the core-genomes, singletons and putatively highly expressed proteins. Additionally, scripts to calculate CA of RSCU and raw codon counts using only R software are provided. 
4. calculate the adaptation index s-tAI and the GC3 content, using RGF (i.e. relative content of tRNA genes for each amino-acid) and the evolutionary distance between core-genomes (which can be obtained from Edgar software), which must be supplied by the user. As a result, the script outputs tables and plots showing the variation of these indices in function of the evolutionary distance. 
5. generate a distance tree based on the tRNA gene content (which must be supplied, e.g. from [GtRNAdb](http://gtrnadb.ucsc.edu/) data base) as in [Novoa *et al.* 2012](https://doi.org/10.1016/j.cell.2012.01.050) 
6. Finally, two plots are created. One showing the change in the codon use frequencies (CUF) for each codon and the corresponding *w*, and the other, a heatmap of the difference in CUF between: a. the initial and the most ancestral core-genomes, b. the initial core-genome and the putatively highly expressed genes (PHE).
----------

## Installation

For the correct function of the pipeline all the required software (see below) must be installed following the author instructions. Additionally, all the programs must be included in the linux $PATH.   

### Requirements
- The software needs G. Olsen software which can be found in [link](http://www.life.illinois.edu/gary/programs/codon_usage.html). The installation instructions are clearly provided by the author. The software must be in the linux $PATH to work correctly. Some of the G. Olsen libraries required minor modifications to run without error on ubuntu. (gjocodonlib.pm: on line 1495 ncpu was replaced by nproc; additionally, on freqs_2_nj_tree_linux, line 119 was replaced by system( 'echo "y" | phylip neighbor' );). Modified files are provided with the software.
- a local installation of [codonw](http://codonw.sourceforge.net/).
- [phylip](http://evolution.genetics.washington.edu/phylip.html) package for the tRNA tree generation (can be installed from linux repository. ~$sudo apt install phylip).
- [R-project] (https://www.r-project.org/) software with the following libraries: stringr, ggplot2, ggrepel, ggthemes, ggpmisc, data.table, tidyr, ggpubr, pheatmap, ape, devtools, tAI, doParallel, ade4.
- For tAI package installation follow author instructions on [tAI](https://github.com/mariodosreis/tai) github page.
- For the installation of the required R packages, some dependencies may be required.
  - From the linux repository install (for ubuntu):
    `sudo apt install libudunits2-dev libgdal-dev libssl-dev libmagick++-dev`
  - Additional R packages that might be required: httr, gh, usethis, magick, animation
- Svg_stack python script should be also in the $PATH. Download from [svg_stack](https://github.com/astraw/svg_stack)
- For CA analysis analysis using R, additional packages are required: gtools, Factoextra and FactoMiner.

## Running the pipeline - General comments
The scripts are programmed to scan for a set of files for each bacteria. These files must be located on a folder within the installation path. An example of the correct folder architecture is shown below: 
  
```  
────Installation-folder/─┬─Family1-files/  
                         ├─Family2-files/ 
                         ╎
                         ╎
                         └─ ┅
```
                      
The folder for each bacterial family must contain: 
- multifasta files named C(n).fa (for Core genomes), single*.fa (for singletons), GNM.fa (complete genome) each containing the CDS of a coregenome/singletons constructed by sequentially adding a more distant genome as in Lopez JL. *et al.* The sets of C genes for the analysis can be obtained using tools such as [EDGAR](https://edgar.computational.bio.uni-giessen.de/cgi-bin/edgar_login.cgi) or [Get_Homologues/Get_Phylomarkers](https://github.com/vinuesa).
- PHE.fa, putatively highly expressed genes (as defined by Lopez JL. *et al.*)
- dist.txt, which contain the distance between the species, extracted from the species phylogenetic tree. 
- rgf.txt (only needed for I, index calculation), is a file containing a normalized frequency based on the number of tRNAs for each codon/amino-acid. The order of the *codons* must follow the example file provided.
- trna.txt, the number of copies of tRNA. Can be obtained from [GtRNAdb](http://gtrnadb.ucsc.edu/). The order of the *codons* must follow the example file provided.

## Running the scripts in the correct order
The pipeline contains scripts which generate files that can be required on following steps. It is therefore recommended to run the scripts in order.
For the principal pipeline run `./calculate_ALL.sh -y` 
 
### 1. Correspondence analysis (CA) and Weighted CA (WCA)
CA (Correspondence analysis) of RSCU is calculated using codonw software.
To run the first analysis use the following scripts:
  
`./calculate_modals2.sh` -> calculates modal codon use frequencies using G.Olsen software, and generates a representative DNA sequence  (uses seq2AvgAA.pl and modal2seq2.pl perl scripts).
  
`./coa_GNM.sh` -> Calculates CA on RSCU, using CodonW software, and generates several plots in R. Requires single*.fa file, containing singleton gene CDS, and GNM.fa (complete genome CDS). Additionally, CA Axis are reoriented in the C1-Cn direction.  
  
The following scripts require trna.txt (see below) and modal sequences calculated with calculate_modals2.sh.  

`coa_CCounts_GNM.sh`->  Calculates CA on raw codon counts using G. Olsen codon count script and R software. Generates genes and codon plots. On the codons plot, only the codons which present the most significant increase in usage from C1 to PHE are labeled and colored.  
`coa_GNM_R.sh` -> Calculates CA on RSCU using seq2rscu.pl script and R software. Generates genes and codon plots. Codons are labeled and ploted as with coa_CCounts_GNM.sh.  
`coa_GNM_rv.sh`-> Similar to coa_GNM.sh but with some corrections and plot improvements. Codons are labeled and ploted as with coa_CCounts_GNM.sh.  

#### Additional scripts
`./calculate_modals.sh` -> calculates modal codon use frequencies using G.Olsen software, and generates a representative DNA sequence  (uses seq2AvgAA.pl and modal2seq2.pl perl scripts). This script doesn't take account of amino-acid composition.
  
`./coa.sh` -> Calculates CA on RSCU, using CodonW software, and generates several plots in R showing the evolutionary traits of codon bias.
  
`./coa_ws.sh` -> Calculates CA on RSCU, using CodonW software, and generates several plots in R. Requires single*.fa file containing singleton gene CDS.
  
`./wca.sh` -> Similar to coa.sh, but calculates Weighted CA (WCA) on RSCU, using R.
  
`./coa_GNM_PNG.sh` -> Calculates CA on RSCU, using CodonW software, and generates several plots in R. Requires single*.fa file, containing singleton gene CDS, and GNM.fa (complete genome CDS). Outputs a png file instead of svg vectorial graphic.

### 2. Calculate Index vs Evolutionary distance plots
`./Index.sh` -> Calculates I index (described on Lopez JL. *et al.*) for all the genes and generates an average I*C1->n* vs evolutionary distances plot.
  
`./Index_Modal.sh` -> Calculates I index from modal codon use frequencies and generates a I*C1->n* vs evolutionary distances plot.
  
`./GC3.sh` -> Calculates GC3 for all the genes and generates an average GC3*C1->n* vs evolutionary distances plot.
  
`./GC3_Modal.sh` -> Calculates GC3 index using modal codon use frequencies and generates an average GC3*C1->n* vs evolutionary distances plot.

#### Summary images 1
`./make_summary_images.sh` -> Makes summary images of the previous plots for rapid visual inspection of the results.

### tRNA Adaptation Index (tAI) calculation
A second step in the analysis is the calculation of the tRNA Adaptation Index. This index, defined by dos Reis M. (2004), and improved by Sabi *et al.* (2014).
   
>The tRNA adaptation index (tAI) is a measure of the level of co-adaptation between the set of tRNA genes and the codon usage bias of protein-coding genes in a given genome. STOP and methionine codons are ignored. The standard genetic code is assumed. dos Reis. M. 2004.
   
For the calculation of tAI, ws (relative adaptiveness values) are required. ws calculation depends on the parameter *Sij* which represents the efficiency of Watson and Crick / wobble base interaction. *Sij* values were estimated by dos Reis, as those values that maximize the correlation between tAI and protein abundance. A Nelder-Mead algorithm was used.
Sabi *et al.* 2014 calculated *Sij* by maximizing, using a Hill-climbing algorithm, the correlation between tAI and DCBS (Directional codon bias score, Sabi *et al.* 2014), and a software package (matlab based) was developed to estimate *Sij* for any organism. This approach, allows the calculation of species specific *Sij* without the need of proteomic data.
In this work, we used a similar approach to estimate *Sij* using R-project software. The main difference is that our software uses the optim library and search for the *Sij* that maximizes the tAI vs DCBS correlation by a Nelder-Mead algorithm. Another difference is that dos Reis tAI R package eliminates MET condons prior to normalization by max(W). Tuller's implementation doesn't seem to do that.
Additional scripts to calculate tAI using Sabi *et al.* 2014 average bacterial *Sij* is provided.   
  
The following scripts require: 
-tRNAs files (trna.txt, with two columns-> codon trna separated by "space character") with the number of copies of all the tRNAs present in the bacteria. This file must include all the 64 codons, STOPs too, as specified by dos Reis *et al.* 2004.
-C1*.fa file containing CDS. 
  
#### Estimation of species *Sij*.
`./calculate_sopt_DCBS_GNM_f.sh -y -i niter` -> Estimates *Sij* from the correlation between tAI and DCBS using the complete genome CDS. niter: number of random start points. 
  
`./tAi_Modal_g.sh` -> Calculates tAI using Sopt-DCBS-GNM for modal sequences from Ci->Cj, and generates a plot of tAI vs evolutionary distance.
  
`./tAi_Modal_ws_St_Sdr.sh`-> Calculates tAI for modal sequences from single->Ci->Cj using the average *Sij* from Sabi and Tuller 2014, and generates a plot of tAI vs evolutionary distance.

##### Other optional analysis
`./calculate_sopt_DCBS_GNM2_f.sh -y -i niter` -> Estimates *Sij* from the correlation between tAI and DCBS using the complete genome CDS. niter: number of random start points. Incorporates to tAI calculation the U:U wobble base interaction.
  
`./tAi_Modal_g2.sh` -> Calculates tAI using Sopt-DCBS-GNM for modal sequences from Ci->Cj, and generates a plot of tAI vs evolutionary distance. Incorporates to tAI calculation the U:U wobble base interaction.

###### Some more...
Estimation of *Sij* from tAI vs Nc-adjusted correlation
`./calculate_sopt_folder.sh -y -i 100` -> Estimates *Sij* from the correlation between tAI and Nc-adjusted (DosReis).
  
Estimation of *Sij* from tAI vs DCBS using the C1 CDS
`./calculate_sopt_DCBS_f.sh -y -i 100` -> Estimates *Sij* from the correlation between tAI and DCBS using the Core Genome (C1).
  
`./calculate_sopt_DCBS_f_ws.sh -y -i 100` -> Estimates *Sij* from the correlation between tAI and DCBS using C1 and singeltons.
  
`./tAi_Modal_ws.sh` -> Calculates tAI using Sopt-DCBS from the previous analysis.
  
`./calculate_sopt_DCBS_f_ws2.sh -y -i 100` -> Estimates *Sij* from the correlation between tAI and DCBS using C1 and singeltons, and incorporates de U:U wobble base interaction in the tAI calculation.
  
`./tAi_Modal_Stuller_Sdosreis.sh` -> Calculates tAI for modal sequences from Ci->Cj using the average *Sij* from Sabi and Tuller 2014, and generates a plot of tAI vs evolutionary distance.

##### tAI Summary images
`./make_tAi_images.sh` -> Makes tAI summary images for rapid visual inspection
`./make_tAi_images2.sh` -> Makes tAI (U:U wobble base interaction) summary images for rapid visual inspection

### tRNA copy number distance tree
This script uses the tRNA copy number and distribution to calculate a distances tree, as in Novoa *et al.* 2012. Requires the Phylip package for NJ tree calculation.
Requires:
- a raw txt file containing the number of tRNA-copies for each codon and for all species (Required columns tab-separated: SP	Name GCA	GCG ...).
- a rename.txt file containing the name as figures on the tree and the desired name on the leaves separated by one space (e.g. "A.gere" "Actinomyces_gerencseriae_DSM_6844_NZ_AUBN01000031")
  
`calculate_tRNA_dist.sh tRNA_copy_file` -> Calculates a distance tree based on Novoa *et al.*

### Heatmaps of codon usage from Ci->Cj->PHE and plots of frequencies of codon use vs evolutionary distance.
The following scripts generate a series of plots for the observation of codon bias changes with evolutionary distance, for each codon.
  
`./Cdelta_Heatmap2.sh` -> Generates a heatmap of relative adaptiveness values (ws), and of the difference of CUF between C1 and Cn, and between Cn and PHE (ΔCn-C1 and ΔPHE-Cn respectively). Allows the easy visualization of CUF changes from C1 to Cn and if there is further adaptation from Cn to PHE.
  
`./Cdelta_Heatmap2t.sh` -> The same as before, but with *Sij* taken from Sabi *et al.* 2014.
  
`./Cdelta_Heatmap2_U.sh` -> The same as the first, but tAi values are calculated taking in account U:U wobble base interaction.	
  
`./Cdelta_plot.sh` -> Generates a plot showing the evolution of the frequency of synonymous codon use from C1 to Cn, for every codon. Bars represents the tRNA copy number. Between brackets is the ws (calculated without the normalization by the maximum value)
  
`./Cdelta_plot_F3.sh` -> Heatmap of relative adaptiveness values (ws) and Delta C1->n. Wiht PHE, Singletons and color codes. Use to generate a figure like the one on Lopez et al. (not publish yet...)
  
`./Cdelta_plot_ws.sh` -> The same as Cdelta_plot.sh, but includes as first point the singletons.
   
`./Cdelta_plot_U.sh`    -> The same as before, but tAi values are calculated taking in account U:U wobble base interaction.
  
`./Cdelta_plot_wt.sh` -> The same as Cdelta_plot.sh, but for *Sij* taken from Sabi *et al.* 2014.

#### Additional scripts
`./Cdelta_Heatmap.sh` -> Generates a heatmap of Delta C1->n, PHE-RSCU, and tRNA copy number for each species. 

## Accessory tools and files

### Bash scripts
`./calculate_fc_all_folder.sh` -> Calculates codon counts for all .fa files in different folders.
`./calculate_fcounts_folder.sh` -> Calculates codon counts for all C1*.fa files in different folders.
`./calculate_gfcounts_folder.sh` -> Calculates codon counts for all GNM.fa files in different folders.
`./calculate_fmcounts_folder.sh` -> Calculates codon counts for modal sequences in different folders.
`./calculate_DCBS_f.sh` -> Calculates the Directional Codon Bias for all .fa files in different folders.

### Perl scripts
`./modal2seq2.pl` -> Calculates a coding sequence from the modal codon use calculated by G.Olsen software. Requires amino-acid composition file.
`./modal2seq.pl` -> Calculates a coding sequence from the modal codon use calculated by G.Olsen software. Doesn't take in account the amino-acid composition of the sequence.
`./seq2AvgAA.pl` -> Calculates the Amino-acid frequencies from sequences.
`./seq2DCBS.pl` -> Calculates the Directional Codon Bias from sequences.
`./seq2rscu.pl` -> Calculates the Relative Synonymous Codon Usage from sequences.
`./fas2fa.sh` -> changes the fas to fa extension, required by the scripts.

### headers, files used to order the data, and files used in perl scripts
G.Olsen -> freq_head.txt / chead.txt  
For codon order -> codOrder.txt, codons.txt, codons2.txt

### _ws variables
If the singletons for the reference genome are available, using _ws variants of scripts, it is posible to incorporate the singletons in the coa and tAI analysis. requires a "singletons*.fa" archive. 

### _t // _dr variables
The _t / -dr is used when the Sij were taken from Sabi 2014 o DosReis papers.


