#!/bin/bash

clear


SKIP=1
WOPHE=1

show_help()
{
echo "
        Usage: coa_GCCounts_GNM.sh [-y] [-h] [-w]

        -y Skips command confirmation
        -w Make calcualtions of Core-genes with PHE substracted
        -h Show this Help
"

echo -e '\nCompute Corresppondence analysis of Codon counts.\nThis script scans all the child directories for CC1 and singleton fasta files and computes the COA using R software.\n'

exit 1
}

while getopts ywh option; do
    case "${option}" in
        y) SKIP=0;;
        w) WOPHE=0;;
        h) show_help;;
    esac
done


if [ ${SKIP} != 1 ]
then
    echo -e 'Beginning calculations.\n'
else
    echo -e 'PRESS ANY KEY TO CONTINUE.\n'
    read -n 1 -s
fi

ls -d */ | while read D; do 
    echo Processing "${D}"
    cd "./${D}/"
    compile_codon_counts < GNM.fa > GNM.counts
    sed -e 's/  / /g' GNM.counts > freq.counts
    cut -f 1 freq.counts > freq.count
    cut -f 2 freq.counts > freq.rname
    paste -d' ' freq.rname freq.count > freq.counts
    rm freq.count freq.rname
	mv freq.counts GNM.counts
	#Get ID of genes in each core set #
	ls *.fa* | grep 'modal.freq.fa' -v | grep 'woPHE' -v | grep '[Gg][Nn][Mm]' -v | while read F; do
		SET=$(echo ${F} | sed -e 's/.fa$//')
		sed -e 's/___*/|/' -e 's/   */|/' -e 's/.*locus_tag=\([^]]*\)\].*/>\1/' ${F} > ${F}ed
		grep '>' ${F}ed | cut -d'|' -f1 | sed -e 's/^>//' > ${SET}.genes
		echo ${SET} | sed -e 's/_[RrNn][eE][FfwW]//' > ${SET}.genes2
		cat ${SET}.genes2 ${SET}.genes > ${SET}.genes1
		mv ${SET}.genes1 ${SET}.genes
		rm ${F}ed ${SET}.genes2 
	done
	#SETS=$(ls *.fa | grep 'woPHE' -v | grep '[Gg][Nn][Mm]' -v | grep '[Pp][Hh][Ee]' -v | grep '[Ss][iI][Nn][gG][lL]' -v | wc -l)
	paste -d'\t' $(ls *.genes) > coreGenes.txt
	rm *.genes

	if [[ ! -f fmdata.txt ]] ; then
		echo 'Calculating modal frequencies counts.'
		cat $(ls *.modal_freq | sort -V) | sed -e 's/,/ /g' -e 's/\t/ /g' -e 's/|/ /g' > CC_modal_freqs.txt
		ls *.modal_freq | sort -V | sed -e 's/_\{0,2\}\(woPHE\).*/:\1/' -e 's/_.*//' -e 's/\..*//' -e 's/:/_/' > CCindex.txt
		paste -d' ' CCindex.txt CC_modal_freqs.txt > fmdata.txt
		rm CC_modal_freqs.txt CCindex.txt
	fi
		
	if [[ ! -f trna.txt ]] ; then
		echo 'File "trna.txt" required.'
		cd ..
		continue
	fi

    cd ..
    Rscript --vanilla graph_coaCCounts_GNM.R "./${D}" ${WOPHE}
done
