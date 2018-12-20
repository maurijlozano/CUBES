#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: tAi_Modal_g.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nComputes tAi Adaptaion Index of Modal frequencies.\nThis script scans all the child directories for C1 and singleton .modal_freq files and computes de adaptation index tAi.\nGenerates a tAi vs Distance plot.\n\n'

exit 1
}

while getopts yh option; do
    case "${option}" in
        y) SKIP=0;;
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


ls -d */ | while read F; do 
    echo Processing "${F}"
    cd "./${F}/"
    
    if [[ ! -f s_opts_DCBS_GNM.txt ]] ; then
        echo 'File "s_opts_DCBS_GNM.txt" is not there, skipping.'
        cd ..
        continue
    fi

    if [[ ! -f modal_ws.counts ]] ; then
		cat $(ls *.modal_freq.fas | sort -V) > modal_ws.counts
		compile_codon_counts  < modal_ws.counts > modal_ws.count
		sed -e 's/  / /g' modal_ws.count > modal_ws.counts
		cut -f 1 modal_ws.counts > modal_ws.count
		echo -e "$(ls *.modal_freq.fas | sort -V)" | sed -e 's/.modal_freq.fas//' -e 's/[_ ][RrNn][Ee][WwFf]$//' -e 's/__/_/' > modal_ws.rname
		paste -d' ' modal_ws.rname modal_ws.count | sed -e 's/\^\n//' > modal_ws.counts
		rm *.count modal_ws.rname
	fi

    cd ..

    Rscript --vanilla Tai_modal_g.R "./${F}"

done
