#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: CCdelta_heatmap2.sh [-h]

        -h Show this Help
"
echo -e '\nCalculates FUC_CCmax-FUC_CCmin and generates a heatmap including Ws (calculated without normalization by the maximum value.\n\n'
exit 1
}

while getopts h option; do
    case "${option}" in
       h) show_help;;
    esac
done




ls -d */ | while read F; do 
    echo Processing "${F}"
	FE=$(echo "${F}" | sed -e 's/[_ ][_ ][_ ]/_/' | sed -e 's/[_ ][_ ]/_/')
    FS=$(echo "${FE}" | sed -e 's/\(^[A-Za-z0-9]*\)[_ ]\(.\)[A-Za-z0-9]*[_ ]\([A-Za-z0-9]*\).*/\1_\2.\3/')
    FN=$(echo "${FS}" | sed -e 's/\///')

    cd "./${F}/"
    
    if [[ ! -f fmdata.txt ]] ; then
		echo 'Calculating modal frequencies.'

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
    Rscript --vanilla heatmap2.r "./${F}" "${FN}"

done
