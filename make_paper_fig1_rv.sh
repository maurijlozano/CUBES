#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: make_paper_fig1.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"
echo -e '\nmake_paper_fig1.\n\nThis script Compiles the most important figures for visual inspection.\nGenerates a single .svg file containting the plots for figure1.\n\n'
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
    FE=$(echo "${F}" | sed -e 's/[_ ][_ ][_ ]/_/' | sed -e 's/[_ ][_ ]/_/')
    FS=$(echo "${FE}" | sed -e 's/\(^[A-Za-z0-9]*\)[_ ]\(.\)[A-Za-z0-9]*[_ ]\([A-Za-z0-9]*\).*/\1_\2.\3/')
    FN=$(echo "${FS}" | sed -e 's/\///')

    if [[ ! -f "./${F}/"CA_genes_GNM_rv.svg || ! -f "./${F}/"CA_codons_GNM_rv.svg ]] ; then
        echo 'You need to run coa.sh for CA analysis first.'
        continue
    fi

    svg_stack.py --direction=h --margin=30 "./${F}/"CA_genes_GNM_rv.svg "./${F}/"CA_codons_GNM_rv.svg  > "${FN}"_fig1_rv.svg
    inkscape -z -e "${FN}"_fig1_rv.png "${FN}"_fig1_rv.svg
done

