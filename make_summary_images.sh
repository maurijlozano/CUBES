#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: make_summary_images.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nMake_summary_images.\n\nThis script Compiles the most important figures for visual inspection.\nGenerates a single .svg file containting the graphics.\n\n'
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

    if [[ ! -f "./${F}/"CA_genes_GNM.svg || ! -f "./${F}/"CA_codons_GNM.svg ]] ; then
        echo 'You need to run coa.sh for CA analysis first.'
        if [[ ! -f "./${F}/"IvsD_modal.svg || ! -f "./${F}/"GC3vsD_modal.svg ]] ; then
               echo 'You need to run Index_Modal.sh | GC3_modal.sh first.'
               continue
        fi
        svg_stack.py --direction=h --margin=30 "./${F}/"IvsD_modal.svg "./${F}/"GC3vsD_modal.svg  > "${FN}".svg
        continue
    else
        if [[ ! -f "./${F}/"IvsD_modal.svg || ! -f "./${F}/"GC3vsD_modal.svg ]] ; then
               echo 'You need to run Index_modal.sh | GC3_modal.sh first.'
        svg_stack.py --direction=h --margin=30 "./${F}/"CA_genes_GNM.svg "./${F}/"CA_codons_GNM.svg  > "${FN}".svg
               continue
        fi
    fi

    svg_stack.py --direction=h --margin=30 "./${F}/"CA_genes_GNM.svg "./${F}/"CA_codons_GNM.svg  > 1.svg
    svg_stack.py --direction=h --margin=30 "./${F}/"IvsD_modal.svg "./${F}/"GC3vsD_modal.svg  > 2.svg
    svg_stack.py --direction=v --margin=50 1.svg 2.svg  > "${FN}".svg
    
    if [[ -f "./${F}/"CA_genes2_GNM_woPHE.svg && -f "./${F}/"IvsD_modal_woPHE.svg ]] ; then
    svg_stack.py --direction=h --margin=30 "./${F}/"CA_genes_GNM.svg "./${F}/"CA_genes2_GNM_woPHE.svg "./${F}/"CA_codons_GNM.svg  > 3.svg
    svg_stack.py --direction=h --margin=30 "./${F}/"IvsD_modal.svg "./${F}/"IvsD_modal_woPHE.svg  > 4.svg
    svg_stack.py --direction=h --margin=30 "./${F}/"GC3vsD_modal.svg "./${F}/"GC3vsD_modal_woPHE.svg > 5.svg
    svg_stack.py --direction=v --margin=50 3.svg 4.svg 5.svg > "${FN}"2.svg
	rm 3.svg 4.svg 5.svg
    fi
    rm 1.svg 2.svg 

done
