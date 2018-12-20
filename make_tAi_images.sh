#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: make_tAi_images.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\make_tAi_images.\n\nThis script Compiles the most important figures for visual inspection.\nGenerates a single .svg file containting the plots involving tAi.\n\n'
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

    if [[ ! -f "./${F}/"Tai_vs_DCBS_GNM.svg || ! -f "./${F}/"tAi_modal_GNM.svg ]] ; then
        echo 'You need to run tAi analysis first.'
        continue
    fi
    
    if [[ ! -f "./${F}/"tAi_modal_woPHE_GNM.svg ]] ; then
                echo 'You need woPHE modal frequences.'
    svg_stack.py --direction=h --margin=30 "./${F}/"Tai_vs_DCBS_GNM.svg "./${F}/"tAi_modal_GNM.svg > "${FN}"3.svg
        continue
    fi

    svg_stack.py --direction=h --margin=30 "./${F}/"Tai_vs_DCBS_GNM.svg "./${F}/"tAi_modal_woPHE_GNM.svg > "${FN}"3.svg
done
