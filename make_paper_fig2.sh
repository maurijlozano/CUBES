#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: make_paper_fig2.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"
echo -e '\nmake_paper_fig2.\n\nThis script Compiles the most important figures for visual inspection.\nGenerates a single .svg file containting the plots for figure2.\n\n'
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

    if [[ ! -f "./${F}/"tAi_modal_GNM.svg || ! -f "./${F}/"GC3vsD_modal_ws.svg ]] ; then
        echo 'You need to run tAI analysis first.'
        continue
    fi

    svg_stack.py --direction=h --margin=30 "./${F}/"tAi_modal_GNM.svg "./${F}/"GC3vsD_modal_ws.svg > "${FN}"_fig2.svg
    
done

