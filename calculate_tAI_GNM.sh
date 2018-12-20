#!/bin/bash
clear
SKIP=1

show_help()
{
echo "
        Usage: calculate_TAI_GNM.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"
}

while getopts yh option; do
    case "${option}" in
        y) SKIP=0;;
        h) show_help;;
    esac
done

    ./calculate_sopt_DCBS_GNM_f.sh -y -i 100
    ./tAi_Modal_g.sh -y
	./Cdelta_plot.sh -y
	./Cdelta_Heatmap2.sh -y

	./calculate_sopt_DCBS_GNM2_f.sh -y -i 100
    ./tAi_Modal_g2.sh -y
	./Cdelta_plot2.sh -y
	./Cdelta_Heatmap2_2.sh -y

echo 'Done!!!'
