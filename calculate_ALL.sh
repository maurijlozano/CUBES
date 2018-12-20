#!/bin/bash
clear
SKIP=1

show_help()
{
echo "
        Usage: calculate_ALL.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"
exit 1
}

while getopts yh option; do
    case "${option}" in
        y) SKIP=0;;
        h) show_help;;
    esac
done

./calculate_modals2.sh -y
./coa_ws.sh -y
./make_paper_fig1.sh -y
./GC3_Modal.sh -y
./calculate_sopt_DCBS_GNM2_f.sh -y -i 100
./tAi_Modal_g2.sh -y
./make_paper_fig2.sh -y
./calculate_tRNA_dist.sh -y
./Cdelta_plot2.sh -y
./Cdelta_heatmap2_2.sh -y
./make_paper_fig4.sh -y

echo 'Done!!!'
