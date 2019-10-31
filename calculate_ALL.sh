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
./coa_GNM.sh -y
./make_paper_fig1.sh -y
./GC3_Modal.sh -y
./calculate_sopt_DCBS_GNM_f.sh -y -i 100
./tAi_Modal_g.sh -y
./make_paper_fig2.sh -y
./Cdelta_plot_F3.sh -y
./Cdelta_heatmap2.sh -y
./make_paper_fig4.sh -y

echo 'Done!!!'
