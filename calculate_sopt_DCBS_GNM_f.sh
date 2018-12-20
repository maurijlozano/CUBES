#!/bin/bash
clear

SKIP=1
Ite=25
show_help()
{
echo "
        Usage: calculate_sopt_DCBS_GNM_f.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
        -i Iterations of the optimization algorithm (25 default)
"

echo -e '\nCompute counts from Fasta files.\nThis script scans the directorie for GNM.fa files and computes DCBS.\nThen, optimizes s for tai calculation.\n\n'

exit 1
}

while getopts ":yhi:" option; do
    case "${option}" in
        y) SKIP=0;;
        h) show_help
            exit 1
            ;;
        i) Ite=$OPTARG;;
        \?) echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done



if [[ $# -eq 0 ]] ; then
    echo 'Parameters required, run -h for help'
    exit 1
fi


if [ ${SKIP} != 1 ]
then
    echo -e 'Beginning calculations.\n'
else
    echo -e 'PRESS ANY KEY TO CONTINUE.\n'
    read -n 1 -s
fi


#generates DCBS table. 

ls -d */ | while read D; do

	if [  ! -f "${D}/genome.freq.counts" ]
	then
	./calculate_gfcounts_folder.sh -y
	fi

    echo Processing "${D}"
    perl seq2DCBS.pl "${D}GNM.fa"
    cd "./${D}/"
    mv $(echo "GNM.fa" | sed -e 's/.fa/.DCBS/') GNM_DCBS.txt
    cd ..
    Rscript --vanilla tai_DCBS_g_f.R "./${D}" "${Ite}"
done
