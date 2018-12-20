#!/bin/bash
clear

SKIP=1
Ite=25
show_help()
{
echo "
        Usage: calculate_sopt_DCBS_f.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
        -i Iterations of the optimization algorithm (25 default)
"

echo -e '\nCompute counts from Fasta files.\nThis script scans the directorie for fasta files and computes DCBS.\nThen, optimizes s for tai calculation.\n\n'

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


if [  ! -f "${D}/freq.counts" ]
then
./calculate_fcounts_folder.sh -y
fi

#generates DCBS table. 

ls -d */ | while read D; do
    echo Processing "${D}"
    cd "./${D}/"
    ls C1*.fa | grep 'woPHE' -v | grep '[Cc]1[0-9]' -v | while read F; do
    Fname=$(echo "${F}" | sed 's/.fa.\?$//')
    echo Processing "${F} fasta file"
    cd ..
    perl seq2DCBS.pl "${D}${F}"
    cd "./${D}/"
    mv $(echo ${F} | sed -e 's/.fa/.DCBS/') DCBS.txt
    done
    cd ..

    Rscript --vanilla tai_DCBS_f.R "./${D}" "${Ite}"
done
