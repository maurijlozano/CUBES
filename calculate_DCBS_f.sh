#!/bin/bash
clear

SKIP=1
Ite=25
show_help()
{
echo "
        Usage: calculate_DCBS.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nCompute counts from Fasta files.\nThen the script scans the directorie for fasta files and computes DCBS.\n\n'
exit 1
}

while getopts "yh" option; do
    case "${option}" in
        y) SKIP=0;;
        h) show_help
            exit 1
            ;;
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
    echo Processing "${D}"
    cd "./${D}/"
    ls *.fa | while read F; do
    Fname=$(echo "${F}" | sed 's/.fa.\?$//')
    echo Processing "${F} fasta file"
    cd ..
    perl seq2DCBS.pl "${D}${F}"
    cd "./${D}/"
    done
    cd ..
done
