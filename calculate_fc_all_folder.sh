#!/bin/bash
clear

SKIP=1
show_help()
{
echo "
        Usage: calculate_fc_all_folder.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help

"

echo -e '\nCompute counts from Fasta files.\nThis script scans the directorie for fasta files and computes counts using G. Olsen software.\n\n'

exit 1
}

while getopts ":yh" option; do
    case "${option}" in
        y) SKIP=0;;
        h) show_help
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

#generates the codon counts. 

ls -d */ | while read D; do
    echo Processing "${D}"
    cd "./${D}/"
    ls *.fa | while read F; do
    Fname=$(echo "${F}" | sed 's/.fa.\?$//')
    echo Processing "${F} fasta file"
    compile_codon_counts  < "${F}" > "${Fname}.count"
    sed -e 's/  / /g' "${Fname}.count" > ${Fname}.counts
    cut -f 1 ${Fname}.counts > ${Fname}.count
    cut -f 2 ${Fname}.counts > ${Fname}.rname
    paste -d' ' ${Fname}.rname ${Fname}.count > ${Fname}.counts
    rm *.count ${Fname}.rname
    done
cd ..

done
