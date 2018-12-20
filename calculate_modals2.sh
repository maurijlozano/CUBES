#!/bin/bash
clear

SKIP=1

show_help()
{
echo "
        Usage: calculate_modals2.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"
echo -e '\nCompute modals frequencies from Fasta files.\nThis script scans all the child directories for fasta files and computes the modal frequencies using G. Olsen software.\nThen it generates a DNA sequence with modal codon usage, taking in account the aminoacidic frequencies.\n\n'

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

ls ./*/*.fa | sed 's/.fa$//' | while read F; do 
echo Processing "${F}.fa"
sed -e 's/-//g' "${F}.fa" > "${F}.faed"
mv "${F}.faed" "${F}.fa"
compile_codon_counts  < "${F}.fa" | modal_usage_from_counts -a $(nproc) > "${F}.modal_freq"
perl seq2AvgAA.pl "${F}.fa"
perl modal2seq2.pl "${F}.modal_freq" "${F}.aaav"
done
