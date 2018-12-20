#!/bin/bash
clear

SKIP=1

show_help()
{
echo "
        Usage: modals2seq_f.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nCompute modal sequences from modal frequencies for all files in children folders.\nThis script scans all the child directories for fasta files and modal frequencies, and calculates modal sequences (Requieres precomputed modal_freqs).\n\n'

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
	if [[ ! -f "${F}.modal_freq" ]] ; then
		echo 'File ${F}.modal_freq is not there, skipping.'
		cd ..
		continue
	fi

	echo Processing "${F}.fa"
	perl seq2AvgAA.pl "${F}.fa"
	perl modal2seq2.pl "${F}.modal_freq" "${F}.aaav"
done
