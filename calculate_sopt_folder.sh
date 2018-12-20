#!/bin/bash
clear

SKIP=1
Ite=25
show_help()
{
echo "
        Usage: calculate_sopt_folder.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
        -i Iterations of the optimization algorithm (25 default)
"
echo -e '\nCompute counts from Fasta files.\nThis script scans the directorie for fasta files and computes counts using G. Olsen software.\nThen, using codonW a table of Effective Number of Codons (ENc), G+C content of gene (all 3 codon positions) and GC of synonymous codons 3rd positions is generated.\n\n'

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

#generates the Highly expressed protein gc3, Nc etc table and codon counts.

ls -d */ | while read D; do
    echo Processing "${D}"
    cd "./${D}/"
    ls C1*.fa | grep 'woPHE' -v | grep '[Cc]1[0-9]' -v | while read F; do
    Fname=$(echo "${F}" | sed 's/.fa.\?$//')
    echo Processing "${F} fasta file"
    compile_codon_counts  < "${F}" > "${Fname}.count"
    codonw "${F}" "${Fname}.tab" "${Fname}.blk" -enc -gc -gc3s -L_aa -nomenu -silent -nowarn
    sed -e 's/  / /g' "${Fname}.count" > freq.counts
    sed -e 's/\r//' -e 's/ //g' -e 's/\t$//g' "${Fname}.tab" > "C1.table" #| grep -v '*'
    cut -f 1 freq.counts > freq.count
    cut -f 2 freq.counts > freq.rname
    paste -d' ' freq.rname freq.count > freq.counts
    rm *.count *.blk freq.rname *.tab
    done
cd ..

Rscript --vanilla tai_f.R "./${D}" "${Ite}"

done
