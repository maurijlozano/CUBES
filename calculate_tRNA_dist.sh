#!/bin/bash

#die function, pass the arguments of the function to stderr and exits
die () {
    echo >&2 "$@"
    exit 1
}

show_help()
{
echo "
        Usage: calculate_tRNA_dist.sh tRNA_copy_file [-h] 

        -h Show this Help
"
echo -e '\nCalcualtes pairwise distances based on tRNA gene content.\n

Distance(A,B) = 1 - { Sum_i [(2*min(ntRNA_ai,ntRNA_bi)/Sum(ntRNA_ai+ntRNA_bi))] / n} \n

Where A and B are vectors containing information of tRNA gene copies for each codon (i).\n

Next, calculates a distance matrix for all the species and uses the phylip beighbor package to calculate a Neighbour Joining distance tree.\n\n\n\n'

exit 1
}

while getopts h option; do
    case "${option}" in
         h) show_help;;
    esac
done

FILE=${1}

clear


[ "$#" -eq 1 ] || die "You have to provide the tRNA copies file"

if [ -f ${FILE} ]; then
    Rscript --vanilla calculate_tRNA_dist.R "./${1}"

    sed -e 's/"//g' row.mat > row.mated
    paste row.mated dist.mat > dist.mated
    cat len dist.mated > infile

    rm row.mat row.mated dist.mat dist.mated len
    if [ -f outfile ]; then
       rm outfile
    fi
    if [ -f outtree ]; then
       rm outtree
    fi
    
phylip neighbor << ANSWERS
y
ANSWERS

    Rscript --vanilla rename_tips.R "${1}"

else
    echo "The file: ${FILE} doesn't exist"
fi
