#!/bin/bash
clear

SKIP=1
Ite=25
show_help()
{
echo "
        Usage: calculate_sopt_DCBS_f_ws.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
        -i Iterations of the optimization algorithm (25 default)
"
echo -e '\nCompute counts from Fasta files, including singleton genes.\nThis script scans the directorie for fasta files and computes DCBS.\nThen, optimizes s for tai calculation.\n\n'

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
    echo Processing "${D}"
    cd "./${D}/"
    ls C1*.fa single*.fa | grep 'woPHE' -v | grep '[Cc]1[0-9]' -v | while read F; do
        Fname=$(echo "${F}" | sed 's/.fa.\?$//')
        echo Processing "${F} fasta file"
          
        compile_codon_counts  < "${F}" > "${Fname}.count"
        sed -e 's/  / /g' "${Fname}.count" > ${Fname}.counts
        cut -f 1 ${Fname}.counts > ${Fname}.count
        cut -f 2 ${Fname}.counts > ${Fname}.rname
        paste -d' ' ${Fname}.rname ${Fname}.count > ${Fname}.counts
        rm *.count ${Fname}.rname
        cd ..
        perl seq2DCBS.pl "${D}${F}"
        cd "./${D}/"
    done
    tail -n +2 $(ls single*.DCBS) > single.tmp
    cat $(ls C1*.counts single*.counts | grep 'woPHE' -v | grep '[Cc]1[0-9]' -v) > freq_ws.counts
    cat $(ls C1*.DCBS single.tmp | grep 'woPHE' -v | grep '10' -v ) > DCBS_ws.txt
    rm C1*.DCBS single*.DCBS single.tmp
    cd ..

    Rscript --vanilla tai_DCBS_f_ws.R "./${D}" "${Ite}"
done
