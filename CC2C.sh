#!/bin/bash
clear

SKIP=1

show_help()
{
echo "
        Usage: CC2CandPHE2RTRP.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nChange names from CCn to Cn and PHE to RTRP.\n\n'
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
ls -d */ | while read D; do
    echo Processing "${D}"
    cd "./${D}/"
	ls *.fa *.modal_freq *.fas *.counts *.aaav | grep -P '^_*CC.*|_*PHE.*' | while read F; do
		echo Processing "${F}"
		FE=$(echo "${F}" | sed -e 's/CC/C/' |  sed -e 's/_\{0,2\}sPHE/_woPHE/' | sed -e 's/[_ ][_ ]/_/')
		if [ ${F} != ${FE} ]
		then
			mv "${F}" "${FE}"
		fi
	done
	ls *.modal_freq.fas | while read F; do
		echo Processing "${F}"
		sed -e 's/_[Rr][Ee][Ff]$//' -e 's/_[Nn][Ee][Ww]$//'  -e 's/_\{0,2\}sPHE/_woPHE/' -e 's/_\{0,2\}PHE/PHE/' -e 's/\/_\{0,2\}CC/\/C/' "${F}" > "${F}"ed
		mv "${F}"ed "${F}"
	done
	sed -e 's/^CC/C/' dist.txt > dist.txted
	mv dist.txted dist.txt
	cd ..
done
