#!/bin/bash
clear

SKIP=1

show_help()
{
echo "
        Usage: fas2fa.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nRename Fasta .fas to .fa files.\nThis script scans all the child directories for fasta files and renames the files removing spaces and changes the extension to ".fa" (Required for most scripts).\n\n'
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
	ls *.fas | grep -P 'C.*|PHE.*' | grep -P 'modal' -v | sed 's/.fas//' | while read F; do
		echo Processing "${F}.fas"
		mv "${F}.fas" "${F}.fa"
	done
	cd ..
done
