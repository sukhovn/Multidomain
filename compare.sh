#!/bin/bash

COMPPATH=/home/nick/Dropbox/Study/Princeton/NumEq/Multidomain/Old/test_results

for file in $(cd test_results/; find *)
do
	if [ -f $COMPPATH/$file ]
	then
		# echo "Attempting to diff $file"
		if ! [ -z "$(diff $COMPPATH/$file test_results/$file)" ]
		then
			echo "Diffing $file"
			diff $COMPPATH/$file test_results/$file
			echo $'\n\n'
		fi
	fi
done