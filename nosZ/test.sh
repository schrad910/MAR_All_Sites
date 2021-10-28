#!/bin/bash

FA='nosZ_db_accno.fa'


function process_line() {
	IFS=$'\t'
	tmp=($1)
	if [ "$tmp[1]" = "" ]
	then
		$tmp[1]=";;;;"
	fi
	RSTR="s/>${tmp[0]}/>${tmp[1]}/g"
	perl -pi -e $RSTR $FA
}

while IFS='' read -r LinefromFile || [[ -n "${LinefromFile}" ]]; do

    process_line "$LinefromFile"

done < "$1"