#!env bash

path='/home/raphael/Desktop/TME/BimProjet/CAIJava/source/Marine/'

for sub in `ls $path`; do
    for ssub in `ls "${path}${sub}"`; do
	count=`grep -n CDS "${path}${sub}/${ssub}" | wc -l`
	if [ "$count" = "0" ]; then
	    echo "${path}${sub}/${ssub}"
	    sleep 2
	fi
    done
done
