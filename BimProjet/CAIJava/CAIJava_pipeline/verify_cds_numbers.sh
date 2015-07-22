#!env bash

# path='/home/raphael/Desktop/TME/BimProjet/CAIJava/source/Marine/'
path='/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Marine/'

for sub in `ls $path`; do
    for ssub in `ls "${path}${sub}"`; do
	if [[ -d "${path}${sub}" && ("${path}${sub}/${ssub}" = *"gb"* ) ]] ; then
	    count=`grep -n CDS "${path}${sub}/${ssub}" | wc -l`
	    # if [ "$count" != "0" ]; then
	    echo "       ${sub}/${ssub} : $count"
	    # fi
	elif [[ -d "${path}${sub}" && ("${path}${sub}/${ssub}" = *"fragGeneScan.ffn") ]] ; then
	    count=`grep -n '>' "${path}${sub}/${ssub}" | wc -l`
	    # if [ "$count" != "0" ]; then
	    echo "       ${sub}/${ssub} : $count"
	    # fi
	fi
    done
done
