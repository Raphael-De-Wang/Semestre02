#!env bash

basepath='/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Marine/'

hmmfile='/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Pfam/Pfam-A.hmm'

for file in `ls $basepath` ; do
    if [ -d "$basepath$file" ] ; then
	for sfile in `ls "$basepath$file"` ; do
	    if [[ $sfile == *"translate.fasta" ]] ; then
		# hmmscan --domtblout /tmp/Spirodela_polyrhiza.domtblout Pfam/Pfam-A.hmm Marine/Streptophyta/Spirodela_polyrhiza_translate.fasta
		fname="$basepath$file/$sfile"
		echo "hmmscan --domtblout ${fname/_translate.fasta/.domtblout} ${hmmfile} ${fname} > ${fname/_translate.fasta/_scan.output}"
	    fi
	done
    fi
done
