#!env bash

function sub () {
    for i in {30..30} ; do
	fname="/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/AT_arc_metatrans.filtered.fasta.cleanup.len$i"
	echo $fname
	if [ ! -f $fname ] ; then
	    path=`pwd`
	    cd /home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/Scripts/
	    python filter_seq_length.py $i
	    cd $path
	fi
	outtxt="outtxt_$i"
	echo "Metagenomic $i" > $outtxt
	java -Xmx6g -XX:-UseGCOverheadLimit -jar CAIJava.jar $fname -s -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> $outtxt
	python extract_signature.py $outtxt >> "sign_vec_len$i.csv"
	# sleep 60
    done  
}

sub
