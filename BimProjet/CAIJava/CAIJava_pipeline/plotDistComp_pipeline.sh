#!env bash

function sub() {
    for i in {30..30} ; do
	echo "len $i"
	cd /home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/CAIJava_pipeline
	ewvalue="ewvalue_len$i"
	fasta="../source/AT_arc_metatrans.filtered.fasta.cleanup.len$i"
	cais="../output/cais_len$i.lst"
	ex_cais="../output/cais_ewvalue_len$i.lst"
	if [ ! -f $fasta ] ; then
	    echo "filtrer les metagenomes"
	    path=`pwd`
	    cd /home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/Scripts/
	    python filter_seq_length.py $i
	    cd $path
	fi
	echo "signauture mean"
	python signauture_mean.py "$i" > $ewvalue
	echo "calculate cai value by ewvalue"
	python calculate_cai_value.py "$fasta" "$ewvalue" "$ex_cais"
	# echo "calculate cai value"
	# java -Xmx6g -XX:-UseGCOverheadLimit -jar CAIJava.jar $fasta -s -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g > /tmp/outtxt
	# ls -al cais.lst 
	# mv cais.lst $cais
	cd /home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/Scripts/
	python plotDistComparaison.py $i
	# sleep 300
    done 
}

sub
