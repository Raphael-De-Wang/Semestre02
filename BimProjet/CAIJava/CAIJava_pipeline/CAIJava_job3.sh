#! env bash

echo 'WRITE CSV TITLE'
echo '"" "tat" "tgt" "ggt" "tct" "ttt" "tgc" "tag" "taa" "tac" "ttc" "tcg" "tta" "ttg" "tcc" "tca" "gca" "gta" "gcc" "gtc" "gcg" "gtg" "cgt" "gtt" "gct" "gat" "ctt" "cct" "cga" "cgc" "ctc" "aca" "cgg" "ggg" "gga" "ggc" "gag" "acg" "gac" "ccg" "gaa" "acc" "atg" "aag" "aaa" "atc" "aac" "ata" "agg" "cag" "agc" "aga" "cat" "aat" "att" "ctg" "cta" "act" "cac" "tga" "caa" "agt" "cca" "ccc" "tgg"' > sign_vec.csv

basepath='/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Marine/'

for file in `ls $basepath` ; do
    if [ -d "$basepath$file" ] ; then
	for sfile in `ls "$basepath$file"` ; do
	    if [[ $sfile == *".fragGeneScan.ffn" ]] ; then
		fname="$basepath$file/$sfile"
		echo $fname
		echo "*$file-${sfile/.fragGeneScan.ffn/}" > outtxt
		java -Xmx6g -XX:-UseGCOverheadLimit -jar CAIJava.jar $fname -s -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> outtxt
		python extract_signature.py >> sign_vec.csv
	    fi
	done
    fi
done


