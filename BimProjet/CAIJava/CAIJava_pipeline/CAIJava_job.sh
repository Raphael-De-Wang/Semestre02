#! env bash

# echo 'CLEANUP...'
# rm tmpSource/*

# echo 'PARSE GENBANK'
# python parse_genbank.py > source_accession.csv

# echo 'REMOVE BAD GENBANK'
# mv tmpSource/NC_007405*gbk tmpSource/NC_007405.gbk.bak
# mv tmpSource/NC_013710*gbk tmpSource/NC_013710.gbk.bak
# mv tmpSource/NC_016739*gbk tmpSource/NC_016739.gbk.bak
# mv tmpSource/NC_027265*gbk tmpSource/NC_027265.gbk.bak
# mv tmpSource/NC_026126*gbk tmpSource/NC_026126.gbk.bak

echo 'WRITE CSV TITLE'
echo '"" "tat" "tgt" "ggt" "tct" "ttt" "tgc" "tag" "taa" "tac" "ttc" "tcg" "tta" "ttg" "tcc" "tca" "gca" "gta" "gcc" "gtc" "gcg" "gtg" "cgt" "gtt" "gct" "gat" "ctt" "cct" "cga" "cgc" "ctc" "aca" "cgg" "ggg" "gga" "ggc" "gag" "acg" "gac" "ccg" "gaa" "acc" "atg" "aag" "aaa" "atc" "aac" "ata" "agg" "cag" "agc" "aga" "cat" "aat" "att" "ctg" "cta" "act" "cac" "tga" "caa" "agt" "cca" "ccc" "tgg"' > sign_vec.csv

# echo 'WRITE SIGNATURE'
# for gbk in `ls tmpSource/*gbk` ; do
    # java -jar /home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/CAIJava.jar $gbk -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g
#     echo $gbk > outtxt
#     java -jar CAIJava.jar $gbk -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> outtxt
#     python extract_signature.py >> sign_vec.csv
# done

# echo 'WRITE BACILLORIAPHYTA SIGNATURE'
# echo 'Bacillariophyta.fasta' > outtxt
# java -jar CAIJava.jar ../source/Bacillariophyta.fasta -s -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> outtxt 
# python extract_signature.py >> sign_vec.csv

echo 'WRITE THALASSIOSIRA PSEUDONANA SIGNATURE'
echo 'Thalassiosira_pseudonana.gb' > outtxt
java -jar CAIJava.jar ../source/Bacillariophyta/Thalassiosira_pseudonana/Thalassiosira_pseudonana.gb -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> outtxt 
python extract_signature.py >> sign_vec.csv

echo 'WRITE PHAEODACTYLUM_TRICORNUTUM SIGNATURE'
echo 'Phaeodactylum_tricornutum.gb' > outtxt
java -jar CAIJava.jar ../source/Bacillariophyta/Phaeodactylum_tricornutum/Phaeodactylum_tricornutum.gb -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> outtxt 
python extract_signature.py >> sign_vec.csv

echo 'WRITE THALASSIOSIRA OCEANICA SIGNATURE'
echo 'Thalassiosira_oceanica.gb' > outtxt
java -jar CAIJava.jar ../source/Bacillariophyta/Thalassiosira_oceanica/Thalassiosira_oceanica.gb -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> outtxt
python extract_signature.py >> sign_vec.csv

echo 'WRITE SIGNATURE METAGENOMIC'
echo 'metagenomics.gbk' > outtxt
# java -jar CAIJava.jar ../source/AT_arc_metatrans.filtered.fasta.cleanup -s -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> outtxt
java -jar CAIJava.jar ../source/AT_arc_metatrans.filtered.fasta -s -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> outtxt 
python extract_signature.py >> sign_vec.csv

echo 'R PLOT'
R -f heatmap.R >> /dev/null
