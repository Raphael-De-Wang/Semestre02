#! env bash

echo 'WRITE CSV TITLE'
echo '"" "tat" "tgt" "ggt" "tct" "ttt" "tgc" "tag" "taa" "tac" "ttc" "tcg" "tta" "ttg" "tcc" "tca" "gca" "gta" "gcc" "gtc" "gcg" "gtg" "cgt" "gtt" "gct" "gat" "ctt" "cct" "cga" "cgc" "ctc" "aca" "cgg" "ggg" "gga" "ggc" "gag" "acg" "gac" "ccg" "gaa" "acc" "atg" "aag" "aaa" "atc" "aac" "ata" "agg" "cag" "agc" "aga" "cat" "aat" "att" "ctg" "cta" "act" "cac" "tga" "caa" "agt" "cca" "ccc" "tgg"' > sign_vec.csv
bpath='/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Marine/'
# bpath='~/Workspace/Semestre02/BimProjet/CAIJava/source/Marine/'


# flist="Apicomplexa.gb  Bacillariophyta.gb  Chlorophyta.gb  Ciliophora.gb  Cryptophyta.gb  Dinoflagellata.gb  Haptophyta.gb  Rhodophyta.gb  Streptophyta.gb"

# flist="Apicomplexa.gb  Bacillariophyta.gb  Ciliophora.gb  Cryptophyta.gb  Haptophyta.gb  Rhodophyta.gb"

# flist="Apicomplexa.gb  Bacillariophyta.gb  Chlorophyta.gb  Ciliophora.gb  Cryptophyta.gb  Haptophyta.gb  Rhodophyta.gb  Streptophyta.gb"

# for file in $flist ; do
#     echo $file >  outtxt
#     java -Xmx6g -XX:-UseGCOverheadLimit -jar CAIJava.jar "$bpath$file" -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> outtxt 
#     python extract_signature.py >> sign_vec.csv
# done

flist='Apicomplexa/Perkinsus_marinus.gb
       Apicomplexa/Toxoplasma_gondii.gb
       Bacillariophyta/Phaeodactylum_tricornutum.gb
       Bacillariophyta/Thalassiosira_oceanica.gb
       Bacillariophyta/Thalassiosira_pseudonana.gb
       Chlorophyta/Auxenochlorella_protothecoides.gb
       Chlorophyta/Bathycoccus_prasinos.gb
       Chlorophyta/Chlamydomonas_reinhardtii.gb
       Chlorophyta/Chlorella_variabilis.gb
       Chlorophyta/Coccomyxa_subellipsoidea.gb
       Chlorophyta/Helicosporidium.gbff
       Chlorophyta/Micromonas.gb
       Chlorophyta/Micromonas_pusilla_CCMP1545.gb
       Chlorophyta/Monoraphidium_neglectum.gb
       Chlorophyta/Nonlabens_ulvanivorans.gb
       Chlorophyta/Ostreococcus_lucimarinus.gb
       Chlorophyta/Ostreococcus_tauri.gb
       Chlorophyta/Volvox_carteri.gb
       Ciliophora/Ichthyophthirius_multifiliis.gb
       Ciliophora/Oxytricha_trifallax.gb
       Ciliophora/Paramecium_tetraurelia.gb
       Ciliophora/Stylonychia_lemnae.gb
       Ciliophora/Tetrahymena_thermophila.gb
       Cryptophyta/Cryptomonas_Paramecium.gb
       Cryptophyta/Guillardia_theta.gb
       Cryptophyta/Hemiselmis_andersenii.gb
       Haptophyta/Emiliania_huxleyi.gb
       Rhodophyta/Chondrus_crispus.gb
       Rhodophyta/Cyanidioschyzon_merolae.gb
       Rhodophyta/Galdieria_sulphuraria.gb
       Streptophyta/Physcomitrella_patens.gb'

flist='Apicomplexa/Perkinsus_marinus.gb
       Apicomplexa/Toxoplasma_gondii.gb'


for file in $flist ; do
    echo "${file/\//-}" >  outtxt
    java -Xmx6g -XX:-UseGCOverheadLimit -jar CAIJava.jar "$bpath$file" -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> outtxt 
    python extract_signature.py >> sign_vec2.csv
done

# echo '"metagenomics" 1.000 1.000 0.844 0.713 1.000 0.141 0.181 1.000 0.222 0.146 0.117 1.000 0.265 0.164 0.838 1.000 0.884 0.200 0.174 0.089 0.167 0.114 1.000 0.816 1.000 0.180 0.798 0.165 0.021 0.035 0.990 0.016 0.782 1.000 0.145 0.150 0.130 0.212 0.066 1.000 0.098 1.000 0.099 1.000 0.132 0.223 0.786 0.208 0.151 0.105 1.000 1.000 1.000 1.000 0.028 0.152 1.000 0.212 0.180 1.000 1.000 1.000 0.191 1.000' >> sign_vec.csv


