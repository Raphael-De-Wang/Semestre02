#!env bash

species_list="*Chlorophyta-Chlorella_vulgaris
*Chlorophyta-Picochlorum
*Chlorophyta-Trebouxia_gelatinosa
*Ciliophora-Paramecium_biaurelia
*Ciliophora-Paramecium_caudatum
*Ciliophora-Paramecium_sexaurelia
*Ciliophora-Tetrahymena_borealis
*Ciliophora-Tetrahymena_elliotti
*Ciliophora-Tetrahymena_malaccensis
*Rhodophyta-Porphyridium_purpureum
*Streptophyta-Klebsormidium_flaccidum
*Streptophyta-Spirodela_polyrhiza
Apicomplexa-Perkinsus_marinus
Apicomplexa-Toxoplasma_gondii
Bacillariophyta-Phaeodactylum_tricornutum
Bacillariophyta-Thalassiosira_oceanica
Bacillariophyta-Thalassiosira_pseudonana
Chlorophyta-Auxenochlorella_protothecoides
Chlorophyta-Bathycoccus_prasinos
Chlorophyta-Chlamydomonas_reinhardtii
Chlorophyta-Chlorella_variabilis
Chlorophyta-Coccomyxa_subellipsoidea
Chlorophyta-Helicosporidium
Chlorophyta-Micromonas
Chlorophyta-Micromonas_pusilla_CCMP1545
Chlorophyta-Monoraphidium_neglectum
Chlorophyta-Nonlabens_ulvanivorans
Chlorophyta-Ostreococcus_lucimarinus
Chlorophyta-Ostreococcus_tauri
Chlorophyta-Volvox_carteri
Ciliophora-Ichthyophthirius_multifiliis
Ciliophora-Oxytricha_trifallax
Ciliophora-Paramecium_tetraurelia
Ciliophora-Stylonychia_lemnae
Ciliophora-Tetrahymena_thermophila
Cryptophyta-Cryptomonas_Paramecium
Cryptophyta-Guillardia_theta
Cryptophyta-Hemiselmis_andersenii
Haptophyta-Emiliania_huxleyi
Rhodophyta-Chondrus_crispus
Rhodophyta-Cyanidioschyzon_merolae
Rhodophyta-Galdieria_sulphuraria
Streptophyta-Physcomitrella_patens
*Dinoflagellata-Symbiodinium_minutum"

basepath='/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Marine'

# cd /users/Etu9/3404759/Workspace/toolkits/InterProScan5/interproscan-5.13-52.0/

for specie in $species_list ; do
    echo $specie
    IFS='-' read -a array << EOF
$specie
EOF
    if [[ $specie != \** ]] ; then
	specie=${specie/\*/};
	sfname="$basepath/${array[0]}/${array[1]}.gb"
	tfname="$basepath/${array[0]}/${array[1]}_CDS.fasta"
	if [ ! -f $tfname] ; then
	    echo "extract CDS"
	    python gbk_extract_cds.py $sfname
	fi
    else
	tfname="$basepath/${array[0]}/${array[1]}.fsa_nt"
    fi
    echo "interproscan $tfname"
    # ofname="$basepath/${array[0]}/${array[1]}.tsv"
    /users/Etu9/3404759/Workspace/toolkits/InterProScan5/interproscan-5.13-52.0/interproscan.sh -dp -f tsv -appl Pfam --goterms -i $tfname
done
