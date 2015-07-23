#!env bash

species_list="*Chlorophyta-Picochlorum
*Chlorophyta-Trebouxia_gelatinosa
*Ciliophora-Paramecium_biaurelia
*Ciliophora-Paramecium_caudatum
*Ciliophora-Paramecium_sexaurelia
*Ciliophora-Tetrahymena_borealis
*Ciliophora-Tetrahymena_elliotti
*Ciliophora-Tetrahymena_malaccensis
*Streptophyta-Spirodela_polyrhiza
Chlorophyta-Nonlabens_ulvanivorans
Ciliophora-Ichthyophthirius_multifiliis
Ciliophora-Oxytricha_trifallax
Ciliophora-Paramecium_tetraurelia
Ciliophora-Stylonychia_lemnae
Ciliophora-Tetrahymena_thermophila
Cryptophyta-Cryptomonas_Paramecium
Cryptophyta-Guillardia_theta
Cryptophyta-Hemiselmis_andersenii
Rhodophyta-Galdieria_sulphuraria
Streptophyta-Physcomitrella_patens"

species_list="*Chlorophyta-Chlorella_vulgaris
*Rhodophyta-Porphyridium_purpureum
*Streptophyta-Klebsormidium_flaccidum
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
Chlorophyta-Ostreococcus_lucimarinus
Chlorophyta-Ostreococcus_tauri
Chlorophyta-Volvox_carteri
Haptophyta-Emiliania_huxleyi
Rhodophyta-Chondrus_crispus
Rhodophyta-Cyanidioschyzon_merolae"

basepath='/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Marine'
hmmfile='/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Pfam/Pfam-A.hmm'

gcai_file_list=""

CLUSTER_NUM=$1
CLUSTER_ID=$2

LOOP_COUNT=0

echo "HMMSCAN Pfam"

for specie in $species_list ; do
    LOOP_COUNT=$((LOOP_COUNT+1))
    HIT_NUM=$((LOOP_COUNT%CLUSTER_NUM)) 
    if [ "$HIT_NUM" -ne "$CLUSTER_ID" ] ; then
	continue
    fi
    echo $specie
    IFS='-' read -a array << EOF
${specie/\*/}
EOF
    sfname="$basepath/${array[0]}/${array[1]}.gb"
    rfname="$basepath/${array[0]}/${array[1]}_CDS.fasta"
    tfname="$basepath/${array[0]}/${array[1]}_translate.fasta"
    dfname="$basepath/${array[0]}/${array[1]}.domtblout"
    ofname="$basepath/${array[0]}/${array[1]}_scan.output"
    pfname="$basepath/${array[0]}/${array[1]}_domain.fasta"
    gfname="$basepath/${array[0]}/${array[1]}_gcai.csv"
    if [[ $specie != \** ]] ; then
	specie=${specie/\*/};
	if [ ! -f "$dfname" ] ; then 
	    if [ ! -f "$tfname" ] ; then
		if [ ! -f "$rfname" ] ; then
		    echo "extract CDS"
		    python gbk_extract_cds.py $sfname $rfname
		fi
		echo "translate CDS"
		python simple_translate.py $rfname $tfname
	    fi
	    echo "hmmscan"
	    hmmscan --domtblout $dfname $hmmfile $tfname > $ofname
	fi
	if [ ! -f "$pfname" ] ; then
	    echo "Extract Domain Sequences"
	    python extract_domain_sequences.py $dfname $pfname
	fi

    elif [ ! -f "$pfname" ] ; then
	echo "Extract Domain Sequences"
	python extract_domain_sequences.py $dfname $pfname -orf6
    fi

    if [ ! -f "$gfname" ] ; then
	echo "CALCULATE gCAI"
	python calculate_cai_value.py $pfname ewvalues $gfname
    fi
    gcai_file_list="$gcai_file_list $gfname"
done

if [ ! -f domain_gcai_abundance.csv ] ; then 
    echo "Abundance Sum UP"
    python ordonner_domain_gcai_abundance.py $gcai_file_list
fi

echo "PLOT"
python plot_gcais.py

