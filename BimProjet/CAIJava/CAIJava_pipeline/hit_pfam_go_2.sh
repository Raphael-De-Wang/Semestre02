#!env bash

# CLUSTER_NUM=$1
# CLUSTER_ID=$2

basepath='/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Marine'

pipeline_process(){
    species_list=$1
    CLUSTER_NUM=$2
    CLUSTER_ID=$3
    hmmfile='/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Pfam/Pfam-A.hmm'
    gcai_file_list=""
    LOOP_COUNT=0
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
	# if [ ! -f "$gfname" ] ; then
	#     echo "CALCULATE gCAI"
	#     python calculate_cai_value.py $pfname ewvalues $gfname
	# fi
	# gcai_file_list="$gcai_file_list $gfname"
    done
    # if [ ! -f domain_gcai_abundance.csv ] ; then 
    #     echo "Abundance Sum UP"
    #     python ordonner_domain_gcai_abundance.py $gcai_file_list
    # fi
    # echo "PLOT"
    # python plot_gcais.py
}

pipeline_process_simple(){
    species_list_meta="*Streptophyta-Spirodela_polyrhiza
Rhodophyta-Chondrus_crispus
Bacillariophyta-Thalassiosira_oceanica
Haptophyta-Emiliania_huxleyi
Chlorophyta-Micromonas_pusilla_CCMP1545
Bacillariophyta-Phaeodactylum_tricornutum
Chlorophyta-Coccomyxa_subellipsoidea
Chlorophyta-Auxenochlorella_protothecoides
Chlorophyta-Helicosporidium
*Chlorophyta-Chlorella_vulgaris
Chlorophyta-Chlamydomonas_reinhardtii
*Chlorophyta-Trebouxia_gelatinosa
Bacillariophyta-Thalassiosira_pseudonana
Apicomplexa-Perkinsus_marinus
Chlorophyta-Chlorella_variabilis
*Streptophyta-Klebsormidium_flaccidum
Chlorophyta-Volvox_carteri
Rhodophyta-Cyanidioschyzon_merolae
*Rhodophyta-Porphyridium_purpureum
Chlorophyta-Ostreococcus_tauri
Chlorophyta-Ostreococcus_lucimarinus
Chlorophyta-Monoraphidium_neglectum"

    species_list_2="*Ciliophora-Tetrahymena_malaccensis
Ciliophora-Tetrahymena_thermophila
Ciliophora-Ichthyophthirius_multifiliis
*Ciliophora-Tetrahymena_elliotti
*Ciliophora-Tetrahymena_borealis
Rhodophyta-Galdieria_sulphuraria
Cryptophyta-Cryptomonas_Paramecium
*Ciliophora-Paramecium_biaurelia
*Ciliophora-Paramecium_sexaurelia
*Ciliophora-Paramecium_caudatum
Ciliophora-Paramecium_tetraurelia
Cryptophyta-Guillardia_theta
Ciliophora-Oxytricha_trifallax
Ciliophora-Stylonychia_lemnae
Streptophyta-Physcomitrella_patens
Chlorophyta-Nonlabens_ulvanivorans
Cryptophyta-Hemiselmis_andersenii"

    CLUSTER_NUM=$1
    CLUSTER_ID=$2
    gcai_file_list=""
    LOOP_COUNT=0
    for specie in $species_list_2 ; do
	LOOP_COUNT=$((LOOP_COUNT+1))
	HIT_NUM=$((LOOP_COUNT%CLUSTER_NUM)) 
	if [ "$HIT_NUM" -ne "$CLUSTER_ID" ] ; then
	    continue
	fi
	echo $specie
	IFS='-' read -a array << EOF
${specie/\*/}
EOF
	pfname="$basepath/${array[0]}/${array[1]}_domain.fasta"
	gfname="$basepath/${array[0]}/${array[1]}_gcai.csv"
	echo "CALCULATE gCAI"
	python calculate_cai_value.py $pfname ewvalues_len30 $gfname
	gcai_file_list="$gcai_file_list $gfname"
    done
    
    echo "Abundance Sum UP"
    python ordonner_domain_gcai_abundance.py $gcai_file_list
    
    echo "PLOT"
    python plot_gcais.py
}

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
Streptophyta-Physcomitrella_patens
*Chlorophyta-Chlorella_vulgaris
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
Chlorophyta-Micromonas_pusilla_CCMP1545
Chlorophyta-Monoraphidium_neglectum
Chlorophyta-Ostreococcus_lucimarinus
Chlorophyta-Ostreococcus_tauri
Chlorophyta-Volvox_carteri
Haptophyta-Emiliania_huxleyi
Rhodophyta-Chondrus_crispus
Rhodophyta-Cyanidioschyzon_merolae"

# pipeline_process "$species_list_2" $1 $2 

calculate_signatures(){
    echo 'WRITE CSV TITLE'
    # echo '"" "tat" "tgt" "ggt" "tct" "ttt" "tgc" "tag" "taa" "tac" "ttc" "tcg" "tta" "ttg" "tcc" "tca" "gca" "gta" "gcc" "gtc" "gcg" "gtg" "cgt" "gtt" "gct" "gat" "ctt" "cct" "cga" "cgc" "ctc" "aca" "cgg" "ggg" "gga" "ggc" "gag" "acg" "gac" "ccg" "gaa" "acc" "atg" "aag" "aaa" "atc" "aac" "ata" "agg" "cag" "agc" "aga" "cat" "aat" "att" "ctg" "cta" "act" "cac" "tga" "caa" "agt" "cca" "ccc" "tgg"' > sign_vec.csv
    species_list=$1
    CLUSTER_NUM=$2
    CLUSTER_ID=$3
    LOOP_COUNT=0
    sign_vec="sign_vec$CLUSTER_ID.csv"
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
	fname="$basepath/${array[0]}/${array[1]}_domain.fasta"
	outtxt="outtxt$CLUSTER_ID"
	echo "$outtxt"
	echo $specie > $outtxt
	java -Xmx6g -XX:-UseGCOverheadLimit -jar CAIJava.jar $fname -s -f /tmp/out.dat -t product -id locus_tag -i 15 -k 3 -g >> $outtxt
	python extract_signature.py $outtxt >> $sign_vec
    done
    cat $sign_vec >> sign_vec.csv
    mv $sign_vec /tmp/
    mv $outtxt /tmp/
}

# calculate_signatures "$species_list_1" $1 $2
# calculate_signatures "$species_list_2" $1 $2

pipeline_process_simple $1 $2