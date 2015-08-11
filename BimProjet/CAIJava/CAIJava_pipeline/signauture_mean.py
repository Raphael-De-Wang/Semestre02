#!env python

import csv
import sys
import numpy as np

# species_list = [ 'Streptophyta-Physcomitrella_patens','Rhodophyta-Galdieria_sulphuraria','Ciliophora-Stylonychia_lemnae','Ciliophora-Oxytricha_trifallax','Chlorophyta-Nonlabens_ulvanivorans','*Ciliophora-Tetrahymena_elliotti','*Ciliophora-Paramecium_sexaurelia','Cryptophyta-Hemiselmis_andersenii','Ciliophora-Paramecium_tetraurelia','*Ciliophora-Paramecium_biaurelia','*Ciliophora-Paramecium_caudatum','Ciliophora-Ichthyophthirius_multifiliis','Cryptophyta-Guillardia_theta','Cryptophyta-Cryptomonas_Paramecium','Ciliophora-Tetrahymena_thermophila','*Ciliophora-Tetrahymena_malaccensis','*Ciliophora-Tetrahymena_borealis']

'''
species_list_len2 = ['Chlorophyta-Nonlabens_ulvanivorans',
'Chlorophyta-Nonlabens_ulvanivorans',
'*Ciliophora-Tetrahymena_elliotti',
'*Ciliophora-Paramecium_sexaurelia',
'Cryptophyta-Hemiselmis_andersenii',
'Streptophyta-Physcomitrella_patens',
'Rhodophyta-Galdieria_sulphuraria',
'Ciliophora-Stylonychia_lemnae',
'Ciliophora-Oxytricha_trifallax',
'*Ciliophora-Tetrahymena_malaccensis',
'Ciliophora-Tetrahymena_thermophila',
'Ciliophora-Paramecium_tetraurelia',
'*Ciliophora-Paramecium_biaurelia',
'*Ciliophora-Paramecium_caudatum',
'Ciliophora-Ichthyophthirius_multifiliis',
'Cryptophyta-Cryptomonas_Paramecium',
'Cryptophyta-Guillardia_theta',
'*Ciliophora-Tetrahymena_borealis']
'''

def load_sign_dict(inum):
    fname="sign_vec_len%d.csv"%inum
    handle = open(fname)
    reader = csv.reader(handle, delimiter=' ')
    sign_dict = {}
    labels = []
    espece_list = []
    for record in reader:
        if len(record[0]) == 0:
            labels = record[1:]
            continue
        sign_dict[record[0]] = np.array([ float(r) for r in record[1:] ])
        espece_list.append(record[0])
    return (sign_dict,labels,espece_list)
    
def calculate_mean(sign_dict,species_list):
    sum = np.zeros(64)
    for s in species_list:
        sum += sign_dict[s]
    moy=sum/len(species_list)
    return moy

def print_ewvalue_1(labels,values):
    for i,l in enumerate(labels):
        print "%s,%f"%(l,values[i])

def print_ewvalue_2(labels,values):
    print ' '.join([ "%f"%m for m in values ])
    
def sign_vec_to_ewvalue(ename):
    sign_dict,labels = load_sign_dict()
    print_ewvalue_1(labels,sign_dict[ename])

def load_espece_names(inum,espece_list):
    with open("cluster_info.txt") as ihandle:
        for line in ihandle:
            v = line.split(',')
            if v[0] == "len_%d"%inum :
                break
    target_list = [ espece_list[int(i)-1] for i in v[1:] ]
    return target_list

inum = float(sys.argv[1])
sign_dict,labels,all_espece_list = load_sign_dict(inum)
# species_list = load_espece_names(inum,all_espece_list)
# species_list = ['*Ciliophora-Paramecium_sexaurelia']
'''
species_list = ['*Ciliophora-Tetrahymena_malaccensis',
'Ciliophora-Tetrahymena_thermophila',
'Ciliophora-Ichthyophthirius_multifiliis',
'*Ciliophora-Tetrahymena_elliotti',
'*Ciliophora-Tetrahymena_borealis',
'Rhodophyta-Galdieria_sulphuraria',
'Cryptophyta-Cryptomonas_Paramecium',
'*Ciliophora-Paramecium_biaurelia',
'*Ciliophora-Paramecium_sexaurelia',
'*Ciliophora-Paramecium_caudatum',
'Ciliophora-Paramecium_tetraurelia',
'Cryptophyta-Guillardia_theta',
'Ciliophora-Oxytricha_trifallax',
'Ciliophora-Stylonychia_lemnae',
'Streptophyta-Physcomitrella_patens',
'Chlorophyta-Nonlabens_ulvanivorans',
'Cryptophyta-Hemiselmis_andersenii']
'''
species_list = ['Metagenomic 30']

moy = calculate_mean(sign_dict,species_list)
# sign_vec_to_ewvalue("metagenomique")
print_ewvalue_1(labels,moy)
# print_ewvalue_2(labels,moy)

