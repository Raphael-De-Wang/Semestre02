#!env python

import csv
import sys
import operator

gcai_dict = {}
abundance_list = []

def simplify_domain_name(dname):
    return dname.split('.')[0]
    

for arg in sys.argv[1:]:
    with open(arg) as gcai_handler:
        reader = csv.reader(gcai_handler, delimiter='\t')
        for record in reader:
            dname = record[0].split('|')[2]
            if cmp(dname[:2],'PF') <> 0:
                raise ValueError("Irregular Domain Name.")
            dname = simplify_domain_name(dname)
            gcai  = float(record[1])
            if not gcai_dict.has_key(dname):
                gcai_dict[dname] = [gcai]
            else:
                gcai_dict[dname].append(gcai)
                
for key,value in gcai_dict.iteritems():
    abundance_list.append([key,len(value)])

abundance_list = sorted(abundance_list, key=operator.itemgetter(1),reverse=True)

csv_handle = open("domain_gcai_abundance.csv",'w')
csv_writer = csv.writer(csv_handle, delimiter='\t')

for dname,abundance in abundance_list:
    dname = simplify_domain_name(dname)
    for gcai in gcai_dict[dname]:
        csv_writer.writerow([dname,gcai,abundance])

csv_handle.close()            
