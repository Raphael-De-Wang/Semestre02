#!env python

import csv

handle = open('Rhodophyta-Chondrus_crispus_sig_vec.csv')
reader = csv.reader(handle,delimiter=' ')
for row in reader:
    labels = row
    break
for row in reader:
    values = row
    break
    
handle.close()
wvalues = {}

for i,l in enumerate(labels):
    if i == 0:
        continue
    # wvalues[l] = values[i]
    print '%s\t ,%s'%(l,values[i])
    
exit()
handle = open('../output/wvalues')

reader = csv.reader(handle,delimiter='\t')
for row in reader:
    # lv = row[:2]+[wvalues[row[1]]]
    lv = [row[1]]+[wvalues[row[1]]]
    line = '\t ,'.join(lv)
    print line
