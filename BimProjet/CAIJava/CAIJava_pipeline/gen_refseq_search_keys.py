#!env python

tgt = '/tmp/txt'

handle = open(tgt)

refseq_list = []

for line in handle:
    # refseq_list.append(line.split()[3])
    print line.split()[3]

# print ' OR '.join(refseq_list)

