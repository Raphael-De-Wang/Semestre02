#!env python

handle = open('/home/raphael/Downloads/genome_result_5.txt')
i = 0
flag = False
for line in handle:
    if '.' in line:
        sep = line.split('.')
        try:
            num = int(sep[0])
            name= sep[1]
        except:
            continue
    if 'Chromosome' in line:
        compl = line.split(':')[-1]
        if 'no data' in line:
            # flag = True
            i += 1
        continue
    if 'Genome ID' in line:
        # if flag :
        # print '    ',num,name.strip(),' - ',compl.strip()
        # else:
        print '    ',num,name.strip(),' - ',compl.strip(),' - ',line.strip()
        # flag = False
        
print i, 'records have no data'
