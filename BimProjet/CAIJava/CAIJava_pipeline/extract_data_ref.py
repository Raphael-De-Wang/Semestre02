def tfunc(fname):
    flag = False
    handle = open(fname)
    for line in handle:
        if 'DNA' in line:
            flag = True
        elif flag is True:
            print line.split()[0]
            flag = False
    handle.close()

fname = '/home/raphael/Downloads/nuccore_result.txt'

tfunc(fname)
