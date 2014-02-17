#!/usr/bin/env python
import os
import re
import sys

template = open(sys.argv[1],'r').readlines()
#print '###########', template

bashfile = open('batch.sh', 'w')
bashfile.write('#! /bin/bash\necho "Starting..."\n')

j = 1
datfiles = [x.split('data.')[1] for x in os.listdir('.') if 'data.' in x]
for datname in datfiles:
    in_name = 'in.{}'.format(datname)
    fout = open(in_name, 'w') 
    bashfile.write('echo "Using file {}" \n'.format(datname))

    # replace occurences of 'tmplate' with desired string
    content = template[:]
    for i in range(len(template)):
        if 'tmplate' in template[i] :
            str_ = '{}'.format(datname)
            content[i] = str_.join(content[i].split('tmplate'))
        fout.write(content[i])
    fout.close()

    bashfile.write('lmp_serial < {} | tee log.{} \n'.format(in_name, datname))
    bashfile.write('echo "({}/{}) done" \n'.format(j, len(datfiles)))
    j += 1

bashfile.write('echo "...aaand we\'re done here!"\n')
bashfile.close()
print 'Best of luck!'
