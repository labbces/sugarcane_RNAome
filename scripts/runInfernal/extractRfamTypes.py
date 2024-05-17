#!/usr/bin/env python3

import argparse

parser= argparse.ArgumentParser(description='extract Rfam types from rfam.cmscan.tblout')
parser.add_argument('--rfamtypes', metavar='rfam-types.txt', type=str, help='unformatted Rfam families and types', required=True)
parser.add_argument('--cmscan', metavar='rfam.cmscan.tblout', type=str, help='cmscan tblout file', required=True)
parser.add_argument('--output', metavar='rfam.cmscan.tblout.types', type=str, help='name for output file', required=True)
args= parser.parse_args()

rfamtypes = args.rfamtypes
cmscan = args.cmscan
output = args.output

RfamTypes = {}

with open(rfamtypes, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        accession = fields[0]
        id = fields[1]
        type = fields[2].split(";") 
                
        if accession not in RfamTypes:
            RfamTypes[accession] = {'id': id, 'types': type}
#print(RfamTypes)

accession_count = {}

with open(cmscan, 'r') as f:
    for line in f:
        fields = line.strip().split()
        if len(fields) > 2:
            accession = fields[2]
            if accession in accession_count:
                accession_count[accession] += 1
            else:
                accession_count[accession] = 1
#print(accession_count)

rnaomeTypes = {}

# count last type for each key
for accession, count in accession_count.items():
    if accession in RfamTypes:
        last_type = RfamTypes[accession]['types'][-1]
        rnaomeTypes[accession] = {'count': count, 'type': last_type}

# concatenate similar types
type_count = {}

for accession, info in rnaomeTypes.items():
    type = info['type']
    count = info['count']
    if type == 'Gene':
        id = RfamTypes[accession]['id']
        type_count[id] = count
    else:
        if type in type_count:
            type_count[type] += count
        else:
            type_count[type] = count

#for type, total_count in type_count.items():
#    print(f'Family type: {type}, {total_count}')
    
with open(output, 'w', newline='\n') as csvfile:
    csvfile.write('Family type,Count\n')
    for tipo, total_count in type_count.items():
        csvfile.write(f'{tipo},{total_count}\n')

print(f"Rfam family type for {cmscan} saved in {output}")