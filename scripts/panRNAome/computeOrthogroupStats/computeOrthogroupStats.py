#!/usr/bin/env python3

import argparse
from genericpath import exists
import os
import matplotlib.pyplot as plt

parser= argparse.ArgumentParser(description='compute stats from orthogroups.tsv')
parser.add_argument('--orthogroupsFile', metavar='orthogroupsFile', type=str, help='file with orthogroups.tsv generated with OrthoFinder2', required=True)
parser.add_argument('--numberSpecies', metavar='numberSpecies', type=int, help='Number of species', required=True)
parser.add_argument('--suffixOut', metavar='suffixOut', type=str, help='Suffix used to create output files', required=True)
args= parser.parse_args()


#Set grlobal variables
orthogroupsFile=args.orthogroupsFile
numberSpecies=args.numberSpecies
suffixOut=args.suffixOut
outputClassificationTableFile='panTranscriptomeClassificationTable_'+ args.suffixOut + '.tsv'
outputDistributionSizeOrthogroupsTableFile='panTranscriptomeDistributionSizeOrthogroupsTable_'+ args.suffixOut + '.tsv'
outputStatisticsOrthogroupsTableFile='panTranscriptomeStatisticsOrthogroups_'+ args.suffixOut + '.tsv'
histogramsOutPNG='panTranscriptomeDistributionSizeOrthogroupsTable_'+ args.suffixOut + '.png'
histogramsOutSVG='panTranscriptomeDistributionSizeOrthogroupsTable_'+ args.suffixOut + '.svg'
numberCoreOrthogroups=0
numberSoftCoreOrthogroups=0
numberAccessoryOrthogroups=0
numberExclusiveOrthogroups=0
numberCoreProteins   =0
numberSoftCoreProteins=0
numberAccessoryProteins=0
numberExclusiveProteins=0
averageProteinsPerOrthogroup=0
averageProteinsPerOrthogroupDict={}
coreOrthogroups={}
softCoreOrthogroups={}
accessoryOrthogroups={}
exclusiveOrthogroups={}

def compositionOrthogroup (data, classification,og,fh):
    for line in data:
        for id in line.split(','):
            id=id.replace(' ','')
            if id != '':
                fh.write(f'{classification}\t{og}\t{id}\n')

def distributionSizeOrthogroups (data,fh):
    for og in data:
        fh.write(f'{og}\t{data[og]}\n')

#Process orthogroupsFile
if os.path.isfile(orthogroupsFile):
    with open(orthogroupsFile, "r") as file, open(outputClassificationTableFile, "w") as outClass, open(outputDistributionSizeOrthogroupsTableFile, "w") as outDist:
        for line in file:
            line=line.rstrip()
            if line.startswith('Orthogroup'):
                header=line.split('\t')
            else:
                numberSpeciesInOrthogroup=0
                numberProteinsInOrthogroup=0
                fields=line.split('\t')
                for i in range(1,len(fields)):
                    if fields[i] != '':
                        numberSpeciesInOrthogroup=numberSpeciesInOrthogroup+1
                        numberProteinsInOrthogroup=numberProteinsInOrthogroup+len(fields[i].split(','))
                        # print(f'{fields[0]}\t{header[i]}\t{fields[i]}')
                averageProteinsPerOrthogroup=numberProteinsInOrthogroup/numberSpeciesInOrthogroup
                averageProteinsPerOrthogroupDict[fields[0]]=averageProteinsPerOrthogroup
                if numberSpeciesInOrthogroup == numberSpecies:
                    #Hard-core groups
                    # print(f'Hard-core OG:{fields[0]}')
                    numberCoreOrthogroups=numberCoreOrthogroups+1
                    numberCoreProteins=numberCoreProteins+numberProteinsInOrthogroup
                    coreOrthogroups[fields[0]]=numberProteinsInOrthogroup
                    compositionOrthogroup(fields[1:numberSpecies+1],"Hard-core",fields[0],outClass)
                if numberSpeciesInOrthogroup >= numberSpecies*0.9:
                    #Soft-core groups
                    # print(f'Soft-core OG:{fields[0]}')
                    numberSoftCoreOrthogroups=numberSoftCoreOrthogroups+1
                    numberSoftCoreProteins=numberSoftCoreProteins+numberProteinsInOrthogroup
                    softCoreOrthogroups[fields[0]]=numberProteinsInOrthogroup
                    compositionOrthogroup(fields[1:numberSpecies+1],"Soft-core",fields[0],outClass)
                if numberSpeciesInOrthogroup > 1 and numberSpeciesInOrthogroup < numberSpecies*0.9 :
                    #Accesory groups
                    # print(f'Accesory-core OG:{fields[0]}')
                    numberAccessoryOrthogroups=numberAccessoryOrthogroups+1
                    numberAccessoryProteins=numberAccessoryProteins+numberProteinsInOrthogroup
                    accessoryOrthogroups[fields[0]]=numberProteinsInOrthogroup
                    compositionOrthogroup(fields[1:numberSpecies+1],"Accessory",fields[0],outClass)
                if numberSpeciesInOrthogroup == 1:
                    #Exclusive groups
                    # print(f'Exclusive OG:{fields[0]}')
                    numberExclusiveOrthogroups=numberExclusiveOrthogroups+1
                    numberExclusiveProteins=numberExclusiveProteins+numberProteinsInOrthogroup
                    exclusiveOrthogroups[fields[0]]=numberProteinsInOrthogroup
                    compositionOrthogroup(fields[1:numberSpecies+1],"Exclusive",fields[0],outClass)
        distributionSizeOrthogroups(averageProteinsPerOrthogroupDict,outDist)
        plt.figure(figsize=[10,8])
        plt.hist(averageProteinsPerOrthogroupDict.values(), bins=100)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Average number of proteins per orthogroup',fontsize=15)
        plt.ylabel('Frequency',fontsize=15)
        plt.savefig(histogramsOutSVG, format="svg")
        plt.savefig(histogramsOutPNG, format="png")


numberCoreProteins=sum(coreOrthogroups.values())

with open(outputStatisticsOrthogroupsTableFile, "w") as outStats:
    outStats.write(f'Number of groups present in all ({numberSpecies}) species:\t{numberCoreOrthogroups}\n')
    outStats.write(f'Number of proteins present in core groups:\t{numberCoreProteins}\n')
    outStats.write(f'Number of groups present in 90% of the species ({numberSpecies*0.9}):\t{numberSoftCoreOrthogroups}\n')
    outStats.write(f'Number of proteins present in soft-core groups:\t{numberSoftCoreProteins}\n')
    outStats.write(f'Number of accesory groups:\t{numberAccessoryOrthogroups}\n')
    outStats.write(f'Number of proteins present in accessory groups:\t{numberAccessoryProteins}\n')
    outStats.write(f'Number of exclusive groups:\t{numberExclusiveOrthogroups}\n')
    outStats.write(f'Number of proteins present in exclusive groups:\t{numberExclusiveProteins}\n')
