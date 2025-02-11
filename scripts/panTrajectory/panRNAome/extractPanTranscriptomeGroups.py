#!/usr/bin/env python3

import pandas as pd
import numpy as np
from math import factorial
import sys
import argparse
import os.path
from plotnine import ggplot, aes, geom_jitter, geom_smooth, theme_bw, scale_y_log10, labs, save_as_pdf_pages, theme

# creating arguments
parser = argparse.ArgumentParser(prog = "extractPanTranscriptomeGroups", description="Extract pan, hard-core, soft-core, accessory, and exclusives orthogroups from GeneCount.tsv file generated by OrthoFinder2", add_help=True)
parser.add_argument("-i", dest="input_file", metavar="<Orthogroups.GeneCount.tsv>", help="OrthoFinde2r GeneCount output", required=True)
parser.add_argument("-p", dest="min_fraction", type=float, help="Minimum fraction of genotypes to be present to consider a group as sort-core", default=0.9)

# getting arguments
args = parser.parse_args()
input_file = args.input_file
min_fraction=args.min_fraction
output_file = "Pan-Transcriptome_Size_" + str(min_fraction) + "." + os.path.basename(input_file)
output_pan = "Pan-Transcriptome_Size_" + str(min_fraction) + "." + os.path.basename(input_file)
output_core = "Core-Transcriptome_Size_" + str(min_fraction) + "." + os.path.basename(input_file)
output_accessory = "Accessory-Transcriptome_Size_" + str(min_fraction) + "." + os.path.basename(input_file)
output_exclusive = "Exclusive-Transcriptome_Size_" + str(min_fraction) + "." + os.path.basename(input_file)
outfigure_Groups_pdf="Pan-Transcriptome_Trajectory_Groups_" + str(min_fraction) + "." + os.path.basename(input_file) + ".pdf"
outfigure_pdf="Pan-Transcriptome_Trajectory_" + str(min_fraction) + "." + os.path.basename(input_file) + ".pdf"
outfigure_Groups_png="Pan-Transcriptome_Trajectory_Groups_" + str(min_fraction) + "." + os.path.basename(input_file) + ".png"
outfigure_Genes_pdf="Pan-Transcriptome_Trajectory_Genes_" + str(min_fraction) + "." + os.path.basename(input_file) + ".pdf"
outfigure_Genes_png="Pan-Transcriptome_Trajectory_Genes_" + str(min_fraction) + "." + os.path.basename(input_file) + ".png"

# setting recursion limit to 1000000
sys.setrecursionlimit(10**6)

data = pd.read_csv(input_file, delimiter='\t', header=0, index_col=0)

print(f'The dimensions of the input dataset are: {data.shape}')
if 'Total' in data.columns:
    sys.stderr.write("The input data has a column names \"Total\", this is not expected and can generate some problems. Column will be removed!")
    data.drop('Total', axis=1, inplace=True)
    print(f'The dimensions of the input dataset, after removing the Total column, are: {data.shape}')

print(f'The dimensions of the input dataset, after removing the Total column, are: {data.shape}')

def sample_random_selection(data, samples):
    my_sample = data.sample(number_genotypes, axis='columns', replace=False)
    
    my_selection = sorted(list(my_sample.columns))
    if my_selection in samples:
        sample_random_selection(data, samples)
    else:
        samples.append(my_selection)
    return(my_sample, samples)

with open(output_file, "w") as write_output_file:
    write_output_file.write("NumberGroups\tNumberGenes\tClassification\tNumberGenotypes\tSample\n")

    for number_genotypes in range(1, data.shape[1]):
        max_n_sample = 20
        samples = []

        max_number_of_samples = int(factorial(data.shape[1]) / factorial((data.shape[1] - number_genotypes)))
        if max_number_of_samples < max_n_sample:
            max_n_sample = max_number_of_samples

        for n_sample in range(0, max_n_sample):
            pan_transcriptome_size = 0
            genes_pan = 0
            hard_core_transcriptome_size = 0
            soft_core_transcriptome_size = 0
            accessory_transcriptome_size = 0
            exclusive_transcriptome_size = 0
            genes_hard_core = 0
            genes_soft_core = 0
            genes_accessory = 0
            genes_exclusive = 0
            my_sample,samples = sample_random_selection(data, samples)
            my_sample = np.array(my_sample)

            for orthogroup in my_sample:
                number_of_genotypes_in_orthogroups = 0
                if sum(orthogroup) > 0:
                    genes_pan += sum(orthogroup)
                    pan_transcriptome_size += 1
                    for genotype in orthogroup:
                        if genotype > 0:
                            number_of_genotypes_in_orthogroups += 1

                    proporcao = number_of_genotypes_in_orthogroups / number_genotypes
                    if proporcao >= min_fraction:
                        genes_soft_core += sum(orthogroup)
                        soft_core_transcriptome_size += 1
                        if number_of_genotypes_in_orthogroups == number_genotypes:
                            genes_hard_core += sum(orthogroup)
                            hard_core_transcriptome_size += 1
                    elif number_of_genotypes_in_orthogroups > 1 and proporcao < min_fraction:
                        genes_accessory += sum(orthogroup)
                        accessory_transcriptome_size += 1
                    elif number_of_genotypes_in_orthogroups == 1:
                        genes_exclusive += sum(orthogroup)
                        exclusive_transcriptome_size += 1
                    else:
                        print("Case not contemplated, check input\n")

            write_output_file.write(str(pan_transcriptome_size)+ "\t"
                                    + str(genes_pan) + "\t"
                                    + "pan" + "\t"
                                    + str(number_genotypes) + "\t"
                                    + str(n_sample) + "\n"
                                    + str(hard_core_transcriptome_size) + "\t"
                                    + str(genes_hard_core) + "\t"
                                    + "hard-core" + "\t"
                                    + str(number_genotypes) + "\t"
                                    + str(n_sample) + "\n"
                                    + str(soft_core_transcriptome_size) + "\t"
                                    + str(genes_soft_core) + "\t"
                                    + "soft-core" + "\t"
                                    + str(number_genotypes) + "\t"
                                    + str(n_sample) + "\n"
                                    + str(accessory_transcriptome_size) + "\t"
                                    + str(genes_accessory) + "\t"
                                    + "accessory" + "\t"
                                    + str(number_genotypes) + "\t"
                                    + str(n_sample) + "\n"
                                    + str(exclusive_transcriptome_size) + "\t"
                                    + str(genes_exclusive) + "\t"
                                    + "exclusive" + "\t"
                                    + str(number_genotypes) + "\t"
                                    + str(n_sample) + "\n")


data = pd.read_csv(output_file, delimiter="\t")
title_text_groups="Pan-Transcriptome Trajectory " + str(min_fraction) + " -- Groups"
title_text_genes ="Pan-Transcriptome Trajectory " + str(min_fraction) + " -- Genes"

fig1=ggplot(data) + aes(x="NumberGenotypes", y="NumberGroups", colour="Classification") + theme_bw() + geom_jitter() + geom_smooth(method="loess") + scale_y_log10() + labs(title=title_text_groups, x="Number of Genotypes", y="Number of Groups") + theme(figure_size=(20, 13))
fig1.save(outfigure_Groups_png, dpi=600, height=13, width=20, units = 'in')

fig2=ggplot(data) + aes(x="NumberGenotypes", y="NumberGenes", colour="Classification") + theme_bw() + geom_jitter() + geom_smooth(method="loess") + scale_y_log10() + labs(title=title_text_genes, x="Number of Genotypes", y="Number of Genes") + theme(figure_size=(20, 13))

fig2.save(outfigure_Genes_png, dpi=600, height=13, width=20, units = 'in')

plots=[fig1, fig2]
save_as_pdf_pages(plots,filename=outfigure_pdf)
