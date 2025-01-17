#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Reformat genotype names in the ClassificationTable.tsv generated by computeOrthogroupStats.py",
                    epilog="This script replaces underscores (_) with hyphens (-) in genotype names within the ClassificationTable.tsv file")
parser.add_argument("-i", dest="infile", metavar="ClassificationTable.tsv", required=True)
parser.add_argument("-o", dest="outfile", metavar="reformatedClassificationTable.tsv", required=True)
parser.add_argument("-s", dest="separator", metavar="[ - or _ ]", help="Specify the desired separator for genotype names (use - for hyphen or _ for underscore)", required=True)

args = parser.parse_args()
infile = args.infile
outfile = args.outfile
separator = args.separator

print("Detected parameters:","\n")
print(f"infile: {infile}")
print(f"outfile: {outfile}")
print(f"separator: {separator}","\n")

# Dictionary for mapping incorrect genotypes to correct ones (used by replace_genotypes_with_hyphen)
genotype_mapping = {
    "CP74_2005": "CP74-2005",
    "F1_Bulk_1": "F1-Bulk-1",
    "F1_Bulk_2": "F1-Bulk-2",
    "FN95_1702": "FN95-1702",
    "GT96_167": "GT96-167",
    "GXU_34140": "GXU-34140",
    "GXU_34176": "GXU-34176",
    "KQ08_2850": "KQ08-2850",
    "KQB07_23863": "KQB07-23863",
    "KQB07_23990": "KQB07-23990",
    "KQB07_24619": "KQB07-24619",
    "KQB07_24739": "KQB07-24739",
    "KQB08_32953": "KQB08-32953",
    "KQB09_20432": "KQB09-20432",
    "KQB09_20620": "KQB09-20620",
    "KQB09_23137": "KQB09-23137",
    "LCP85_384": "LCP85-384",
    "MT11_610": "MT11-610",
    "QA02_1009": "QA02-1009",
    "QA96_1749": "QA96-1749",
    "QBYN04_26041": "QBYN04-26041",
    "QC02_402": "QC02-402",
    "QN05_1460": "QN05-1460",
    "QN05_1509": "QN05-1509",
    "QN05_1743": "QN05-1743",
    "QN05_803": "QN05-803",
    "QS99_2014": "QS99-2014",
    "SP80_3280": "SP80-3280",
}

# Replace (_) with (-) in genotype names
def replace_genotypes_with_hyphen(data, mapping):
    for line in data:
        parts = line.split("\t")
        old_genotype = parts[2].strip()  # Use the third column directly as the genotype
        for incorrect, correct in mapping.items():
            old_genotype = old_genotype.replace(incorrect, correct)
        parts[2] = old_genotype
        yield "\t".join(parts) + "\n"

# Replace (-) with (_) in genotype names
def replace_genotypes_with_underscore(data, mapping):
    for line in data:
        parts = line.split("\t")
        old_genotype = parts[2].strip().split('_')[0]  # Extracts the genotype part before the underscore
        if old_genotype in mapping:
            # Replace the old genotype with the correct one
            new_genotype = parts[2].replace(old_genotype, mapping[old_genotype])
            parts[2] = new_genotype
        yield "\t".join(parts)

# Invert dictionary for mapping incorrect genotypes to correct ones (used by replace_genotypes_with_underscore) 
def invert_dict(d):
    return {v: k for k, v in d.items()}

with open(infile, "r") as input_file:
    data = input_file.readlines()

if separator == "-":
    replace_function = replace_genotypes_with_hyphen
    mapping = genotype_mapping
    print("Replacing underscore with hyphen ...", "\n")
elif separator == "_":
    replace_function = replace_genotypes_with_underscore
    mapping = invert_dict(genotype_mapping)
    print("Replacing hyphen with underscore ...", "\n")
else:
    print("Invalid separator. Please use - or _.", "\n")
    exit(0)

# Apply genotype replacement
new_data = list(replace_function(data, mapping))

# Write the modified data to the same file or a new file
with open(outfile, "w") as output_file:
    output_file.writelines(new_data)

print("Done!","\n")
