#!/usr/bin/env python3

# Dictionary for mapping incorrect genotypes to correct ones
genotype_mapping = {
    "CP74-2005": "CP74_2005",
    "F1-Bulk-1": "F1_Bulk_1",
    "F1-Bulk-2": "F1_Bulk_2",
    "FN95-1702": "FN95_1702",
    "GT96-167": "GT96_167",
    "GXU-34140": "GXU_34140",
    "GXU-34176": "GXU_34176",
    "KQ08-2850": "KQ08_2850",
    "KQB07-23863": "KQB07_23863",
    "KQB07-23990": "KQB07_23990",
    "KQB07-24619": "KQB07_24619",
    "KQB07-24739": "KQB07_24739",
    "KQB08-32953": "KQB08_32953",
    "KQB09-20432": "KQB09_20432",
    "KQB09-20620": "KQB09_20620",
    "KQB09-23137": "KQB09_23137",
    "LCP85-384": "LCP85_384",
    "MT11-610": "MT11_610",
    "QA02-1009": "QA02_1009",
    "QA96-1749": "QA96_1749",
    "QBYN04-26041": "QBYN04_26041",
    "QC02-402": "QC02_402",
    "QN05-1460": "QN05_1460",
    "QN05-1509": "QN05_1509",
    "QN05-1743": "QN05_1743",
    "QN05-803": "QN05_803",
    "QS99-2014": "QS99_2014",
    "SP80-3280": "SP80_3280",
}

# Function to replace genotypes in the data
def replace_genotypes(data, mapping):
    for line in data:
        parts = line.split("\t")
        old_genotype = parts[2].strip().split('_')[0]  # Extracts the genotype part before the underscore
        if old_genotype in mapping:
            # Replace the old genotype with the correct one
            new_genotype = parts[2].replace(old_genotype, mapping[old_genotype])
            parts[2] = new_genotype
        yield "\t".join(parts)

# Read the content of the file panRNAomeClassificationTable_0.8.tsv
with open("panRNAomeClassificationTable_0.8.tsv", "r") as input_file:
    data = input_file.readlines()

# Apply genotype replacement
new_data = list(replace_genotypes(data, genotype_mapping))

# Write the modified data to the same file or a new file
with open("panRNAomeClassificationTable_0.8_fixedNames.tsv", "w") as output_file:
    output_file.writelines(new_data)
