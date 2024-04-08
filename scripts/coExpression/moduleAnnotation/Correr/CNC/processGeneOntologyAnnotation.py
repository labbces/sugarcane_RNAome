#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Script to process GO universe annotation, panRNAome genes and generate GO universe filtered by ontology for enrichment analysis.")
parser.add_argument("-a", "--annotation", dest="annotation_file", metavar="GO_universe_annotation_list", required=True, help="Path to the GO annotation list")
parser.add_argument("-p", "--panrnaome",  dest="panrnaome_file", metavar="panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv", required=True, help="Path to the panRNAome file")
parser.add_argument("-o", "--output", dest="output_file", metavar="GO_annotations", required=True, help="Prefix to output file")
parser.add_argument("-og", "--ontology", dest="ontology", metavar="BP, CC or MF", required=True, help="filter by ontology")
parser.add_argument("-ppv", dest="ppv_threshold", type=float, metavar="0.7", required=True, help="threshold to filter PPV scores [0.00-1.00]")
args = parser.parse_args()

def process_annotation_file(annotation_file):
    transcript_annotations = {}
    with open(annotation_file, 'r') as f:
        for line in f:
            
            if line.startswith("qpid"):
                continue  # Skip header line
            
            parts = line.split("\t")
            ontology = parts[1]
            if ontology != f"{args.ontology}":
                continue  # Skip lines where the ontology is not "BP"
            
            #print(parts[5])
            ppv = float(parts[5])  # Extract the value from the ARGOT_PPV column
            if ppv <= args.ppv_threshold:
                continue  # Skip lines with PPV value less than or equal to 0.9
            
            transcript = parts[0]
            annotation = parts[2]
            if transcript in transcript_annotations:
                transcript_annotations[transcript].add(annotation)
            else:
                transcript_annotations[transcript] = {annotation}
    return transcript_annotations

out = args.output_file + "_" + args.ontology + "_" + "PPV" + str(args.ppv_threshold) + ".tsv"

def process_transcript_file(transcript_file, transcript_annotations):
    with open(out, 'w') as output:
        with open(transcript_file, 'r') as f:
            for line in f:
                parts = line.split("\t")
                gene = parts[1]
                transcript = parts[2]
                if transcript in transcript_annotations:
                    annotations = transcript_annotations[transcript]
                    annotations_with_go = [f"GO:{a}" for a in annotations]
                    annotation_str = ','.join(annotations_with_go)
                    #print(f"{gene}\t{annotation_str}")
                    output.write(f"{gene}\t{annotation_str}\n")

#annotation_file = 'scripts/coExpression/moduleAnnotation/Correr/CNC/1000GO_universe_annotation_list2'
#panrnaome_file = 'scripts/coExpression/moduleAnnotation/Correr/CNC/1000panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv'
#output_file = 'scripts/coExpression/moduleAnnotation/Correr/CNC/GO_annotations_BP_PPV0.7.tsv'

transcript_annotations = process_annotation_file(args.annotation_file)
process_transcript_file(args.panrnaome_file, transcript_annotations)
print(f"output saved as {out}")