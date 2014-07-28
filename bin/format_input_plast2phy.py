#!/usr/bin/env python

import argparse
import sys

from Bio import SeqIO
from os import path


""" This script pepares input files for the run_plast2phy.py script, by parsing
either DOGMA or NCBI formatted plastid genes.
* DOGMA formatted files are extracted from DOGMA, through the 'Extract Sequeces' 
menu, select either protein coding aa or nt sequences

* NCBI formatted files are extracted from NCBI's Genbank by navigating to the 
plastid nucleotide site and selecting send -> coding sequences -> and select either 
FASTA Protein or FASTA Nucleotide.
"""

def parse_arguments():
    '''Parses arguments from input'''
    parser = argparse.ArgumentParser(description='Converts input fasta files from ' +\
        'either DOGMA or NCBI plastid files, to plast2phy format.')
    parser.add_argument('-i', '--in_fasta', help='Input Fasta file', 
        required=True)
    parser.add_argument('-o', '--out_folder', help='Output folder, usually ' +\
        'input/', required=True)
    parser.add_argument('-s', '--species', help='Name of species (or individual)' +\
        ', seperate words using underscore(_), e.g. Pisum_sativum_cds.fa',required=True)
    parser.add_argument('-f', '--in_format', help='DOGMA (d) or NCBI (n) format of ' +\
        'input file', choices=['d','n'], required=True)
    parser.add_argument('-t', '--type', help='Nucleotide (n) or protein (p) ' +\
        ' sequences in inputfile input file', choices=['n','p'], required=True)
    return(parser)


def read_fasta(fastafile, input_format):
    '''Reads in the input fasta files, renames genes and outputs a tuple with
    gene names and gene sequences.
    '''

    new_gene_names = []
    gene_sequences = []
    if input_format == 'd':
        for gene in SeqIO.parse(fastafile, 'fasta'):
            new_gene_names.append(rename_dogma_gene(gene))
            gene_sequences.append(str(gene.seq))
        if new_gene_names.count('petD') > 1:
            # One problem with DOGMA is that it is sometimes unable to join the
            # 2 petD exons
            fixed_genes = fix_petD(new_gene_names,gene_sequences)
        if new_gene_names.count('rps12') > 1:
            # Anoher problem with DOGMA that sometimes the rps12 exons are
            # translated in opposite directions
            x = 0
            for gene in new_gene_names:
                if gene == 'rps12':
                    new_gene_names.pop(x)
                    gene_sequences.pop(x)
                x += 1
            rps12_contigs = []
            for gene in SeqIO.parse(fastafile, 'fasta'):
                if rename_dogma_gene(gene) == 'rps12':
                    rps12_contigs.append(gene)
            rps12_joined = join_rps12_contigs(rps12_contigs)
            new_gene_names.append(rps12_joined[0])
            gene_sequences.append(rps12_joined[1])
            return(zip(new_gene_names, gene_sequences))
        else:
            return(zip(new_gene_names, gene_sequences))
    elif input_format == 'n':
        for gene in SeqIO.parse(fastafile, 'fasta'):
            new_gene_names.append(rename_ncbi_gene(gene))
            gene_sequences.append(str(gene.seq))
        return(zip(new_gene_names,gene_sequences))


def rename_dogma_gene(gene):
    '''Renames the DOGMA formatted genes'''
    new_gene_name = gene.description.split()[0].strip(':')
    return(new_gene_name)


def fix_petD(new_gene_names,gene_sequences):
    '''This function removes the small exon of petD in DOGMA formats'''
    x = 0
    index = 0
    for gene in zip(new_gene_names,gene_sequences):
        if gene[0] == 'petD':
            if len(gene[1]) == 9:
                index = x
            elif len(gene[1]) == 3:
                index = x
        x+=1
    # print(index)
    new_gene_names.pop(index)
    gene_sequences.pop(index)
    return(zip(new_gene_names,gene_sequences))


def rename_ncbi_gene(gene):
    '''Renames the NCBI formatted genes'''
    new_gene_name = gene.description.split()[1].split('=')[1].strip(']')
    return(new_gene_name)


def join_rps12_contigs(rps12_contigs):
    joined_rps12_gene = []
    if str(rps12_contigs[0].seq)[0:3] == 'ATG':
        joined_rps12_gene.append(str(rps12_contigs[0].seq))
        joined_rps12_gene.append(str(rps12_contigs[1].seq))
    elif str(rps12_contigs[0].seq)[0] == 'M':
        joined_rps12_gene.append(str(rps12_contigs[0].seq))
        joined_rps12_gene.append(str(rps12_contigs[1].seq))
    else:
        joined_rps12_gene.append(str(rps12_contigs[1].seq))
        joined_rps12_gene.append(str(rps12_contigs[0].seq))
    rps12_joined = ['rps12',''.join(joined_rps12_gene)]
    return(rps12_joined)


def write_output(genes, species, seq_type, out_folder):
    '''Writes the output file'''
    name_seq_type = {}
    name_seq_type['n'] = 'cds'
    name_seq_type['p'] = 'aa'
    out_fasta_file = out_folder +  '/' + species + '_' + name_seq_type[seq_type] + '.fsa'
    if path.isdir(out_folder) == False:
        os.mkdir(out_folder)
    if path.isfile(out_fasta_file) == True:
        sys.exit('ERROR: %s is invalid output file name, already a file with that name'\
            % out_fasta_file)
    f = open(out_fasta_file, 'w')
    for gene in genes:
        f.write('>' + gene[0] + '\n' + gene[1] + '\n')
    f.close()
    no_genes = str(len(genes))
    gene_names = []
    for gene in genes:
        gene_names.append(gene[0])
    for gene in set(gene_names):
        if gene_names.count(gene) > 1:
            print "WARNING: %s has more than one copies, probably unjoined exons."%gene
    print no_genes + ' written to %s' % out_fasta_file
    return


def main():
    '''Main function, runs all required functions'''
    parser = parse_arguments()
    args = parser.parse_args()
    gene_sequences = read_fasta(args.in_fasta, args.in_format)
    write_output(gene_sequences, args.species, args.type, args.out_folder)


main()
