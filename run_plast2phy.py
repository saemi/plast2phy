#!/usr/bin/env python

'''Script reads in the coding sequences of plastomes, and
generates a concatenead alignment nexus file for either GARLI or Mr.Bayes
'''

import argparse
import ConfigParser
import os
import shutil
import sys

from Bio import SeqIO
from glob import glob
from random import randint

def parse_arguments():
    '''Parses arguments from input'''
    parser = argparse.ArgumentParser(description='Converts input fasta files from ' +\
        'either DOGMA or NCBI plastid files, to plast2phy format.')
    parser.add_argument('-c', '--conf', help='Config file', required=True)
    return parser


def parse_config_file(config_file):
    '''Parses the config file'''
    # Global variables
    global INPUT_AA
    global INPUT_CDS
    global LIST_OF_SPECIES
    global GENES_TO_EXCLUDE
    global Config
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)
    INPUT_AA = ConfigSectionMap('folders')['input_folder'] + '/'
    INPUT_CDS = ConfigSectionMap('folders')['input_folder'] + '/'
    LIST_OF_SPECIES = get_list_of_species(ConfigSectionMap('species')['species_to_use'],
        int(ConfigSectionMap('species')['number_of_species']))
    GENES_TO_EXCLUDE = ConfigSectionMap('genes')['genes_to_exclude'].\
        replace(' ', '').split()
    return Config


def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print "exception on %s!" % option
            dict1[option] = None
    return dict1


def check_input_files(aa_dir, cds_dir):
    '''Checks the number of files and their formats, the script stops if a
    test fails'''

    check_cds_files(cds_dir)
    check_aa_files(aa_dir)
    print "# Files successfully pass quality check"


def get_list_of_species(species_string, no_species):
    '''reads in list of species'''
    species_list = []
    species_list = species_string.split()
    if len(species_list) == no_species:
        return species_list
    else:
        sys.exit("ERROR: Number of species and species listed doesn't match " +\
            "in the config file, please doubble check")


def check_cds_files(cds_dir):
    '''Function that does a stupid test if sequences are in nt format it is
    only based on the number of unique characters, maybe find something
    better later.
    '''

    list_of_cds = glob(cds_dir + '*')
    for fasta_file in list_of_cds:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = set(str(record.seq))
            if len(sequence) > 4:
                print 'Something wrong with %s' % record.id
                print str(record.id) + '\n' + str(record.seq)
                sys.exit('Check Sequnce, are they really in nucleotide format?')
    return


def check_aa_files(aa_dir):
    '''Function that does a stupid test if sequences are in aa format
    it is only based on the length of the sequence and number of unique
    characters, maybe find something better later.
    '''

    list_of_aa = glob(aa_dir + '*')
    for fasta_file in list_of_aa:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = set(str(record.seq))
            if len(sequence) < 4 and len(record) > 10:
                print 'Something wrong with %s' % record.id
                print str(record.id) + '\n' + str(record.seq)
                sys.exit('Check Sequnce, are they really in amino acid format?')


def create_list_of_genes(regions_to_use):
    '''Generates a set of genes in the cds fasta files, removes genes that
    are to be excluded. returns a sorted tuple'''

    set_of_cds = []
    list_of_genes = []
    if ['cds'] == regions_to_use:
        set_of_cds = get_list_of_cds()
    list_of_genes = set_of_cds
    list_of_genes.sort()
    set_of_genes = list(set(list_of_genes))
    return tuple(set_of_genes)


def get_list_of_cds():
    '''Creatse a list of cds to use'''
    list_of_cds = []
    list_of_all_cds = []
    for species in LIST_OF_SPECIES:
        list_of_all_cds.append(INPUT_CDS + '/' + species + '_cds.fsa')
    for fasta_file in list_of_all_cds:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if str(record.id) not in GENES_TO_EXCLUDE:
                list_of_cds.append(str(record.id))
    set_of_cds = list(set(list_of_cds))
    return set_of_cds


def clean_output_dict(outdir):
    '''Removes any existing file in directory'''

    list_of_files = glob(outdir + '*')
    if len(list_of_files) > 0:
        for i in list_of_files:
            os.remove(i)

def make_output_directory(directory):
    '''Make directories and handles errors'''
    try:
        os.mkdir(directory)
    except OSError:
        sys.exit("Output directory: " + directory + ' already exists' +
            " specify a different direcory")

# Functions that write outputs to disk (not helper functions)

def write_unaligned_genes(outdir_genes, list_of_genes, regions_to_use):
    '''Reads in the DOGMA outputs, changes gene names to species names and
    writes the output to a folder (i.e. input for MSA program)
    '''

    make_output_directory(outdir_genes)
    file_endings = {'cds': '_cds.fsa'}
    input_dir = INPUT_CDS
    file_ending = file_endings[regions_to_use[0]]
    write_unalniged_fasta_files(input_dir, list_of_genes, file_ending,
        outdir_genes)


def write_unalniged_fasta_files(input_dir, list_of_genes, file_ending,
        outdir_genes):
    '''Writes genes to directories'''

    species_list = LIST_OF_SPECIES
    for species in species_list:
        fasta_in = input_dir + species + file_ending
        in_dict = {}
        for contig in SeqIO.parse(fasta_in, 'fasta'):
            if contig.id in in_dict.keys():
                if contig.id in GENES_TO_EXCLUDE:
                    pass
                else:
                    sys.exit("Two instances of the gene " + contig.id +\
                        " it is not in " + \
                        "the list of genes to exclude, doubble check dogma " +\
                        "file, probably unjoined exons") % str(contig.id)
            else:
                in_dict[contig.id] = contig
        for gene in list_of_genes:
            if gene in in_dict.keys():
                gene_record = in_dict[gene]
                f_out = open(outdir_genes + gene +'.fa', 'aw')
                f_out.write('>' + species + '\n' + str(gene_record.seq) + '\n')
                f_out.close()
    return

def check_duplicate_gens(input_fasta):
    records = []
    for gene in SeqIO.parse(input_fasta, 'fasta'):
        records.append(gene.id)
    for gene_id in set(records):
        if records.count(gene_id) > 1:
            shutil.rmtree('output')
            sys.exit("More than one copy of %s, probably two exons "%gene_id +
            "to be joined")
        else:
            pass
    return


def run_mafft_alignment(genes_dir, outdir_align, list_of_genes):
    '''Generates mafft alignment on all the genes'''
    print "# Mafft alignment started"
    from sh import mafft
    make_output_directory(outdir_align)
    for gene in list_of_genes:
        mafft("--auto", genes_dir + gene + '.fa', _out= outdir_align + gene +
            '.aln.fa')
    print "# Mafft alignment complete"


def run_trimal(aln_dir, outdir_trim_align, list_of_genes, trim_option):
    '''Runs trimal'''
    from sh import trimal
    print "# TrimAl started"
    make_output_directory(outdir_trim_align)
    for gene in list_of_genes:
        trimal("-in", aln_dir + gene + '.aln.fa', "-out",
            outdir_trim_align + gene + '.trim.aln.fa', trim_option)
    print "# TrimAl finshed"
    return


def run_modeltest(trim_aln_dir, outidr_model_test, list_of_genes):
    '''Runs jModleTest2 on all alignments'''
    print "# Modeltesting started"
    from sh import run_jmodeltest
    make_output_directory(outidr_model_test)
    for gene in list_of_genes:
        run_jmodeltest(trim_aln_dir + gene + '.trim.aln.fa', outidr_model_test
        + gene + '.jmodeltest.out')
    print "# Modeltesting completed"
    return


def conv_fasta_to_phylip(trim_aln_dir, gene_trees, list_of_genes):
    '''Converst fasta alignments to phylip'''
    import sh
    make_output_directory(gene_trees)
    run = sh.Command("./bin/convertFasta2Phylip.sh")
    for gene in list_of_genes:
        run(trim_aln_dir + gene + '.trim.aln.fa',
            _out=gene_trees + gene + ".phy")
    return


def run_raxml(gene_trees, no_searches, model, outgroup, list_of_genes):
    '''Runs raxml on all alignments'''
    print "# RAXML started"
    print "# RAXML searching for best trees "
    from sh import raxmlHPC as raxml
    os.chdir(gene_trees)
    for gene in list_of_genes:
        raxml("-m", model, "-s", gene +".phy", "-n", gene + '.' + model, "-N",
            no_searches, "-o", outgroup, '-p', '12345')
    os.chdir('../../')
    print "# RAXML searching tree search finshed"
    return


def run_raxml_bstrap(gene_trees, no_bstrap, model, list_of_genes):
    '''Raxml boostrapping '''
    print "# RAXML bootstrapping"
    from sh import raxmlHPC as raxml
    os.chdir(gene_trees)
    for gene in list_of_genes:
        raxml("-m", model, "-p", '12345', "-b", '12345', "-#", no_bstrap,
        "-s", gene +".phy", "-n", gene + '.boot')
    os.chdir('../../')
    return


def raxml_consensus(gene_trees, model, outgroup, list_of_genes):
    '''Generates consensus trees from bootstraps and puts support values on
    best trees'''

    from sh import raxmlHPC as raxml
    os.chdir(gene_trees)
    for gene in list_of_genes:
        raxml("-m", model, "-p", "12345", "-f", "b", "-t", "RAxML_bestTree." +
            gene + '.' + model, "-z", "RAxML_bootstrap." + gene + '.boot', "-n",
            gene + ".cons.tre", "-o", outgroup)
    from sh import cat
    consensus_trees = glob('RAxML_bipartitions.*.cons.tre')
    cat(consensus_trees, _out='all_consensus_gene_trees.tre')
    # NEED TO ADD, ETE2 OR SOMETHING TO PUT NAMES ON TREES
    os.chdir('../../')
    print "# RAXML bootstrap search finished"
    return


def make_concat_nexus(trim_aln_dir, outdir_concat_nex, list_of_genes):
    '''Concatenates the alignments into a nexus file '''
    make_output_directory(outdir_concat_nex)
    nexus_header = generate_nexus_header(list_of_genes, trim_aln_dir)
    nexus_main = generate_nexus_main(trim_aln_dir, list_of_genes)
    nexus_end = generate_nexus_tail()
    nexus_file_wo_models = nexus_header + nexus_main + nexus_end
    write_nexus_file_wo_models(nexus_file_wo_models, outdir_concat_nex)
    return


def write_nexus_file_wo_models(nexus_file_wo_models, outdir_concat_nex):
    out_nexus_file = outdir_concat_nex + 'plast2phy_lathyrus_wo_models.nex'
    f_out = open(out_nexus_file, 'w')
    f_out.write(nexus_file_wo_models)
    f_out.close()
    print '# Nexus file (without models) written to ' + out_nexus_file
    return


# Helper functions for generating nexus file
def generate_nexus_header(list_of_genes, trim_aln_dir):
    ''' Generates this kind of nexus header
    #NEXUS
        begin data;
        dimensions ntax= 12 nchar=374705;
        format datatype=dna missing=? gap =- interleave;
    '''
    no_taxa = int(ConfigSectionMap('species')['number_of_species'])
    no_char = get_total_no_characters(trim_aln_dir, list_of_genes)
    nexus_header = ['#NEXUS', '\n',
        '\t', 'begin data;', '\n',
        '\t', 'dimensions ntax = ', str(no_taxa), ' nchar = ', str(no_char),
        ';', '\n', '\t', 'format datatype = dna missing = ? gap = - interleave;',
        '\n', 'matrix', '\n']
    return ''.join(nexus_header)


def get_total_no_characters(trim_aln_dir, list_of_genes):
    '''Counts the number of characters in all the alignments'''
    no_chars = []
    for gene in list_of_genes:
        f_in = open(trim_aln_dir + gene + '.trim.aln.fa')
        no_char = int(f_in.readline().split()[1])
        no_chars.append(no_char)
        f_in.close()
    return sum(no_chars)


def generate_nexus_main(trim_aln_dir, list_of_genes):
    '''Concatenates the fasta alignment files into a single nexus file, line
    length is set to 80'''
    # no_white_spaces = get_longest_sp_name(LIST_OF_SPECIES) + 1
    line_length = 80
    no_white_spaces = len(max(LIST_OF_SPECIES, key=len)) + 1
    equal_sp_name = add_white_spaces_sp_name(LIST_OF_SPECIES, no_white_spaces)
    concat_alignments = concat_char(trim_aln_dir, list_of_genes,
        LIST_OF_SPECIES)
    nexus_lines = generate_nexus_lines(concat_alignments,
        LIST_OF_SPECIES, equal_sp_name, line_length, no_white_spaces)
    return nexus_lines


def concat_char(trim_aln_dir, list_of_genes, list_of_species):
    '''Returns a dictionary with the species names and concatenated string of
    characters'''

    concat_characters = []
    concat_alingment = []
    for species in list_of_species:
        characters = []
        for gene in list_of_genes:
            fasta_file = trim_aln_dir + gene + '.trim.aln.fa'
            record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
            characters.append(str(record_dict[species].seq).upper())
        concat_characters.append(''.join(characters))
    concat_alingment = dict(zip(list_of_species, concat_characters))
    return concat_alingment


def add_white_spaces_sp_name(list_of_species, no_white_spaces):
    '''Returns a dictionary with species names where the correct no of white
    spaces have been added for the nexus file'''
    old_species_names = list_of_species
    new_species_names = []
    for species in list_of_species:
        new_species_name = ''
        diff = no_white_spaces - len(species)
        new_species_name = species + (diff * ' ')
        new_species_names.append(new_species_name)
    return dict(zip(old_species_names, new_species_names))


def get_no_nexus_lines(concat_alignments, no_white_spaces, line_length):
    '''Estimates how many lines the nexus file will be'''
    species = concat_alignments.keys()
    no_characters = len(concat_alignments[species[0]])
    no_char_per_line = line_length - no_white_spaces
    if no_characters % no_char_per_line == 0:
        return no_characters/no_char_per_line
    else:
        return no_characters/no_char_per_line + 1


def generate_nexus_lines(concat_alignments, list_of_species, equal_sp_name,
    line_length, no_white_spaces):
    '''Generate the main nexus block'''
    no_nexus_lines = get_no_nexus_lines(concat_alignments, no_white_spaces,
        line_length)
    main_nexus_block = []
    no_char_per_line = line_length - no_white_spaces
    for i in range(0, no_nexus_lines):
        nexus_block = []
        for species in list_of_species:
            if i == 0:
                char_start = i * no_char_per_line
                char_end = char_start + no_char_per_line
            elif i == no_nexus_lines - 1:
                char_start = i * no_char_per_line
                char_end = char_start + len(concat_alignments[species]\
                    [char_start:])
            else:
                char_start = i * no_char_per_line + 1
                char_end = char_start + no_char_per_line
            species_line = []
            species_line.append(equal_sp_name[species])
            species_line.append(concat_alignments[species][char_start:char_end])
            species_line.append('\n')
            nexus_block.append(''.join(species_line))
        main_nexus_block.append(''.join(nexus_block) + '\n')
    return ''.join(main_nexus_block)


def generate_nexus_tail():
    '''Genrats the end of the nexus block'''
    nexus_end = [';', '\n', 'end;', '\n']
    return ''.join(nexus_end)


def main():
    parser = parse_arguments()
    args = parser.parse_args()
    parse_config_file(args.conf)
    output_directory = ConfigSectionMap('folders')['output_folder'] + '/'
    print
    print '###################################################################'
    print '# Plast2Phy started'
    print '###################################################################'
    check_input_files(INPUT_AA, INPUT_CDS)
    print '###################################################################'
    print
    '''Reads in dogma out cds fils and converts them to a format that aligner
    can read (e.g. the same gene from all speceies together in file)'''
    make_output_directory(output_directory)
    #
    regions_to_use = ['cds']
    list_of_genes = create_list_of_genes(regions_to_use)
    # print(list_of_genes)
    # Output directories
    outdir_genes = output_directory + '1_genes_unaligned/'
    outdir_align = output_directory + '2_gene_alignments/'
    outdir_trim_align = output_directory + '3_trimmed_alignments/'
    # outidr_model_test = output_directory +'4_model_testing/'
    outdir_gene_trees = output_directory + '4_gene_trees/'
    outdir_concat_nex = output_directory + '5_concat_nexus/'
    # 1. Write genes
    write_unaligned_genes(outdir_genes, list_of_genes, regions_to_use)
    # 2. Generate alignments
    run_mafft_alignment(outdir_genes, outdir_align, list_of_genes)
    # 3. Trim alignments
    trim_option = ConfigSectionMap('alignment')['trimming']
    run_trimal(outdir_align, outdir_trim_align, list_of_genes, trim_option)
    # 4. Runs jModelTest
    # model_test = 'FALSE' # NEED TO ADD CONF PARSER
    # if model_test == 'TRUE':
        # run_modeltest(outdir_trim_align, outidr_model_test, list_of_genes)
    # 4.1 Converts fasta to phylip and shortens species name
    conv_fasta_to_phylip(outdir_trim_align, outdir_gene_trees, list_of_genes)
    # 4.2 Runs raxml search of best tree
    no_searches = ConfigSectionMap('gene_tree_search')['number_of_searches']
    model = ConfigSectionMap('gene_tree_search')['model']
    outgroup = ConfigSectionMap('species')['outgroup']
    run_raxml(outdir_gene_trees, no_searches, model, outgroup, list_of_genes)
    # 5.3 Runs raxmlt bootstrap
    raxml_bstrap = ConfigSectionMap('gene_tree_search')['bootstrap']
    if raxml_bstrap == 'yes':
        no_bstrap = ConfigSectionMap('gene_tree_search')['number_of_bootstraps']
        run_raxml_bstrap(outdir_gene_trees, no_bstrap, model, list_of_genes)
        # 5.4 Makes raxml consensus trees
        raxml_consensus(outdir_gene_trees, model, outgroup, list_of_genes)
        # NEED TO ADD SOMETHING TO CONCATENATE THE TREES AND KEEP NAMES
    # 6.1 Generate concatenated nexus file
    make_concat_nexus(outdir_trim_align, outdir_concat_nex, list_of_genes)


main()
