plast2phy

Pipeline written in python, that extracts protein coding genes from whole plastome assemblies/annotation and uses them in phylogenetic reconstruction. Is is written and tested on Ubuntu GNU/Linux, V.14.04. Should work on previous versions of Ubuntu. 

Dependenceis.

Python modules
* sh (http://amoffat.github.io/sh/) [apt-get install python-pip && sudo pip install sh in ubuntu]
* Biopython [apt-get python-biopython in ubuntu]

Programs/scripts in PATH. Copies are located in the requried_programs folder 
* mafft [http://mafft.cbrc.jp/alignment/software/linux.html]
* trimal v.1.2 [http://trimal.cgenomics.org/]
* raxml v.8.x.x [https://github.com/stamatak/standard-RAxML]
* convertFasta2Phylip.sh, copy to PATH

How to install pipeline:

1) install the Python modules and required programs, listed above. They need to be in the PATH so a simple ln -s from where the programs are compiled to ~/bin/ should work for the programs. Unfortuanetly you might need root accsess to install the Python modules.
 
How to test the pipeline:

1) run ./configure_plast2phy.py. This checks if you have all the required programs and python modules installed

2) cd test_data/

3) mkdir input && ../bin/format_input_plast2phy.py -i plastid_cds_aa/Pisum_sativum_plastidGenome_DOGMA_nucleotide.fa -o input/ -s Pisum_sativum -f d -t n && ../bin/format_input_plast2phy.py -i plastid_cds_aa/Trifolium_glanduliferum_plastidGenome_DOGMA_nucleotide.fa -o input/ -s Trifolium_glanduliferum -f d -t n && ../bin/format_input_plast2phy.py -i plastid_cds_aa/Trifolium_strictum_plastidGenome_DOGMA_nucleotide.fa -s Trifolium_strictum -f d -t n -o input/ && ../bin/format_input_plast2phy.py -i plastid_cds_aa/Trifolium_lupinaster_plastidGenome_DOGMA_nucleotide.fa -o input/ -s Trifolium_lupinaster -f d -t n

4) ../run_plast2phy.py -c plast2phy_test.conf

The pipeline should run without errors and end by printing out: Nexus file (without models) written to output/5_concat_nexus/t2phy_test.nex 

How to run the pipeline:

1) run ./format_input_plast2phy.py for each plastome you want to use in the phylogenetic reconstruction

2) edit the plast2phy.conf, with the location and name of input files and vairous settings.

3) run ./run_plast2phy.py, this reads in the conf file and outputs a concatenated NEXUS matrix. You can also set it to output individual gene trees and so on.

If you have any questions, please email saemundur dot sveinsson (at) gmail dot com. No one except me has used it so far, so it is not unlikely that if you are the first to try it, that it won't work for you. But please let me know
