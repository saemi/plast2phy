plast2phy
=========

Pipeline written in python, that extracts protein coding genes from whole plastome assemblies/annotation and uses them in phylogenetic reconstruction. Is is written and tested on Ubuntu GNU/Linux, V.14.04. Should work on previous versions of Ubuntu. 

Dependenceis.

Python modules
* sh (http://amoffat.github.io/sh/)
* Biopython (python-biopython in ubuntu)

Programs/scripts in PATH
* mafft (apt-get install mafft in ubuntu) [http://mafft.cbrc.jp/alignment/software/linux.html]
* trimal v.1.2 () [http://trimal.cgenomics.org/]
* raxml v.8.x.x() [https://github.com/stamatak/standard-RAxML]
* convertFasta2Phylip.sh (located in ./bin/), copy to PATH


How to run pipeline:

1) run ./configure_plast2phy.py. This checks if you have all the required programs and python modules installed

2) run ./format_input_plast2phy.py for each plastome you want to use in the phylogenetic reconstruction

3) edit the plast2phy.conf, with the location and name of input files and vairous settings.

4) run ./run_plast2phy.py, this reads in the conf file and outputs a concatenated NEXUS matrix. You can also set it to output individual gene trees and so on.

If you have any questions, please email saemundur dot sveinsson (at) gmail dot com. No one except me has used it so far, so it is not unlikely that if you are the first to try it, that it won't work for you. But please let me know
