# InfantMicrobiotaQC

There are 3 independent data sets for this project. The initial processing steps were identical and are listed below. The subsequent processing steps were completed in R with the code available here. Note that these steps begin with OTU picking with QIIME. 


1. OTU picking with QIIME pick_open_reference_otus.py -i InputFileName.fa -o OutputFolderName
2. The output will be an OTU table and a phylogenetic tree.  The file name of the tree is rep_set.tre
The file name of the OTU table is otu_table_mc2.biom or otu_table_mc2_w_tax.biom with the taxonomy table.
3. Convert the .biom to .json format. biom convert -i otu_table_mc2_w_tax.biom -o jsonOTUtax.biom --to-json
4. Read the OTU/taxonomy tables and phylogenetic tree into R for subsequent analyses. 

See R code 
