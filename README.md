## General info
This folder contains selected R scripts used to analyze metagenomic data obtained from dilbit-coated clay beads incubated at Douglas Channel, BC, Canada.
	
## Scripts and files
* getConsensus_v2.R - R script to calculate a consensus gene taxonomy
* example_DIAMOND_file.txt - example DIAMOND blastp file to be used with the getConsensus_v2.R script
* compostional_metrics.R - R script containing modified functions of the "propr" package for calculating compostional metrics (clr, alr and phi) from OTU/ASV count tables
	
## Usage
### getConsensus_v2.R
1) Obtain the MetaErg reference databases
```
get http://ebg.ucalgary.ca/metaerg/db.tar.gz -P $HOME
tar -xvzf $HOME/db.tar.gz
```
2) Run a DIAMOND blastp analysis with the genes to by classified and the DIAMOND database provide by MetaErg
```
diamond blastp -k 0 --quiet -p 4 -q genes_for_classification.faa -d genomedb.dmnd -e 0.00001 -f 6 qseqid sseqid qlen qstart qend sstart send qframe pident bitscore evalue qcovhsp -o genes_for_classification.DIAMOND_out.txt
```
3) Convert the DIAMOND results file into the more common tsv format using sed
```
sed -i '' -e 's/;/     /g' genes_for_classification.DIAMOND_out.txt
```
4) Follow the getConsensus_v2.R script interactively in R Studio

## References
* Buchfink B, Xie C & Huson DH (2015) Fast and sensitive protein alignment using DIAMOND. Nature Methods 12: 59-60.
* Dong X & Strous M (2019) An integrated pipeline for annotation and visualization of metagenomic contigs. Frontiers in Genetics 10.
* Quinn TP, Richardson MF, Lovell D & Crowley TM (2017) propr: An R-package for Identifying Proportionally Abundant Features Using Compositional Data Analysis. Scientific Reports 7: 16252.
