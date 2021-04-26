#Version 2021-04-26
library("dplyr")
library("RSQLite")
library("readr")

getConsensus <- function(diamond.table,
                         top.percent = 5){
  
  #Calculate bitscore cutoff
  max.bitscore <- max(diamond.table$bitscore)
  bitscore.cutoff <- max.bitscore * (1 - (top.percent / 100))
  
  diamond.table.ordered <-  diamond.table[diamond.table$bitscore >= bitscore.cutoff,]
  
  #Take apart the lineage for consensus building
  diamond.table.ordered$kingdom <- sapply(strsplit(diamond.table.ordered$lineage, split = ";"), function(x) x[1])
  diamond.table.ordered$kingdom <- substr(diamond.table.ordered$kingdom, start = 4, stop = nchar(diamond.table.ordered$kingdom))
  
  diamond.table.ordered$phylum <- sapply(strsplit(diamond.table.ordered$lineage, split = ";"), function(x) x[2])
  diamond.table.ordered$phylum <- substr(diamond.table.ordered$phylum, start = 4, stop = nchar(diamond.table.ordered$phylum))
  
  diamond.table.ordered$class <- sapply(strsplit(diamond.table.ordered$lineage, split = ";"), function(x) x[3])
  diamond.table.ordered$class <- substr(diamond.table.ordered$class, start = 4, stop = nchar(diamond.table.ordered$class))
  
  diamond.table.ordered$order <- sapply(strsplit(diamond.table.ordered$lineage, split = ";"), function(x) x[4])
  diamond.table.ordered$order <- substr(diamond.table.ordered$order, start = 4, stop = nchar(diamond.table.ordered$order))
  
  diamond.table.ordered$family <- sapply(strsplit(diamond.table.ordered$lineage, split = ";"), function(x) x[5])
  diamond.table.ordered$family <- substr(diamond.table.ordered$family, start = 4, stop = nchar(diamond.table.ordered$family))
  
  diamond.table.ordered$genus <- sapply(strsplit(diamond.table.ordered$lineage, split = ";"), function(x) x[6])
  diamond.table.ordered$genus <- substr(diamond.table.ordered$genus, start = 4, stop = nchar(diamond.table.ordered$genus))
  
  #Do the consensus
  output.df <- data.frame("gene_id" = unique(diamond.table$gene_id)[1])
  
  grouped <- group_by(.data = diamond.table.ordered, by = kingdom)
  summed <- tally(grouped)
  summed.ordered <- summed[order(-summed$n),]
  if(length(summed.ordered$by) == 1){
    output.df$kingdom <- summed.ordered$by[1]
  } else {output.df$kingdom <- "NULL"}
  
  grouped <- group_by(.data = diamond.table.ordered, by = phylum)
  summed <- tally(grouped)
  summed.ordered <- summed[order(-summed$n),]
  if(length(summed.ordered$by) == 1){
    output.df$phylum <- summed.ordered$by[1]
  } else {output.df$phylum <- "NULL"}
  
  grouped <- group_by(.data = diamond.table.ordered, by = class)
  summed <- tally(grouped)
  summed.ordered <- summed[order(-summed$n),]
  if(length(summed.ordered$by) == 1){
    output.df$class <- summed.ordered$by[1]
  } else {output.df$class <- "NULL"}
  
  grouped <- group_by(.data = diamond.table.ordered, by = order)
  summed <- tally(grouped)
  summed.ordered <- summed[order(-summed$n),]
  if(length(summed.ordered$by) == 1){
    output.df$order <- summed.ordered$by[1]
  } else {output.df$order <- "NULL"}
  
  grouped <- group_by(.data = diamond.table.ordered, by = family)
  summed <- tally(grouped)
  summed.ordered <- summed[order(-summed$n),]
  if(length(summed.ordered$by) == 1){
    output.df$family <- summed.ordered$by[1]
  } else {output.df$family <- "NULL"}
  
  grouped <- group_by(.data = diamond.table.ordered, by = genus)
  summed <- tally(grouped)
  summed.ordered <- summed[order(-summed$n),]
  if(length(summed.ordered$by) == 1){
    output.df$genus <- summed.ordered$by[1]
  } else {output.df$genus <- "NULL"}
  
  output.df$support <- length(diamond.table.ordered$gene_id)
  
  colnames(output.df) <- c("gene_id","tax_kingdom","tax_phylum","tax_class","tax_order","tax_family","tax_genus","support")
  
  #Return consensus lineage
  return(output.df)
}


#Carry out DIAMOND search:
#diamond blastp -k 0 --quiet -p 4 -q genes_for_reannotation_20191210.faa -d /Applications/BioInformatics/MetaErg/databases/db/diamond/genomedb.dmnd -e 0.00001 -f 6 qseqid sseqid qlen qstart qend sstart send qframe pident bitscore evalue qcovhsp -o genes_for_reannotation_20191210.DIAMOND_out.txt
#Replace the "|" with a TAB in the DIAMOND results file
#sed -i '' -e 's/;/     /g' allDegGenes_genes4R_.DIAMOND_reannotated_20210114.txt

#Load the DIAMOND results file
annotation.diamond <- read_tsv(file = "../example_DIAMOND_file.txt", col_names = FALSE)

colnames(annotation.diamond) <- c("gene_id",
                                  "tax_id",
                                  "tax_gene",
                                  "qlen",
                                  "qstart",
                                  "qend",
                                  "sstart",
                                  "send",
                                  "qframe",
                                  "pident",
                                  "bitscore",
                                  "evalue",
                                  "qcovhsp")
#Load the metaerg database                               
con = dbConnect(RSQLite::SQLite(), dbname="metaerg.db")

#Get genome2taxon table as data frame
annotation.taxonomy = dbGetQuery( con,'select * from genome2taxon' )

#Merge link metaerg gene names with their taxonomy
annotation.diamond <- merge(annotation.diamond,annotation.taxonomy, by.x = "tax_id", by.y = "gid", all.x = TRUE)

#We will only need the gene name, taxonomy and biscore
annotations <- annotation.diamond[,c("gene_id","lineage","bitscore")]


#Calculate the consensus taxonomy
growing.df <- data.frame("gene_id" = "DUMMY",
                         "tax_kingdom" = "NULL",
                         "tax_phylum" = "NULL",
                         "tax_class" = "NULL",
                         "tax_order" = "NULL",
                         "tax_family" = "NULL",
                         "tax_genus" = "NULL",
                         "support" = 10)

for(i in 1:length(unique(annotations$gene_id))){
  gene.id <- unique(annotations$gene_id)[i]
  temp.df <- annotations[annotations$gene_id == gene.id, c("gene_id","bitscore","lineage")]
  getConsensus(diamond.table = temp.df, top.percent = 5)
  growing.df <- rbind.data.frame(growing.df,getConsensus(diamond.table = temp.df, top.percent = 5))
}

growing.df$gene_id <- as.character(growing.df$gene_id)
growing.df$tax_kingdom <- as.character(growing.df$tax_kingdom)
growing.df$tax_phylum <- as.character(growing.df$tax_phylum)
growing.df$tax_class <- as.character(growing.df$tax_class)
growing.df$tax_order <- as.character(growing.df$tax_order)
growing.df$tax_family <- as.character(growing.df$tax_family)
growing.df$tax_genus <- as.character(growing.df$tax_genus)

final.annotations <- growing.df[!(growing.df$gene_id == "DUMMY"),]
