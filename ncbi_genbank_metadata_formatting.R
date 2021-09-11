# Mandy Waters
# July 2021
# Script to read in gisaid and SRA data and write out genbank metadata
# Usage (short) = ./ncbi_genbank_metadata_formatting.R -p /path/to/working/dir -t sequencingtype
# Usage (long) = ./ncbi_genbank_metadata_formatting.R -path /path/to/working/dir -seq_type sequencingtype
# Usage notes: only strings with spaces need quotes.  It will throw an error if you don't supply the path and sequence type name or if the sequence type is misspelled.
# Usage example = ./ncbi_genbank_metadata_formatting.R -p /mandy_waters/metadata -t "Oxford Nanopore GridION"

# load packages, but don't print the messages to terminal
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(tidyr))

# set argument / options
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
# vignette: http://www.icesi.edu.co/CRAN/web/packages/optparse/vignettes/optparse.pdf

option_list = list(
  make_option(c("-p", "--path"), action="store", default=NA, type='character',
              help="path to set working directory"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Should the program print extra stuff out? [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Make the program not be verbose."),
  make_option(c("-t", "--seq_type"), action="store", default=NA, type='character',
              help="sequencing type as a string in quotes")
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$v) {
  # print the variables to the terminal
  cat("\nworking directory:\n")
  cat(opt$path)
  cat("\n\nsequencing type:\n")
  cat(opt$seq_type)
  cat("\n\n\n\n")
}

# if all options are filled in then fill in the metadata otherwise error out
if(!is.na(opt$path) & !is.na(opt$seq_type)){
  if(opt$seq_type == "Illumina MiSeq" | opt$seq_type == "Illumina NextSeq" | opt$seq_type == "Oxford Nanopore GridION"){

  # set working directory
  setwd(opt$path)
  date = Sys.Date()
    
    if(opt$seq_type == "Oxford Nanopore GridION"){
    # Reading in ncbi sra and GISAID data
    ncbi_in = ldply(list.files(pattern = "metadata-.*-processed-ok.tsv"), read.delim, header=TRUE)
    ncbi_split = ncbi_in %>% separate(filename, c("accession_short",NA,NA), sep = "[.]", remove = FALSE)
    gisaid_in = ldply(list.files(pattern = "gisaid_hcov-19_.*.tsv"), read.delim, header = TRUE)
    gisaid_split = gisaid_in %>% separate(Virus.name, c(NA,NA,"accession_short", "year"), sep = "/", remove = FALSE)
    gisaid_filtered = gisaid_split[gisaid_split$accession_short %in% ncbi_split$accession_short,]
    ncbi_gisaid_merge = merge(x=ncbi_split, y=gisaid_filtered[,1:5], by = "accession_short", all.x=TRUE)
    
    # Creating an empty matrix to fill in with metadata 
    genbank_metadata = as.data.frame(matrix("",ncol=9,nrow=nrow(ncbi_gisaid_merge)))
    colnames(genbank_metadata) = c("Sequence_ID","country","host","isolate","collection-date","isolation-source","BioSample","SRA","note")
    
    # Filling in columns whose values do not change or will be put in with options
    genbank_metadata$Sequence_ID = ncbi_gisaid_merge$Virus.name
    genbank_metadata$country = "USA"
    genbank_metadata$host = "Homo sapiens"
    genbank_metadata$isolate = ncbi_gisaid_merge$accession_short
    genbank_metadata$`collection-date`= ncbi_gisaid_merge$Collection.date
    genbank_metadata$`isolation-source` = "patient isolate"
    genbank_metadata$BioSample = ncbi_gisaid_merge$biosample_accession
    genbank_metadata$SRA = ncbi_gisaid_merge$accession
    genbank_metadata$note = paste("GISAID_ID:", ncbi_gisaid_merge$Accession.ID, sep = " ")
    
    }else if(opt$seq_type == "Illumina MiSeq"){
      # Reading in ncbi sra and GISAID data
      ncbi_in = ldply(list.files(pattern = "metadata-.*-processed-ok.tsv"), read.delim, header=TRUE)
      ncbi_split = ncbi_in %>% separate(filename, c("accession_cdphe",NA), sep = "__", remove = FALSE)
      ncbi_split$accession_short = paste("CO-CDPHE-",ncbi_split$accession_cdphe,sep="")
      gisaid_in = ldply(list.files(pattern = "gisaid_hcov-19_.*.tsv"), read.delim, header = TRUE)
      gisaid_split = gisaid_in %>% separate(Virus.name, c(NA,NA,"accession_short", "year"), sep = "/", remove = FALSE)
      gisaid_filtered = gisaid_split[gisaid_split$accession_short %in% ncbi_split$accession_short,]
      ncbi_gisaid_merge = merge(x=ncbi_split, y=gisaid_filtered[,1:5], by = "accession_short", all.x=TRUE)
      
      # Creating an empty matrix to fill in with metadata 
      genbank_metadata = as.data.frame(matrix("",ncol=9,nrow=nrow(ncbi_gisaid_merge)))
      colnames(genbank_metadata) = c("Sequence_ID","country","host","isolate","collection-date","isolation-source","BioSample","SRA","note")
      
      # Filling in columns whose values do not change or will be put in with options
      genbank_metadata$Sequence_ID = ncbi_gisaid_merge$Virus.name
      genbank_metadata$country = "USA"
      genbank_metadata$host = "Homo sapiens"
      genbank_metadata$isolate = ncbi_gisaid_merge$accession_short
      genbank_metadata$`collection-date`= ncbi_gisaid_merge$Collection.date
      genbank_metadata$`isolation-source` = "patient isolate"
      genbank_metadata$BioSample = ncbi_gisaid_merge$biosample_accession
      genbank_metadata$SRA = ncbi_gisaid_merge$accession
      genbank_metadata$note = paste("GISAID_ID:", ncbi_gisaid_merge$Accession.ID, sep = " ")
    }else if(opt$seq_type == "Illumina NextSeq"){
      # Reading in ncbi sra and GISAID data
      ncbi_in = ldply(list.files(pattern = "metadata-.*-processed-ok.tsv"), read.delim, header=TRUE)
      ncbi_split = ncbi_in %>% separate(filename, c("accession_cdphe",NA,NA), sep = "[.]", remove = FALSE)
      ncbi_split$accession_short = paste("CO-CDPHE-",ncbi_split$accession_cdphe,sep="")
      gisaid_in = ldply(list.files(pattern = "gisaid_hcov-19_.*.tsv"), read.delim, header = TRUE)
      gisaid_split = gisaid_in %>% separate(Virus.name, c(NA,NA,"accession_short", "year"), sep = "/", remove = FALSE)
      gisaid_filtered = gisaid_split[gisaid_split$accession_short %in% ncbi_split$accession_short,]
      ncbi_gisaid_merge = merge(x=ncbi_split, y=gisaid_filtered[,1:5], by = "accession_short", all.x=TRUE)
      
      # Creating an empty matrix to fill in with metadata 
      genbank_metadata = as.data.frame(matrix("",ncol=9,nrow=nrow(ncbi_gisaid_merge)))
      colnames(genbank_metadata) = c("Sequence_ID","country","host","isolate","collection-date","isolation-source","BioSample","SRA","note")
      
      # Filling in columns whose values do not change or will be put in with options
      genbank_metadata$Sequence_ID = ncbi_gisaid_merge$Virus.name
      genbank_metadata$country = "USA"
      genbank_metadata$host = "Homo sapiens"
      genbank_metadata$isolate = ncbi_gisaid_merge$accession_short
      genbank_metadata$`collection-date`= ncbi_gisaid_merge$Collection.date
      genbank_metadata$`isolation-source` = "patient isolate"
      genbank_metadata$BioSample = ncbi_gisaid_merge$biosample_accession
      genbank_metadata$SRA = ncbi_gisaid_merge$accession
      genbank_metadata$note = paste("GISAID_ID:", ncbi_gisaid_merge$Accession.ID, sep = " ")
    }else{
        cat("ERROR: can't create metdata.  Please check column names.\n\n")
      }
  
  } else {
    cat("ERROR: check the spelling of your sequence type input: sequence type should match 'Illumina MiSeq', 'Illumina NextSeq', or 'Oxford Nanopore GridION'\n\n")
  }

} else {
  cat("ERROR: you didn't specify all varialbles. Please supply all values\n\n") # print error messages
}

write.table(genbank_metadata, file=paste("ncbi_genbank_submission", date, "metadata.tsv", sep="_"), row.names = FALSE, quote = TRUE, sep = '\t')

