 
  ### Updated R Script with Parameter Values
  
  ```r
# Load necessary libraries
library(dplyr)
library(tidyr)
library(purrr)
library(rtracklayer)  # For reading GTF files

# Function to create directories
create_directories <- function() {
  dir.create("7-MetagTbl", showWarnings = FALSE)
  dir.create("7-MetagTbl/Reports", showWarnings = FALSE)
}

# Function to read in and process data
process_metagene_data <- function(SRAList, Names, Params) {
  
  # Create necessary directories
  create_directories()
  
  # Assigning parameters from Params
  rlmin <- as.integer(Params$ReadLenMiN)  # Minimum read length
  rlmax <- as.integer(Params$ReadLenMaX)  # Maximum read length
  Span <- as.integer(Params$MetagSpan)    # Number of nucleotides before/after start/stop codons
  Mapping <- Params$Mapping                # Mapping strategy (5' or 3')
  dataNorm <- Params$Normalised            # "raw" or "rpm" normalization
  MappedTwice <- Params$MappedTwice        # Whether to include reads mapped once or twice
  MetagThreshold <- as.integer(Params$MetagThreshold)  # Minimum number of reads
  MetagSpancorr <- as.integer(Params$MetagSpancorr)  # Corrected span for metagenomics
  GeneRawMeanThr <- as.numeric(Params$GeneRawMeanThr)  # Gene raw mean threshold
  CodonRawMeanThr <- as.numeric(Params$CodonRawMeanThr)  # Codon raw mean threshold
  
  # For the normalized data, adjust thresholds
  if (dataNorm == "rpm") {
    norm_factor <- Params$NormalizationFactor  # Adjust this parameter based on your RPM data
    GeneRpmMeanThr <- GeneRawMeanThr / norm_factor
    CodonRpmMeanThr <- CodonRawMeanThr / norm_factor
  }
  
  columns <- as.character(seq(rlmin, rlmax)) # Column names (from min read length to max read length)
  columns <- c(columns, "sum")              # Add "sum" column
  
  rlrange <- paste(rlmin, rlmax, sep="-")    # Range of read lengths for filenames
  LogFileName <- paste0("7-MetagTbl/Reports/MetagTabl_", Mapping, "-End_", rlrange, "_iv_log.txt")
  LOG_FILE <- file(LogFileName, "wt")
  
  for (iN in Names) {
    # Initialize counters
    cf1 <- cr1 <- cf2 <- cr2 <- 0  
    report <- paste("\nName:", iN)
    cat(report, "\n")
    writeLines(report, LOG_FILE)
    
    # Generate file names based on the name and parameters
    fn_body <- paste(iN, Mapping, "-End_", rlrange, sep="_")
    outf_start <- paste0("7-MetagTbl/", fn_body, "_", dataNorm, "_Start_iv_Meta_Sum.txt")
    outf_stop <- paste0("7-MetagTbl/", fn_body, "_", dataNorm, "_Stop_iv_Meta_Sum.txt")
    infile_h5 <- paste0("6-AssignRaw/", fn_body, "_iv.h5")
    
    # Initialize empty data frames for start and stop codons
    meta_start_dff <- data.frame(matrix(ncol = length(columns), nrow = 2*Span + 1))
    meta_start_dfr <- data.frame(matrix(ncol = length(columns), nrow = 2*Span + 1))
    meta_stop_dff <- data.frame(matrix(ncol = length(columns), nrow = 2*Span + 1))
    meta_stop_dfr <- data.frame(matrix(ncol = length(columns), nrow = 2*Span + 1))
    
    colnames(meta_start_dff) <- colnames(meta_start_dfr) <- colnames(meta_stop_dff) <- colnames(meta_stop_dfr) <- columns
    
    # Annotation file (adjust with actual GTF file)
    tabixfile <- rtracklayer::import("0-References/genome.gtf.gz", format = "GTF")  # GTF file for gene annotations
    Threshold <- MetagThreshold  # Threshold for gene inclusion
    
    # Adjust for RPM if normalization is set to "rpm"
    if (dataNorm == "rpm") {
      BamName <- paste0("5-Aligned/", iN, ".bam")  # BAM file
      Threshold <- raw_metag_threshold_to_rpm(BamName, Threshold, Params)  # Function to adjust threshold for RPM
    }
    
    report <- paste("Threshold for Metagene", Threshold, dataNorm)
    cat(report, "\n")
    writeLines(report, LOG_FILE)
    
    # Data frames for forward and reverse strands
    Fkey <- ifelse(dataNorm == "rpm", "For_rpm", "For_raw")
    Rkey <- ifelse(dataNorm == "rpm", "Rev_rpm", "Rev_raw")
    df_f <- readHDF5File(infile_h5, Fkey)  # Read forward strand data
    df_r <- readHDF5File(infile_h5, Rkey)  # Read reverse strand data
    
    for (Chr in unique(tabixfile$seqnames)) {
      # Get data for each chromosome
      dfc_f <- df_f[df_f$Chr == Chr, columns]
      dfc_r <- df_r[df_r$Chr == Chr, columns]
      
      # Loop through each GTF feature (start/stop codons)
      for (gtf in tabixfile[tabixfile$seqnames == Chr]) {
        # Process Stop codons on Forward strand
        if (gtf$type == 'stop_codon' && gtf$strand == '+') {
          df <- dfc_f[gtf$start - Span : gtf$start + Span, ]
          if (nrow(df) < Threshold) next  # Check for sufficient data
          meta_stop_dff <- meta_stop_dff + df  # Add to sum
          cf2 <- cf2 + 1
        }
        
        # Process Stop codons on Reverse strand
        if (gtf$type == 'stop_codon' && gtf$strand == '-') {
          df <- dfc_r[gtf$end - Span - 1 : gtf$end + Span - 1, ]
          if (nrow(df) < Threshold) next  # Check for sufficient data
          meta_stop_dfr <- meta_stop_dfr + df  # Add to sum
          cr2 <- cr2 + 1
        }
        
        # Process Start codons on Forward strand
        if (gtf$type == 'start_codon' && gtf$strand == '+') {
          df <- dfc_f[gtf$start - Span : gtf$start + Span, ]
          if (nrow(df) < Threshold) next  # Check for sufficient data
          meta_start_dff <- meta_start_dff + df  # Add to sum
          cf1 <- cf1 + 1
        }
        
        # Process Start codons on Reverse strand
        if (gtf$type == 'start_codon' && gtf$strand == '-') {
          df <- dfc_r[gtf$end - Span - 1 : gtf$end + Span - 1, ]
          if (nrow(df) < Threshold) next  # Check for sufficient data
          meta_start_dfr <- meta_start_dfr + df  # Add to sum
          cr1 <- cr1 + 1
        }
      }
    }
    
    # Summing up the data
    meta_start_sum <- meta_start_dff + meta_start_dfr
    meta_stop_sum <- meta_stop_dff + meta_stop_dfr
    
    # Save to files
    write.table(meta_start_sum, outf_start, sep = "\t", row.names = TRUE, col.names = TRUE)
    write.table(meta_stop_sum, outf_stop, sep = "\t", row.names = TRUE, col.names = TRUE)
    
    report <- paste("Around START:", cf1 + cr1, "included", outf_start)
    cat(report, "\n")
    writeLines(report, LOG_FILE)
    
    report <- paste("Around STOP:", cf2 + cr2, "included", outf_stop)
    cat(report, "\n")
    writeLines(report, LOG_FILE)
  }
  
  # Close the log file
  close(LOG_FILE)
}

# Set up the parameter values (directly from your input)
Params <- list(
  ReadLenMiN = 25,         # Minimum read length
  ReadLenMaX = 35,         # Maximum read length
  MetagSpan = 60,          # Number of nucleotides before/after start/stop in metagene profiles
  Mapping = "5",           # Mapping strategy: "5" for 5' end
  Normalised = "rpm",      # Normalization method: "rpm"
  MappedTwice = "No",      # Include reads mapped once or twice
  MetagThreshold = 30,     # Minimum reads required for inclusion
  MetagSpancorr = 90,      # Corrected span for metagenomic profile
  GeneRawMeanThr = 0.3,    # Gene raw mean threshold (raw counts per 100 nt)
  CodonRawMeanThr = 1.6,   # Codon raw mean threshold (raw counts per

