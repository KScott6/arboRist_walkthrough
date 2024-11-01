# arboRist  -  semi-automatic NCBI nucleotide data retrieval and curation (plus additional iqtree2 partitioned analysis file preparation)

# main options
base_dir <- "path.expand('~')" # change to whereever you want your project directory stored
project_name <- "Blackwellomyces_tree" # name your project
taxa_of_interest <- c("Blackwellomyces", "Flavocillium") # specify the taxa you want to grab accessions for
regions_to_include <- c("RPB2", "TEF", "ITS") # specify the regions you want to grab accessions for
min_region_requirement <- 3 # choose the minimum number of regions a strain needs to have in order to be included (defaults to the number of regions you selected)


# optional: provide full path to custom metadata file
my_lab_sequences <- "" # each entry should have "Accession", strain", "sequence", "organism", and "gene" metadata provided at minimum

# additional options
search_options <- "NOT (Contig[All Fields]) NOT (scaffold[All Fields]) NOT (genome[All Fields]) AND (biomol_genomic[PROP] AND (100[SLEN]:5000[SLEN])"
max_acc_per_taxa <- "max"
ncbi_api_key <- "" # HIGHLY recommended to provide your NCBI API key (increases number requests/second allowed)
metadata_categories_keep <- c("GBSeq_locus", "GBSeq_length","GBSeq_strandedness","GBSeq_moltype",
                              "GBSeq_update.date","GBSeq_create.date","GBSeq_definition",
                              "GBSeq_accession.version","GBSeq_project","GBSeq_taxonomy",
                              "GBSeq_sequence","GBReference_title","GBSeq_feature.table",
                              "_title","_journal","ref_id","pubmed")
acc_to_exclude <- c("") # optionally remove selected accessions from metadata files



# automatic project folder setup
setup_project_structure <- function(base_dir,
                                    subdirs = c("intermediate_files", "metadata_files", "results_files", "temp_files", "multifastas",
                                                "multifastas/aligned_fastas", "multigene_tree", "multigene_tree/prep")) {
  project_dir <- file.path(base_dir, project_name)
  
  if (!dir.exists(project_dir)) {
    dir.create(project_dir)
  }
  
  for (dir in subdirs) {
    full_path <- file.path(project_dir, dir)
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = TRUE)
    }
  }
  setwd(project_dir)
}



get_sleep_duration <- function() {
  if (!is.null(ncbi_api_key)) {
    return(0.3)
  } else {
    return(0.5)
  }
}



load_required_packages <- function(packages = c("rentrez", "stringr", "plyr", "dplyr", "withr", "XML", 
                                                "data.table", "tidyr", "phylotools", "scales", 
                                                "purrr", "readr", "phytools", "RColorBrewer", 
                                                "maps", "ggplot2", "tidygeocoder", "treeio", 
                                                "ggtree", "ggrepel", "taxize", "Biostrings")) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}



# to retrieve accessions for individual taxa names
fetch_accessions_for_taxon <- function(taxon, max_acc = max_acc_per_taxa) {
  cat("Searching term:", taxon, "\n")
  filters <- paste("(", taxon, "[All Fields] AND Fungi[Organism])",")", sep="")
  # 10,000 accessions is the limit before the webhistory method needs to be involved
  search <- entrez_search(db="nucleotide", term=filters, retmax=9999)
  
  if ((length(search$ids)) == 9999) {
    cat(taxon, " has >10,000 NCBI accessions. Using webhistory to search, this may take a while.\n")
    large_search <- entrez_search(db="nucleotide", term=filters, use_history=TRUE)
    temp_filename <- paste0("./temp_files/temp_file_accessions_from_", taxon, ".txt")
    
    total_accession_count <- large_search[["count"]]
    max_acc <- ifelse(max_acc == "max", total_accession_count, as.numeric(max_acc))
    cat(total_accession_count, "accessions available for", taxon, ". A maximum of", max_acc, "accessions will be pulled.\n")
    
    for (seq_start in seq(1, max_acc, 50)) {
      recs <- entrez_fetch(db="nuccore", web_history=large_search$web_history, rettype="acc", retmax=50, retstart=seq_start)
      cat(recs, file=temp_filename, append=TRUE)
      cat(seq_start + 49, "accessions recorded\r")
    }
    
    large_temp_df <- read.table(temp_filename)
    colnames(large_temp_df) <- c("Accession")
    large_temp_df$genus <- taxon
    if (file.exists(temp_filename)) file.remove(temp_filename)
    cat("Accession retrieval for", taxon, "successful, temporary file removed\n\n")
    return(large_temp_df)
  } else {
    if (length(search$ids) <= 300) {
      summary <- entrez_summary(db="nuccore", id=search$ids)
    } else {
      summary <- list()
      index <- split(seq(1, length(search$ids)), ceiling(seq_along(seq(1, length(search$ids)))/300))
      for (p in index) {
        summary[p] <- entrez_summary(db="nuccore", id=search$ids[p])
      }
      class(summary) <- c("esummary_list", "list")
    }
    
    tempdf <- data.frame(Accession=unname(extract_from_esummary(summary, "caption")))
    tempdf$genus <- taxon
    cat("Search complete for", taxon, "\n")
    return(tempdf)
  }
}


# to iterate through all taxa names provided
get_accessions_for_all_taxa <- function(taxa_list, max_acc_per_taxa) {
  taxa_frame_acc <- list()
  
  for (i in seq_along(taxa_list)) {
    tryCatch({
      term <- taxa_list[i]
      tempdf <- fetch_accessions_for_taxon(term, max_acc = max_acc_per_taxa)
      
      # writing individual files (.csv) with accessions separate by taxa
      outfile_name <- paste0("./intermediate_files/Accessions_for_", term, ".csv")
      write.csv(tempdf, outfile_name, row.names=FALSE, quote=FALSE)
      taxa_frame_acc[[i]] <- tempdf
      
    }, error = function(e) {
      Sys.sleep(get_sleep_duration())
      cat("ERROR:", conditionMessage(e), "\n")
    })
  }
  
  # making file (.csv) with ALL retrieved accessions
  accession_list <- do.call(rbind, taxa_frame_acc)
  write.csv(accession_list, "./intermediate_files/all_pulled_accessions.csv", row.names=FALSE)
  cat("All accessions retrieved and saved.\n")
}


fetch_metadata_for_accession <- function(accession) {
  # getting data from NCBI nucleotide database using rentrez
  out.xml <- entrez_fetch(db="nuccore", id=accession, rettype="xml")
  list.out <- xmlToList(out.xml)
  accession_dfs <- lapply(list.out, data.frame, stringsAsFactors = FALSE)
  all_metadata_df <- accession_dfs[[1]]
  
  # filtering metadata to retain selected categories
  select_metadata_df <- all_metadata_df[,grep(paste(metadata_categories_keep, collapse = "|"), 
                                              x=names(all_metadata_df))]
  
  # separating metadata into basic and feature table categories
  basic_info_df <- select_metadata_df[,grep("GBSeq_locus|GBSeq_length|GBSeq_strandedness|GBSeq_update.date|GBSeq_create.date|GBSeq_definition|GBSeq_accession.version|GBSeq_project|GBSeq_organism|GBSeq_taxonomy|GBSeq_sequence", 
                                            x=names(select_metadata_df))]
  feature_table_df <- select_metadata_df[,grep("GBSeq_feature.table", x=names(select_metadata_df))]
  
  # transforming feature table columns containing "_name" as new column names with associated values
  name_cols <- which(grepl("_name", colnames(feature_table_df)))
  feature_transformed_df <- data.frame(feature_table_df[,name_cols + 1])
  colnames(feature_transformed_df) <- feature_table_df[,name_cols]
  
  # combining basic and transformed feature table data
  metadata_entry <- cbind(basic_info_df, feature_transformed_df)
  names(metadata_entry) <- sub("GBSeq_", "", names(metadata_entry))
  setnames(metadata_entry, old = "locus", new = "Accession")
  setnames(metadata_entry, old = "definition", new = "accession_title")
}

retrieve_ncbi_metadata <- function(project_name) {
  accession_list <- read.csv(paste0("./intermediate_files/all_pulled_accessions.csv"), header = TRUE)
  metadata_database_list <- list()
  
  for (i in 1:nrow(accession_list)) {
    tryCatch({
      accession <- accession_list$Accession[i]
      cat("Processing Accession:", accession, "\n")
      
      metadata_entry <- fetch_metadata_for_accession(accession)
      metadata_database_list[[i]] <- metadata_entry
      
      # Print small summary for the current accession metdata
      cat("#", i, "For nucleotide Accession:", metadata_entry$Accession, "retrieved information including:", "\n",
          "\t", "Species name:", metadata_entry$organism, "\n",
          "\t", "Isolation source:", metadata_entry$isolation_source, "\n",
          "\t", "Host information:", metadata_entry$host, "\n",
          "\t", "Nucleotide accession titled:", metadata_entry$accession_title, "\n")
    }, error = function(e) {
      Sys.sleep(ifelse(!is.null(ncbi_api_key), 0.3, 0.5))
      cat("ERROR:", conditionMessage(e), "\n")
    })
  }
  
  # Write metadata databases to CSV
  metadata_database <- do.call(rbind.fill, metadata_database_list)
  write.csv(metadata_database, paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, ".csv"), row.names = FALSE)
  cat("Metadata and publication data retrieved and saved for project:", project_name, "\n")
}



# if error on reading in custom file, add new line to end of custom metadata file (need to fix this)
merge_metadata_with_custom_file <- function(project_name) {
  metadata_file_path <- paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, ".csv")
  metadata_database <- read.csv(metadata_file_path, header = TRUE)
  custom_sequences <- read.csv(my_lab_sequences, header = TRUE, fill = TRUE)
  
  # Checking to make sure the custom file has the required columns (Accession, strain, sequence, organism, gene)
  if (!all(c("Accession", "strain", "sequence", "organism", "gene") %in% colnames(custom_sequences))) {
    stop("Custom file must contain at least 'accession', 'strain', 'sequence', 'organism', and 'gene' columns.")
  }
  
  # Merging the metadata with the custom sequences using rbind.fill
  merged_data <- rbind.fill(metadata_database, custom_sequences)
  
  # Writing the merged data to a new file
  write.csv(merged_data, paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, ".csv"), row.names = FALSE)
}



# Note: there's going to be issues if you are trying to incorporate "complex" sequences with >1 gene name (e.g. LSU/ITS1/ITS2/SSU). I'm still working on solving this. 
# Using previously fetched accessions, retrieve metadata and publication information for each accession
curate_metadata <- function(project_name) {
  metadata_file_path <- paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, ".csv")
  accession_list <- read.csv(metadata_file_path, header = TRUE)
  
  # curating strain names, making strain.standard. Prioritizing names from specimen_voucher, strain, isolate, Accession from most to least priority
  accession_list <- accession_list %>%
    mutate(strain.standard = coalesce(specimen_voucher, strain, isolate, Accession))
  
  # stripping special characters and whitespaces from strain.standard names - I see many strains where they differ based only on punctuation
  remove_char <- c("\\>", "\\<", "\\s+", ":", ";", "_", "-", "\\.", "\\(", "\\)", "&", "\\|")
  accession_list$strain.standard <- str_remove_all(accession_list$strain.standard, paste(remove_char, collapse = "|"))
  
  # making a separate strain name that indicates if the strain is TYPE (aka if it has anything recorded in the type_material column)
  accession_list$strain.standard.type <- ifelse(!is.na(accession_list$type_material), 
                                                paste(accession_list$strain.standard, ".TYPE", sep = ""), 
                                                accession_list$strain.standard)
  
  # Gene name standardization
  accession_list$gene.standard <- accession_list$gene
  accession_list <- accession_list %>%
    mutate(gene.standard = case_when(
      str_detect(gene.standard, regex('tef-*\\d*', ignore_case = TRUE)) ~ "TEF",
      str_detect(gene.standard, regex('EF1-alpha', ignore_case = TRUE)) ~ "TEF",
      str_detect(gene.standard, regex('ef1a', ignore_case = TRUE)) ~ "TEF",
      str_detect(gene.standard, regex('b-tub', ignore_case = TRUE)) ~ "BTUB",
      str_detect(gene.standard, regex('TUB2', ignore_case = TRUE)) ~ "BTUB",
      str_detect(gene.standard, regex('RBP2', ignore_case = TRUE)) ~ "RPB2",
      str_detect(gene.standard, regex('RPB2', ignore_case = TRUE)) ~ "RPB2",
      str_detect(gene.standard, regex('RBP1', ignore_case = TRUE)) ~ "RPB1",
      str_detect(gene.standard, regex('RPB1', ignore_case = TRUE)) ~ "RPB1",
      str_detect(gene.standard, regex('act[16]', ignore_case = TRUE)) ~ "actin",
      TRUE ~ gene.standard
    ))
  
  # Product name standardization
  accession_list$product.standard <- accession_list$product
  accession_list <- accession_list %>%
    mutate(product.standard = case_when(
      str_detect(product.standard, regex('elongation factor 1', ignore_case = TRUE)) ~ "TEF",
      str_detect(product.standard, regex('18S', ignore_case = TRUE)) ~ "SSU",
      str_detect(product.standard, regex('16S', ignore_case = TRUE)) ~ "SSU",
      str_detect(product.standard, regex('26S', ignore_case = TRUE)) ~ "LSU",
      str_detect(product.standard, regex('28S', ignore_case = TRUE)) ~ "LSU",
      str_detect(product.standard, regex('actin beta', ignore_case = TRUE)) ~ "actin",
      str_detect(product.standard, regex('beta-tubulin', ignore_case = TRUE)) ~ "BTUB",
      str_detect(product.standard, regex('licensing\\D*7\\D*', ignore_case = TRUE)) ~ "MCM7",
      str_detect(product.standard, regex('polymerase II larg[est]*', ignore_case = TRUE)) ~ "RPB1",
      str_detect(product.standard, regex('polymerase II second largest', ignore_case = TRUE)) ~ "RPB2",
      str_detect(product.standard, regex('small subunit ribosomal', ignore_case = TRUE)) ~ "SSU",
      str_detect(product.standard, regex('^large[st]* subunit ribosomal', ignore_case = TRUE)) ~ "LSU",
      str_detect(product.standard, regex('internal transcribed', ignore_case = TRUE)) ~ "ITS",
      TRUE ~ product.standard
    ))
  
  # Accession title standardization
  accession_list$acc_title.standard <- accession_list$accession_title
  accession_list <- accession_list %>%
    mutate(acc_title.standard = case_when(
      str_detect(acc_title.standard, regex('actin beta', ignore_case = TRUE)) ~ "actin",
      str_detect(acc_title.standard, regex('beta-*\\s*tubulin', ignore_case = TRUE)) ~ "BTUB",
      str_detect(acc_title.standard, regex('licensing\\D*7\\D*', ignore_case = TRUE)) ~ "MCM7",
      str_detect(acc_title.standard, regex('RPB1', ignore_case = TRUE)) ~ "RPB1",
      str_detect(acc_title.standard, regex('RPB2', ignore_case = TRUE)) ~ "RPB2",
      str_detect(acc_title.standard, regex('internal transcribed.*,*\\s*complete sequence\\;', ignore_case = TRUE)) ~ "ITS",
      str_detect(acc_title.standard, regex('translation elongation', ignore_case = TRUE)) ~ "TEF",
      TRUE ~ acc_title.standard
    ))
  
  # Final region name assignment
  accession_list <- accession_list %>%
    mutate(region.standard = coalesce(gene.standard, product.standard, acc_title.standard))
  
  # Creating the custom FASTA headers
  accession_list$org_name <- accession_list$organism
  accession_list$org_name <- gsub('\\s+', '\\.', accession_list$org_name)
  accession_list$fasta.header <- paste(">", accession_list$org_name, "_", accession_list$strain.standard, sep = "")
  accession_list$fasta.header.type <- paste(">", accession_list$org_name, "_", accession_list$strain.standard.type, sep = "")
  
  # Save the curated metadata
  write.csv(accession_list, paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, "_curated.csv"), row.names = FALSE)
}



filter_metadata <- function(project_name, taxa_of_interest, acc_to_exclude = NULL) {
  accession_list <- read.csv(paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, "_curated.csv"), header = TRUE)
  # Getting reported genus name
  accession_list$complex_name <- accession_list$organism
  accession_list <- accession_list %>% 
    separate(complex_name, c('listed.genus', 'V1', 'V2', 'V3'), sep = " ", extra = 'merge')
  accession_list <- subset(accession_list, select = -c(V1, V2, V3))
  # filtering based on 'taxa_of_interest'
  accession_list <- accession_list[accession_list$listed.genus %in% taxa_of_interest, ]
  # Optional filtering of specified accessions
  if (!is.null(acc_to_exclude) && length(acc_to_exclude) > 0 && any(acc_to_exclude != "")) {
    accession_list <- accession_list[!accession_list$Accession %in% acc_to_exclude, ]
  }
  # Handling invalid UTF-8 characters - accents from NCBI?
  accession_list[] <- lapply(accession_list, function(x) iconv(x, from = "UTF-8", to = "UTF-8", sub = "byte"))
  
  write.csv(accession_list, 
            paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, "_curated_filtered.csv"), 
            row.names = FALSE)
}



select_regions <- function(project_name, min_region_requirement = length(regions_to_include)) {
  accession_list <- read.csv(paste0("./metadata_files/all_accessions_pulled_metadata_", project_name, "_curated_filtered.csv"), header = TRUE)
  
  # filter so only accessions with region.standard values in 'regions_to_include' are kept
  multifasta_prep <- accession_list[accession_list$region.standard %in% regions_to_include, ]
  write.csv(multifasta_prep, paste0("./metadata_files/selected_accessions_metadata_", project_name, ".csv"), row.names = FALSE)
  
  # creating a simplified dataframe
  multifasta_prep_complete <- subset(multifasta_prep, select = c(strain.standard.type, organism, Accession, region.standard))
  
  # generate a complete region attendance sheet - useful for publishing
  complete_region_attendance <- multifasta_prep_complete %>%
    pivot_wider(
      names_from = "region.standard",
      values_from = "Accession", values_fn = list)
  complete_region_attendance <- complete_region_attendance %>%
    mutate(across(all_of(regions_to_include), ~replace(., lengths(.) == 0, NA)))
  
  # filtering duplicate accessions for a clean region attendance sheet
  multifasta_prep_select <- distinct(multifasta_prep_complete, strain.standard.type, region.standard, .keep_all= TRUE)
  select_region_attendance <- multifasta_prep_select %>%
    pivot_wider(
      names_from = "region.standard",
      values_from = "Accession")
  
  # filter strains by minimum region inclusion threshold
  select_region_attendance_filtered <- select_region_attendance %>%
    mutate(total=rowSums(!is.na(select(., -strain.standard.type, -organism)))) %>%
    filter(total >= min_region_requirement) %>%
    select(-total)
  
  write.csv(select_region_attendance_filtered, paste0("./metadata_files/Region_attendance_sheet_", project_name, ".csv"), row.names = FALSE)
}



create_multifastas <- function(project_name, output_dir = "./multifastas/") {
  # Reading in filtered region data
  select_region_attendance_filtered <- read.csv(paste0("./metadata_files/Region_attendance_sheet_", project_name, ".csv"), header = TRUE)
  filtered_accessions <- data.frame(Accession = unlist(select_region_attendance_filtered[,-c(1:2)]), row.names = NULL)
  
  # Reading in curated metadata
  accession_list <- read.csv(paste0("./metadata_files/all_accessions_pulled_metadata_",project_name,"_curated_filtered.csv"), header = TRUE)
  filtered_accessions_metadata <- merge(filtered_accessions, accession_list, by = "Accession", all.x = TRUE)
  
  # prep data for fasta file creation
  multifasta_prep_simple <- subset(filtered_accessions_metadata, 
                                   select = c(strain.standard, organism, Accession, region.standard, fasta.header.type, sequence))
  region_dfs <- split(multifasta_prep_simple, with(multifasta_prep_simple, region.standard), drop = TRUE)
  
  # double check output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # creating separate multifasta file, per region
  for (i in seq_along(region_dfs)) {
    locus_df <- region_dfs[[i]]
    locus_df <- locus_df[order(locus_df$fasta.header.type, decreasing = FALSE), ]
    region.name <- unique(locus_df$region.standard)
    seqs_fasta <- c(rbind(locus_df$fasta.header.type, locus_df$sequence))
    filename <- paste0(output_dir, region.name, ".fasta")
    write.table(seqs_fasta, filename, row.names = FALSE, quote = FALSE, col.names = FALSE)
    
    message("Created multifasta for region: ", region.name)
  }
  cat("Multifasta files created in ", output_dir)
}


# To cncatenate aligned FASTA sequences by region
concatenate_fasta_sequences <- function(project_name) {
  # Set directory and regions to include
  fasta_dir <- paste0("./multifastas/aligned_fastas/")
  fasta_files <- list.files(fasta_dir, pattern = "\\.fas$", full.names = TRUE)
  
  # Check if all regions in 'regions_to_include' have a corresponding FASTA file
  matched_files <- fasta_files[sapply(regions_to_include, function(region) any(grepl(region, fasta_files)))]
  if (length(matched_files) < length(regions_to_include)) {
    warning("Some regions in 'regions_to_include' do not have matching FASTA files.")
  }
  
  aligned_sequences <- list()
  # looping over each FASTA (.fas) file for the regions of interest
  for (region_file in matched_files) {
    region_name <- str_match(basename(region_file), "(\\w+)\\.fasta")[,2]
    sequences <- readDNAStringSet(region_file)
    region_df <- data.frame(strain = names(sequences), sequence = as.character(sequences), stringsAsFactors = FALSE)
    colnames(region_df)[2] <- region_name  # Name the sequence column by region
    aligned_sequences[[region_name]] <- region_df
  }
  
  # merging all regions by strain
  all_regions_df <- purrr::reduce(aligned_sequences, full_join, by = "strain")
  
  # concatenating sequences from each region for each strain, handling missing data
  concat_sequences <- all_regions_df %>%
    mutate(concat_sequence = apply(select(., -strain), 1, function(row) {
      paste(row[!is.na(row)], collapse = "")
    })) %>%
    select(strain, concat_sequence)
  
  # writing concatenated sequences to a new FASTA file
  fasta_out_path <- paste0("./multigene_tree/prep/aligned_concatenated_sequences.fasta")
  writeXStringSet(DNAStringSet(setNames(concat_sequences$concat_sequence, concat_sequences$strain)), fasta_out_path)
  
  cat("Concatenated FASTA file saved to:", fasta_out_path, "\n")
}



create_nexus_file <- function() {
  # list of IQ-TREE log files
  log_dir = "./multifastas/aligned_fastas/"
  output_file = output_file = "./multigene_tree/prep/multigene_nexus.nexus"
  logfile_list <- list.files(log_dir, pattern = "*.iqtree", full.names = TRUE)
  nexus_data_list <- list()
  
  # Loop over each iqtree2 logfile (single gene) to extract region information, sequence length, and model
  for (i in seq_along(logfile_list)) {
    # extracting region name from file name
    region_name_match <- str_match(logfile_list[[i]], "(\\w+)_aligned")
    region_name <- region_name_match[2]
    logfile <- read_file(logfile_list[[i]])
    
    # getting sequence length information
    seq_length_match <- str_match(logfile, "Input\\sdata:\\s(\\d+)\\ssequences\\swith\\s(\\d+)\\snucleotide\\ssites")
    seq_length <- as.numeric(seq_length_match[3])
    
    # getting best-fit model according to BIC
    sg_model_match <- str_match(logfile, "Best-fit\\smodel\\saccording\\sto\\sBIC:\\s(.+)\\n")
    sg_model <- sg_model_match[2]
    
    # Storing info in a data frame for each region
    nexus_data_list[[region_name]] <- data.frame(region_num = i, sg_model = sg_model, seq_length = seq_length)
  }
  
  nexus_data <- do.call(rbind, nexus_data_list)
  # calculating character partition sets
  nexus_data$partition_begin <- cumsum(c(1, head(nexus_data$seq_length, -1)))
  nexus_data$partition_end <- cumsum(nexus_data$seq_length)
  # writig the Nexus file header
  write.table("", file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("#nexus\nbegin sets;\n", file = output_file, append = TRUE)
  
  # writing charset definitions for each region
  for (i in 1:nrow(nexus_data)) {
    cat("\tcharset part", nexus_data$region_num[i], " = ", nexus_data$partition_begin[i], "-", nexus_data$partition_end[i], ";\n",
        file = output_file, sep = "", append = TRUE)
  }
  
  # writing character partitions
  cat("\tcharpartition mine = ", file = output_file, sep = "", append = TRUE)
  for (i in 1:nrow(nexus_data)) {
    cat(nexus_data$sg_model[i], ":part", nexus_data$region_num[i], if (i < nrow(nexus_data)) ", " else ";",
        file = output_file, sep = "", append = TRUE)
  }
  # finish nexus file
  cat("\nend;", file = output_file, append = TRUE)
  
  cat("Nexus file created at:", output_file, "\n")
}



# wrapper for data retrieval functions
ncbi_data_fetch <- function(project_name, max_acc_per_taxa = 'max', taxa_of_interest = NULL) {
  # getting list(s) of accessions
  get_accessions_for_all_taxa(taxa_of_interest, max_acc_per_taxa)
  # fetching metadata for accessions
  retrieve_ncbi_metadata(project_name)
}


# wrapper for all curation functions
data_curate <- function(project_name, taxa_of_interest = NULL, acc_to_exclude = NULL, min_region_requirement = length(regions_to_include)) {
  # merge metadata with custom file
  merged_data <- merge_metadata_with_custom_file(project_name)
  # curate metadata
  curated_data <- curate_metadata(project_name)
  # filter metadata based on taxa and accession exclusions
  filtered_data <- filter_metadata(project_name, taxa_of_interest, acc_to_exclude)
  # select based on specific regions
  selected_regions_data <- select_regions(project_name, min_region_requirement)
  # create multi-FASTA files
  create_multifastas(project_name)
}


#wrapper for multigene tree prep functions
multigene_tree_prep <- function(project_name) {
  concatenate_fasta_sequences(project_name)
  create_nexus_file()
}
