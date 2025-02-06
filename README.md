The arboRist script was created in order to automate the process of searching the NCBI nucleotide database for nucleotide accessions and curating the associated metadata. 
Additionally, you can use this script to also produce input for an iqtree2 paritioned (multi-gene) analysis.

The NCBI database contains metadata of varying quality. This script attempts to standize the most common variations I have observed, however, there will almost certainly be unique adjustments needed for each analysis depending on your needs. I encourage you to read the data curation code and make adjustments that suit your purposes.

In order to use this script, please specify your desired options:

# Step 0: Install/load required packages

required_packages <- c("rentrez", "stringr", "plyr", "dplyr", "withr", "XML", 
                       "data.table", "tidyr", "phylotools", "scales", 
                       "purrr", "readr", "phytools", "RColorBrewer", 
                       "maps", "ggplot2", "tidygeocoder", "treeio", 
                       "ggtree", "ggrepel", "taxize", "Biostrings")
installed_packages <- required_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(required_packages[!installed_packages])
}
invisible(lapply(required_packages, library, character.only = TRUE))

# Step 1: Specify project options

#### Main Options
> base_dir <- "path.expand('~')" # change to whereever you want your project directory stored
> 
> project_name <- "Blackwellomyces_tree" # name your project
> 
> taxa_of_interest <- c("Blackwellomyces", "Flavocillium") # specify the taxa you want to grab accessions for
> 
> regions_to_include <- c("RPB2", "TEF", "ITS") # specify the regions you want to grab accessions for
> 
> min_region_requirement <- 3 # choose the minimum number of regions a strain needs to have in order to be included (defaults to the number of regions you selected)

#### optional: provide full path to a custom metadata file
> my_lab_sequences <- "" # each entry should have "Accession", strain", "sequence", "organism", and "gene" metadata provided at minimum

#### additional options

Here, I specified that I was only interested in accessions that were not genomes, scaffolds, or contigs, and were between 100-5000bp long.
I also HIGHLY recommended to provide your NCBI API key (increases number requests/second allowed).

> search_options <- "NOT (Contig[All Fields]) NOT (scaffold[All Fields]) NOT (genome[All Fields]) AND (biomol_genomic[PROP] AND (100[SLEN]:5000[SLEN])"
> 
> max_acc_per_taxa <- "max"
> 
> ncbi_api_key <- "" 
> 
> metadata_categories_keep <- c("GBSeq_locus", "GBSeq_length","GBSeq_strandedness","GBSeq_moltype",
                              "GBSeq_update.date","GBSeq_create.date","GBSeq_definition",
                              "GBSeq_accession.version","GBSeq_project","GBSeq_taxonomy",
                              "GBSeq_sequence","GBReference_title","GBSeq_feature.table",
                              "_title","_journal","ref_id","pubmed")
> 
> acc_to_exclude <- c("") # optionally remove selected accessions from metadata files

# Step 2: Set up project file structure and load packages 
> setup_project_structure()
> load_required_packages()


# Step 3:  Retrieve data from NCBI (accessions and metadata)
The metadata retrieval step in particular can take quite a long time, depending on your search terms. If you have your NCBI API key set, the script will retrieve the metadata at a rate of ~ 1 accession / 0.5 seconds. 

> ncbi_data_fetch(project_name, max_acc_per_taxa, taxa_of_interest)


# Step 4:  Curate and filter the data, then produce multifastas for each region you specified
> data_curate(project_name, taxa_of_interest, acc_to_exclude, min_region_requirement)


At this point the data-retrieval steps in arboRist are complete. You may continue this pipeline to run a iqtree2 partitioned analysis (multi-gene tree)


# Step 5:  Separate alignment of region fastas
At this point, arboRist does not automatically align your sequences for you. 
Do your own alignment (output aligned fasta format, include gaps). Save aligned files to: "./multifastas/aligned_fastas"
Eventually I'll get multiple formats supported.


# Step 6:  Run iqtree2 on single-gene multifasta files
(with modelfinder option)
(example provided is not run in R, you need to download iqtree2 separately and run on Terminal)

Like so:

> cd ./multifastas/aligned_fastas
> 
> for fasta in ./*.fas; do iqtree2 -redo -s $fasta -m MFP+MERGE -B 1000; done

iqtree2 citaiton:
Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams, M.D., Von Haeseler, A. and Lanfear, R., 2020. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Molecular biology and evolution, 37(5), pp.1530-1534.

Note: please cite the following paper if you use the model selection:
S. Kalyaanamoorthy, B.Q. Minh, T.K.F. Wong, A. von Haeseler, and L.S. Jermiin (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. Nat. Methods, 14:587–589. DOI: 10.1038/nmeth.4285

# Step 7: Prepare input files for iqtree2 partitioned analysis for multigene tree

This script automatically prepares the files as specified in the iqtree2 manual: http://www.iqtree.org/doc/Advanced-Tutorial#partitioned-analysis-for-multi-gene-alignments

> multigene_tree_prep(project_name)

If you use the partitioned models please cite:
O. Chernomor, A. von Haeseler, and B.Q. Minh (2016) Terrace aware data structure for phylogenomic inference from supermatrices. Syst. Biol., 65:997-1008. https://doi.org/10.1093/sysbio/syw037

# Step 8:  Running partitioned model analysis in iqtree2
(example provided is run on Terminal)

> cd ./multigene_tree
> 
> iqtree2 -redo -s ./prep/aligned_concatenated_sequences.fasta -p ./prep/multigene_nexus.nexus -b 1000 --prefix multigene

