# FoldDelay_script.R
# Set the locale to English
Sys.setenv(LANGUAGE = "en")
# Suppress warning messages
options(warn = -1)
# First of all, install required packages
list.of.packages <- c("ggplot2", "bio3d", "dplyr", "tidyr", "stringr", "tibble","ggpubr","jsonlite","optparse","curl")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = c(
  binary = "https://packagemanager.rstudio.com/all/__linux__/focal/latest",
  source = "https://packagemanager.rstudio.com/all/latest",
  CRAN = "https://cloud.r-project.org"
))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(bio3d))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(curl))
code3 <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
           "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
           "TYR", "VAL","H1S","H2S","TPO","PTR","SEP","HYP")
code1 <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K",
           "M", "F", "P", "S", "T", "W", "Y", "V","H","H","T","Y","S","P")
# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to a text file containing UniProt identifiers", metavar = "character"),
  make_option(c("-k", "--keep"), type = "logical", default = FALSE,
              help = "Set to TRUE if keep all contacts. Otherwise, only the most C-terminal partner per residue is kept", metavar = "BOOLEAN"),
  make_option(c("-t", "--transrate"), type = "numeric", default = 5,
              help = "Set the translation rate (aa/s)", metavar = "NUMBER"),
  make_option(c("-d", "--distance"), type = "numeric", default = 6,
              help = "Distance cutoff (Å) within which a contact can be formed", metavar = "NUMBER")
)
# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
transrate <- opt$transrate
distance <- opt$distance
keep <- opt$keep
# Convert defaults from factory server
if (!is.null(opt$input) && opt$input == "None") {
  opt$input <- NULL
}
# If no file is provided give an error
if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("You must provide a text file containing UniProt identifier.", call. = FALSE)
}
#Read input file
uniprot_ids <- read.table("~/Downloads/uniprot_ids.txt", quote="\"", comment.char="")
#Loop through all the structures
fold_delay_list <- vector("list", length = nrow(uniprot_ids)) 
for (i in 1:nrow(uniprot_ids)){
  temp_uid <- toupper(trimws(uniprot_ids$V1[i]))
  v_found <- NA_integer_
  file_name <- af_url <- NULL
  # try v20..v1; bump upper bound if needed
  for (vv in 20:1) {
    test_name <- sprintf("AF-%s-F1-model_v%d.pdb", temp_uid, vv)
    test_url  <- paste0("https://alphafold.ebi.ac.uk/files/", test_name)
    ok <- tryCatch({
      h <- curl::new_handle(
        nobody = TRUE,
        customrequest = "HEAD",
        followlocation = TRUE,
        useragent = "af-fetch/1.0"
      )
      res <- curl::curl_fetch_memory(test_url, handle = h)
      res$status_code >= 200 && res$status_code < 400
    }, error = function(e) FALSE)
    if (ok) {
      file_name <- test_name;
      af_url <- test_url;
      v_found <- vv;
      break
    }
  }
  if (is.na(v_found)) {
    message(sprintf("No AlphaFold PDB found for %s (v20..v1).", uid), call. = FALSE)
    next
  }
  message("Fetching structure from: ", af_url)
  res <- curl::curl_fetch_memory(af_url)
  # sanity check (keep your own checks too)
  if (is.null(res$content) || length(res$content) == 0) {
    message("Empty download: ", af_url, call. = FALSE)
    next
  }
  tmp <- tempfile(fileext = ".pdb")
  writeBin(res$content, tmp)
  pdb <- bio3d::read.pdb(
    tmp,
    maxlines = -1,
    multi = FALSE,
    rm.insert = FALSE,
    rm.alt = TRUE,
    ATOM.only = FALSE,
    hex = FALSE,
    verbose = TRUE
  )
  unlink(tmp)  # delete temp file immediately
  pdb_ca <- pdb$atom
  pdb_ca <- subset(pdb_ca, elety == "CA", select = c(5,6,7,9,10,11,13))
  for (k in 1:length(code3)){
    pdb_ca$resid <- gsub(code3[k],code1[k],pdb_ca$resid,ignore.case=TRUE)
  }
  length=max(pdb_ca$resno)
  colnames(pdb_ca)[c(1,3)] <- c("Amino_acid","Index")
  pdb_ca <- subset(pdb_ca, select = c(1,3))
  cat("File processing complete.\n")
  # Calculate contact map
  chimaerarainbow <- colorRampPalette(c("blue","cyan","green","yellow","red"))
  ref.cont <- cmap(pdb, dcut=distance, scut=1,verbose=FALSE)
  data = as.data.frame(ref.cont)
  names(data) <- seq(1,length)
  data$from = seq(1,length)
  simple <- gather(data,"to","value",-from)%>%
    mutate(color = rep(chimaerarainbow(length),length))%>%
    na.omit() %>%
    filter(value !=0) %>%
    select("to","from","color")
  colnames(simple)[2] <- "Index"
  simple <- merge(simple, pdb_ca)
  colnames(simple) <- c("index_1","Index","color","aa_1")
  simple <- merge(simple, pdb_ca)
  colnames(simple) <- c("index_2","index_1","color","aa_1","aa_2")
  simple <- subset(simple, select = c(2,4,1,5,3))
  simple$index_1 <- as.numeric(simple$index_1)
  simple$index_2 <- as.numeric(simple$index_2)
  simple <- simple[order(simple$index_1,simple$index_2),]
  simple$distance <- simple$index_2 - simple$index_1
  simple <- simple %>%
    arrange(index_1)
  # Donwload PAE 
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-",temp_uid,"-F1-predicted_aligned_error_v",v_found,".json")
  res <- curl::curl_fetch_memory(url)
  # sanity check
  if (is.null(res$content) || length(res$content) == 0) {
    message("Empty PAE download: ", url)
    next
  }
  tmp_json <- tempfile(fileext = ".json")
  writeBin(res$content, tmp_json)
  json_data <- jsonlite::fromJSON(tmp_json)
  json_matrix <- json_data$predicted_aligned_error[[1]]
  unlink(tmp_json)  # delete temp file immediately
  simple$flag <- "No"
  # Add domain information
  # Extract the 'cath_domains' list
  url <- paste0("https://ted.cathdb.info/api/v1/uniprot/summary/",temp_uid,"?skip=0&limit=100")
  h <- new_handle()
  handle_setheaders(h, accept = "application/json")
  con <- curl(url, handle = h)
  json_data <- fromJSON(readLines(con, warn = FALSE))
  close(con)
  domains <- json_data$data
  domains <- domains[order(as.numeric(sub(".*TED", "", domains$ted_id))), ]
  #Initialise variables
  simple$Cath_label <- NA
  simple$domain_id <- NA
  simple$Cath_label_2 <- NA
  simple$domain_id_2  <- NA
  simple$Scope <- "–"
  # If there are no domains, just keep domain_ranges empty and continue
  if (is.null(domains) || length(domains) == 0) {
    domain_ranges <- data.frame(
      domain_id   = integer(0),
      start       = integer(0),
      end         = integer(0),
      Cath_label  = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    # Parse information on domains
    domains <- subset(domains, select = c(5,15))
    expand_chopping_ranges <- function(domains) {
      expanded <- do.call(rbind, lapply(1:nrow(domains), function(i) {
        chopping_str <- domains$chopping[i]
        Cath_label <- domains$cath_label[i]
        domain_id <- i
        chopping_str <- gsub("[\u2013\u2014]", "-", chopping_str)
        if (is.na(chopping_str) || chopping_str == "") return(NULL)
        parts <- unlist(strsplit(chopping_str, "_"))
        rows <- lapply(parts, function(p) {
          bounds <- as.numeric(unlist(strsplit(p, "-")))
          if (length(bounds) != 2 || any(is.na(bounds))) return(NULL)
          data.frame(
            domain_id = domain_id,
            start = bounds[1],
            end = bounds[2],
            Cath_label = Cath_label,
            stringsAsFactors = FALSE
          )
        })
        rows <- Filter(Negate(is.null), rows)
        if (length(rows) == 0) return(NULL)
        do.call(rbind, rows)
      }))
      return(expanded)
    }
    domain_ranges <- expand_chopping_ranges(domains)
  }
  # Add domain information to dataframe
  for (j in seq_len(nrow(domain_ranges))) {
    in_range_1 <- simple$index_1 >= domain_ranges$start[j] &
      simple$index_1 <= domain_ranges$end[j]
    in_range_2 <- simple$index_2 >= domain_ranges$start[j] &
      simple$index_2 <= domain_ranges$end[j]
    simple$Cath_label[in_range_1] <- domain_ranges$Cath_label[j]
    simple$domain_id[in_range_1]  <- domain_ranges$domain_id[j]
    simple$Cath_label_2[in_range_2] <- domain_ranges$Cath_label[j]
    simple$domain_id_2[in_range_2]  <- domain_ranges$domain_id[j]
  }
  simple$Cath_label <- sub(",.*", "", simple$Cath_label) #Removes cases with multiple labels
  # Determine interaction scope
  both_in_domain <- !is.na(simple$domain_id) & !is.na(simple$domain_id_2)
  simple$Scope[both_in_domain & simple$domain_id == simple$domain_id_2] <- "intra"
  simple$Scope[both_in_domain & simple$domain_id != simple$domain_id_2] <- "inter"
  # Flag interactions with an expected position error (based on PAE) > 6 Å
  for(j in 1:nrow(simple)){
    temp_from <- simple$index_1[j]
    temp_to <- simple$index_2[j]
    temp_pae <- json_matrix[temp_from,temp_to]
    if (temp_pae > distance) {
      simple$flag[i] <- "Yes"
    }
  }
  simple <- subset(simple, flag == "No")
  # Table for FD per-residue profile
  simple$time <- simple$distance / transrate
  simple$uid <- temp_uid
  # Keep only most C-terminal partner per residue 
  if (keep == FALSE){
    simple <- simple %>%
      group_by(index_1) %>%
      slice_max(order_by = index_2, n = 1, with_ties = FALSE) %>%
      ungroup()
  }
  fold_delay_nopae <- subset(simple, select = c(14,1:4,6,13,8,9,12))
  colnames(fold_delay_nopae) <- c("UniProt id","Index 1","Aa 1","Index 2","Aa 2","NFD (aa)","NFD (s)","Cath label","Domain id","Scope")
  fold_delay_list[[i]] <- fold_delay_nopae
  rm(pdb, simple)
  gc()
} 
fold_delay_list <- bind_rows(fold_delay_list)
write.csv(
  fold_delay_list,
  sep = ";",
  file = "folddelay.csv",
  row.names = FALSE,
  quote = FALSE
)
cat("Output file generated.\n")