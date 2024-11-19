input_fasta <- as.character(commandArgs(trailingOnly = TRUE)[1])
input_gff <- as.character(commandArgs(trailingOnly = TRUE)[2])
input_bam <- as.character(commandArgs(trailingOnly = TRUE)[3])
sample <-as.character(commandArgs(trailingOnly = TRUE)[4])
threshold <- as.numeric(commandArgs(trailingOnly = TRUE)[5])
N.threads <- as.integer(commandArgs(trailingOnly = TRUE)[6])

suppressMessages({
  suppressWarnings({
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    
    if (!require("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", quiet = TRUE)
    }
    
    required_packages <- c("dplyr", "foreach", "doParallel", "data.table")
    for(package in required_packages) {
      if(!require(package, character.only = TRUE, quietly = TRUE)) {
        install.packages(package, quiet = TRUE)
      }
    }
    
    if(!require("rtracklayer", quietly = TRUE)) {
      BiocManager::install("rtracklayer", quiet = TRUE)
    }
  })
})

# Load all packages with suppressed messages
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  library(rtracklayer)
  library(foreach)
})


### Load input data
gff <- rtracklayer::import(input_gff)
gff <- as.data.frame(gff)
interesting_features <- c("CDS")
feature_types <- unique(gff$type)
interesting_list <- intersect(interesting_features, feature_types)
gff <- gff[gff$type %in% interesting_list, ]
list_chr <- unique(gff$seqnames)



duplicated_columns <- names(gff)[duplicated(names(gff))]
gff <- gff[, !names(gff) %in% duplicated_columns]

# Group by locus_tag and filter the row with the maximum width
gff <- gff %>%
  group_by(locus_tag) %>%
  filter(end - start == max(end - start))
gff <- as.data.frame(gff)

#gff$locus_tag <- gsub(".*locus_tag=([^;]+);.*", "\\1", gff$locus_tag)

#### Subset convergent and divergent pairs of genes
#gff <- gff[gff$seqnames == list_chr[1],]
gff <- gff[, c("seqnames", "start", "end", "strand", "type", "locus_tag", "gene")]
convergent_genes <- list()
for (i in 1:(nrow(gff) - 1)) {
  if (gff$strand[i] == "+" & gff$strand[i + 1] == "-") {
    convergent_genes <- append(convergent_genes, list(c(i, i + 1)))
  }
}

convergent_genes <- gff[unlist(convergent_genes), ]



divergent_genes<- list()
for (i in 1:(nrow(gff) - 1)) {
  if (gff$strand[i] == "-" & gff$strand[i + 1] == "+") {
    divergent_genes <- append(divergent_genes, list(c(i, i + 1)))
  }
}
divergent_genes <- gff[unlist(divergent_genes), ]



### Get plus and minus genes of interest
genes_of_interest <- rbind(convergent_genes, divergent_genes)
#genes_of_interest <- na.omit(genes_of_interest)
plus_genes <- unique(genes_of_interest$locus_tag[genes_of_interest$strand == "+"])
minus_genes <- unique(genes_of_interest$locus_tag[genes_of_interest$strand == "-"])

num_cores <- N.threads
##### Coverage #####
depth_plus <- read.table(paste0("output/coverage_data/", sample, "_plus_depth.txt"),
                         sep = "\t", header = FALSE)




depth_minus <- read.table(paste0("output/coverage_data/", sample, "_minus_depth.txt"),
                          sep = "\t", header = FALSE)


depth_plus <- depth_plus[!duplicated(depth_plus), ]
depth_minus <- depth_minus[!duplicated(depth_minus), ]

random_genes <- gff %>% 
  filter(strand == "+") %>% 
  sample_n(50)


get_coverage <- function(gene_id){
  
  
  start <- gff[gff$locus_tag == gene_id, ]$start
  end <- gff[gff$locus_tag == gene_id, ]$end
  median_plus <- median(depth_plus[start:end, 3])
  median_minus <- median(depth_minus[start:end, 3])
  
  return(list(median_plus = median_plus, median_minus = median_minus))
  
}

coverage_results <- lapply(random_genes$locus_tag, get_coverage)

# Extract only the median_plus values
median_plus_values <- sum(sapply(coverage_results, function(x) x$median_plus))
median_minus_values <- sum(sapply(coverage_results, function(x) x$median_minus))

if (median_plus_values < median_minus_values) {
  # Swap the values of depth_plus and depth_minus
  temp <- depth_plus
  depth_plus <- depth_minus
  depth_minus <- temp
}


#### Transcript annotation


# Load required libraries

# Clean up any existing connections
if (exists("cl")) {
  try(stopCluster(cl), silent = TRUE)
  rm(cl)
}
closeAllConnections()

# Define the main processing function
calculate_transcript_info <- function(gene_id, depth_table) {
  locus_tag <- gff[gff$locus_tag == gene_id, ]$locus_tag
  gene_name <- gff[gff$locus_tag == gene_id, ]$gene
  strand <- gff[gff$locus_tag == gene_id, ]$strand
  gene_tss <- gff[gff$locus_tag == gene_id, ]$start
  gene_tts <- gff[gff$locus_tag == gene_id, ]$end
  
  # Check if coordinates exist in depth_table
  if (gene_tss > nrow(depth_table) || gene_tts > nrow(depth_table)) {
    return(NULL)
  }
  
  # Extract coverage for the gene region
  coverage <- depth_table[gene_tss:gene_tts, "V3"]
  
  # Check for NA values
  na_percentage <- sum(is.na(coverage)) / length(coverage) * 100
  if(na_percentage > 30) {
    return(NULL)
  }
  
  coverage <- coverage[!is.na(coverage)]
  
  # Calculate gene coverage - only consider non-zero values
  gene_coverage <- if(length(coverage[coverage > 0]) > 0) {
    round(median(coverage[coverage > 0]))
  } else {
    0
  }
  print(threshold)
  cov_threshold <- threshold * gene_coverage
  
  # Rest of the function remains the same...
  is_below_threshold1 <- sapply(coverage[1:min(6, length(coverage))], function(x) x < cov_threshold)
  
  if (sum(is_below_threshold1) == min(6, length(coverage)) | gene_coverage == 0) {
    transcript_start <- gene_tss
  } else {
    upstream_pos <- depth_table[1:gene_tss, ]
    upstream_pos <- upstream_pos[rev(seq_len(nrow(upstream_pos))), ]
    
    upstream_pos$below_threshold <- ifelse(upstream_pos$V3 < cov_threshold, 1, 0)
    upstream_pos$consecutive_count <- cumsum(upstream_pos$below_threshold)
    
    split_df <- upstream_pos[upstream_pos$consecutive_count >= 3, ]
    transcript_start <- if(nrow(split_df) > 0) split_df[1, "V2"] else gene_tss
  }
  
  is_below_threshold2 <- sapply(coverage[(length(coverage) - min(5, length(coverage)-1)):length(coverage)], 
                                function(x) x < cov_threshold)
  
  if (sum(is_below_threshold2) == min(6, length(coverage)) | gene_coverage == 0) {
    transcript_end <- gene_tts
  } else {
    downstream_pos <- depth_table[gene_tts:nrow(depth_table), ]
    downstream_pos$below_threshold <- ifelse(downstream_pos$V3 < cov_threshold, 1, 0)
    downstream_pos$consecutive_count <- cumsum(downstream_pos$below_threshold)
    
    split_df <- downstream_pos[downstream_pos$consecutive_count >= 3, ]
    transcript_end <- if(nrow(split_df) > 0) split_df[1, "V2"] else gene_tts
  }
  
  return(data.frame(
    locus_tag = locus_tag,
    gene_id = gene_name,
    transcript_start = transcript_start,
    transcript_end = transcript_end,
    gene_coverage = gene_coverage
  ))
}


cl <- makeCluster(num_cores)
registerDoParallel(cl)

print("Annotating TUs of plus genes")
# Run the parallel computation
results_list_plus <- foreach(gene_id = plus_genes, .combine = "rbind", .packages = c("dplyr")) %dopar% {
  calculate_transcript_info(gene_id, depth_plus)
}
print("Annotating TUs of minus genes")

results_list_minus <- foreach(gene_id = minus_genes, .combine = "rbind", .packages = c("dplyr")) %dopar% {
  calculate_transcript_info(gene_id, depth_minus)
}

# Stop parallel backend
stopCluster(cl)

transcripts_df <- rbind(results_list_plus, results_list_minus)


######## OVERLAP ANNOTATION #######


df_convergent <- convergent_genes %>%
  left_join(transcripts_df %>%
              select(transcript_start, transcript_end, gene_coverage, locus_tag),
            by = "locus_tag") %>%
  distinct() %>%
  mutate(Type = "",
         Overlapping_cov_plus = "",
         Overlapping_cov_minus = "")


# Modified loop with error handling
for (i in seq(1, nrow(df_convergent) - 1, 2)) {
  tryCatch({
    transcript_end <- df_convergent[i, ]$transcript_end
    transcript_start <- df_convergent[i + 1, ]$transcript_start
    
    if (is.na(transcript_end) || is.na(transcript_start)) {
      next
    }
    
    if (!is.numeric(transcript_end) || !is.numeric(transcript_start)) {
      next
    }
    
    if ((transcript_end - transcript_start) > 20) {
      df_convergent[i, ]$Type <- "Excludon"
      df_convergent[i + 1, ]$Type <- "Excludon"
      range2 <- transcript_end
      range1 <- transcript_start
      
      df_convergent[i, ]$Overlapping_cov_plus <- round(mean(depth_plus[range1:range2, 3]))
      df_convergent[i, ]$Overlapping_cov_minus <- round(mean(depth_minus[range1:range2, 3]))
    }
  }, error = function(e) {
    
  })
}

df_convergent_Excludons <- df_convergent[df_convergent$Type == "Excludon", ]

###############################################
###############################################

df_divergent <- divergent_genes %>%
  left_join(transcripts_df %>%
              select(transcript_start, transcript_end, gene_coverage, locus_tag),
            by = "locus_tag") %>%
  distinct() %>%
  mutate(Type = "",
         Overlapping_cov_plus = "",
         Overlapping_cov_minus = "")



# Modified loop with error handling for divergent genes
for (i in seq(1, nrow(df_divergent) - 1, 2)) {
  tryCatch({
    transcript_end <- df_divergent[i, ]$transcript_end
    transcript_start <- df_divergent[i + 1, ]$transcript_start
    
    if (is.na(transcript_end) || is.na(transcript_start)) {
      next
    }
    
    if (!is.numeric(transcript_end) || !is.numeric(transcript_start)) {
      
      next
    }
    
    if ((transcript_end - transcript_start) > 20) {
      df_divergent[i, ]$Type <- "Excludon"
      df_divergent[i + 1, ]$Type <- "Excludon"
      range2 <- transcript_end
      range1 <- transcript_start
      
      df_divergent[i, ]$Overlapping_cov_plus <- round(mean(depth_plus[range1:range2, 3]))
      df_divergent[i, ]$Overlapping_cov_minus <- round(mean(depth_minus[range1:range2, 3]))
    }
  }, error = function(e) {
    
  })
}

df_divergent_Excludons <- df_divergent[df_divergent$Type == "Excludon", ]


message("Number of Excludons:\n",
        "  Convergent: ", nrow(df_convergent_Excludons)/2, "\n",
        "  Divergent: ", nrow(df_divergent_Excludons)/2)


if (!dir.exists("output/Excludons")) {
  dir.create("output/Excludons")
}

write.csv(df_convergent_Excludons, file = paste0("output/Excludons/Convergent_Excludons_", sample, ".csv"))
write.csv(df_divergent_Excludons, file = paste0("output/Excludons/Divergent_Excludons_", sample, ".csv"))




