# Loading required packages
suppressMessages(library(ShortRead))  
library(tidyverse)   
suppressMessages(library(yaml))
suppressMessages(library(gridExtra))

# Customized functions----
# source("Static/R/common.R")

resultDir <- file.path("Analysis", "Results")
dir.create(resultDir, showWarnings = FALSE, recursive = TRUE) 

resultDir1 <- file.path("Results", "QC")
dir.create(resultDir1, showWarnings = FALSE, recursive = TRUE) 

# Loading the config file----
config <- yaml::yaml.load_file("config.yaml")

# Create a study design by using the config.yaml-file
studyDesign <- data.frame()
for (i in 1:length(config$Samples)) {
  studyDesign <- rbind(studyDesign,
                       data.frame(filename=unlist(config$Samples[i])))
}
coverage_filter <- function(df) {
  df %>%
    mutate(coverage_filter = case_when(
      VariantLevel < 0.02 ~ "cov_out2pct",
      VariantLevel <= 0.05 & VariantLevel >=0.02 & Coverage < 33275 ~ "cov_out33275",
      VariantLevel <= 0.2 & VariantLevel >=0.05 & Coverage < 6386 ~ "cov_out6386",
      VariantLevel <= 0.35 & VariantLevel >0.2 & Coverage < 1344 ~ "cov_out1344",
      VariantLevel <= 0.5 & VariantLevel >0.35 & Coverage < 624 ~ "cov_out624",
      VariantLevel <= 0.65 & VariantLevel >0.5 & Coverage < 336 ~ "cov_out336",
      VariantLevel <= 0.8 & VariantLevel >0.65 & Coverage < 181 ~ "cov_out181",
      VariantLevel <= 1 & VariantLevel >0.8 & Coverage < 84 ~ "cov_out84",
      TRUE ~ "in"))
}
rCRS <- read_file(paste0(ref_dir,"/reference/rCRS.fasta")) %>% 
  unlist() %>% 
  str_remove("chrM") %>%
  str_remove(">") %>%
  str_remove_all("\n") %>%
  str_remove_all("\r")%>%
  str_remove_all("\r") |>
  str_remove("rCRS")
mitopore_blacklist <- str_locate_all(rCRS, 
                                     "CCCCCC|AAAAAA|TTTTTT|GGGGGG") %>%
  as.data.frame() %>% 
  rowwise() %>%
  mutate(mitopore_blacklist = list(start:end)) %>%
  unnest(mitopore_blacklist) %>% 
  pull(mitopore_blacklist) %>%
  c(3107)
vcf_file <- read_delim("Analysis/Results/result1.txt") %>% 
  coverage_filter() %>% 
  dplyr::filter(!Pos %in% mitopore_blacklist, !Filter == "STRAND_BIAS") %>%
  dplyr::filter(coverage_filter == "in")

disease_df <- readxl::read_xlsx(paste0(ref_dir,"/static/Lists MITOMAP/MutationsCodingControl.xlsx"), skip = 1) %>% 
  dplyr::select(Position, `Plasmy Reports(Homo/Hetero)`, NucleotideChange, Disease) %>%
  separate(NucleotideChange, into = c("REF", "ALT"), sep = "-") %>%
  separate(`Plasmy Reports(Homo/Hetero)`, into = c("Homoplasmy", "Heteroplasmy"), sep = "/") %>%
  mutate(POS = Position, DiseaseStatus = Disease) %>%
  dplyr::select(-Position, -Disease)

rtRNA_df <- readxl::read_xlsx(paste0(ref_dir,"/static/Lists MITOMAP/MutationsRNA MITOMAP Foswiki.xlsx"), skip = 1) %>% 
  dplyr::select(Position, Allele, Disease, Status, Homoplasmy, Heteroplasmy) %>%
  separate(Allele, into = c("REF", "ALT"), sep = "\\d+") %>%
  mutate(POS = Position, DiseaseStatus_tr = Disease) %>%
  dplyr::select(-Position, -Disease)

poly_df <- readxl::read_xlsx(paste0(ref_dir,"/static/Lists MITOMAP/Polymorphisms MITOMAP Foswiki.xlsx")) %>% 
  dplyr::select(pos, ref, alt, aachange) %>%
  mutate(POS = pos, REF = ref, ALT = alt, `Status` = "polymorphism") %>%
  dplyr::select(-pos, -ref, -alt)

orgdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

check <- ensembldb::genes(orgdb, ensembldb::listColumns(orgdb, "is_circular"))
mt_pos <- as_tibble(check) %>% dplyr::filter(seqnames == "MT")
mt_genome <- ensembldb::genes(orgdb, filter = AnnotationFilter::SeqNameFilter("MT")) %>% as_tibble()
mt_anno <- mt_genome %>%
  rowwise() %>%
  mutate(POS = list(start:end)) %>%
  unnest(POS) %>%
  dplyr::select(POS, gene_name, gene_biotype)



rect_positions <- mt_genome %>% distinct(gene_name, start, end, gene_biotype) %>% mutate(gene_biotype = case_when(is.na(gene_biotype) ~ "unknown",
                                                                                                                  str_detect(gene_biotype, "tRNA") ~ "tRNA",
                                                                                                                  str_detect(gene_biotype, "rRNA") ~ "rRNA",
                                                                                                                  TRUE ~ gene_biotype),
                                                                                         start = as.double(start),
                                                                                         end = as.double(end))
vcf_file <- vcf_file  %>%
  full_join(mt_anno, by=c("Pos"="POS")) %>%
  full_join(poly_df, by=c("Pos"="POS", "Ref"="REF", "Variant"="ALT")) %>%
  full_join(rtRNA_df, by=c("Pos"="POS", "Ref"="REF", "Variant"="ALT")) %>%
  full_join(disease_df, by=c("Pos"="POS", "Ref"="REF", "Variant"="ALT")) %>%
  mutate(Status = case_when(is.na(Status.x) ~ "no polymorphism", TRUE ~ Status.x)) %>%
  mutate(trans = factor(str_c(Ref, ">", Variant)))

writexl::write_xlsx(vcf_file %>% dplyr::filter(!is.na(ID)), "Analysis/Results/disease_table.xlsx")

vcf_file %>%
  dplyr::filter(!is.na(ID)) %>%
  ggplot(data = ., mapping = aes(Pos, VariantLevel)) +
  geom_point(aes(color = Status), size = 6, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill=alpha("orange", 0.5)),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.text = element_text(size = 26, face = "bold"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 45, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 16)) +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.02, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
  geom_rect(data = filter(rect_positions, gene_biotype == "protein_coding"),
            aes(xmin = start, xmax = end, ymin = -.2, ymax = -.10),
            fill = "cornflowerblue", color = "black", inherit.aes = F) +
  geom_text(data = filter(rect_positions, gene_biotype == "protein_coding"),
            aes(x = (start+end)/2, y = -.3, label = str_remove(gene_name, "MT-")), inherit.aes = F, angle = 90, size = 7, fontface="italic") +
  geom_rect(data = filter(rect_positions, !gene_biotype == "protein_coding"),
            aes(xmin = start, xmax = end, ymin = -.175, ymax = -.125),
            fill = "red", color = "black", inherit.aes = F) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), limits = c(-0.35,1)) +
  scale_x_continuous(limits = c(0, 16569), breaks = (c(0, 5000, 10000, 15000, 16569))) +
  scale_color_manual(values = c("polymorphism" = "limegreen", "no polymorphism" = "lightcoral")) +
  #scale_shape_manual(values = c("indel" = 4, "SNV" = 20)) +
  ylab("Allele Frequency\n") +
  xlab("\nnt position in mt Genome") +
  ggrepel::geom_text_repel(data = vcf_file %>%
                             dplyr::filter(!Status == "polymorphism", VariantLevel >= 0.05, !gene_name == "non-coding") %>%
                             distinct(), aes(label = gene_name), size = 6, fontface="italic") +
  if(vcf_file %>% pull(ID) %>% unique() %>% length() > 1){facet_wrap(~ID, ncol=1)}
  
nr_samples <- vcf_file %>%
  pull(ID) %>%
  unique() %>%
  length()

ggsave('Results/mitomap.png',
       height = if(nr_samples > 1) {
         if(nr_samples>5) {
           45
         } else{nr_samples*9}
       } else{9}, width = 19)



# Raw sequence review----
processQCFastq <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file, widthIds=FALSE)
  c(
    reads = formatC(length(fastq), big.mark=","),
    mbs = formatC(round(sum(width(fastq)) / 1000 / 1000, digits=1), big.mark=","),
    min = min(width(fastq)),
    max = max(width(fastq)),
    mean = round(mean(width(fastq)), digits=1),
    median = round(median(width(fastq)), digits=0),
    qval = round(mean(alphabetScore(fastq) / width(fastq)), digits=1),
    gc = round(mean(letterFrequency(sread(fastq), "GC")  / width(fastq)) * 100, digits=1)
  )
}

data <- lapply(row.names(studyDesign), processQCFastq)
qcData <- data.frame(data)
colnames(qcData) <- row.names(studyDesign)
qcData$INFO <-row.names(qcData)
qcData <- qcData[,c(ncol(qcData),1:(ncol(qcData)-1))]
df <- dplyr::tibble(qcData)
df_t <- gt::gt(df)
gt::gtsave(df_t, "Results/QC/sample_summary.html")

# png("Results/QC/sample_summary.png",  width=600, height = 300)       # Export png
# grid.table(qcData)
# dev.off()


flagstatTargets <- file.path("Analysis", "flagstat", 
                             paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE),".txt",sep=""))

loadFlagstat <- function(file) {
  x <- read.table(file, header=FALSE, sep=" ", fill=NA)
  x[1:8,1]
}

flagstatRes <- data.frame(matrix(unlist(lapply(flagstatTargets, loadFlagstat)), ncol=length(flagstatTargets)), stringsAsFactors = FALSE)
colnames(flagstatRes) <- rownames(studyDesign)
rownames(flagstatRes) <- c("Total read/QC passed", "Primary","Secondary", "Supplementary", "Duplicates", "Primary Duplicated","Total Mapped", "Primary Mapped")


# flagstatRes[nrow(flagstatRes)+1,] <- as.numeric(gsub(",","",t(qcData)[, "reads"]))
# rownames(flagstatRes)[6] <- "nreads"
# 
# getVal <- function(word) {
#   sum(as.numeric(unlist(strsplit(word, "/"))))
# }
# 
# zreads <- unlist(lapply(flagstatRes["read mappings", ], getVal)) -
#   unlist(lapply(flagstatRes["Secondary", ], getVal)) -
#   unlist(lapply(flagstatRes["Supplementary", ], getVal)) - 
#   unlist(lapply(flagstatRes["Duplicates", ], getVal)) 

flagstatRes[nrow(flagstatRes)+1,] <- round(as.numeric(flagstatRes["Total Mapped", ]) / as.numeric(flagstatRes["Total read/QC passed", ]) * 100, digits = 2)
rownames(flagstatRes)[9] <- "%mapping"

flagstatRes$INFO <-row.names(flagstatRes)
flagstatRes <- flagstatRes[,c(ncol(flagstatRes),1:(ncol(flagstatRes)-1))]

df <- dplyr::tibble(flagstatRes)
df_t <- gt::gt(df)
gt::gtsave(df_t, "Results/QC/mtDNA_mapping_summary.html")


# png("Results/QC/mtDNA_mapping_summary.png",  width=600, height = 300)     # Export PDF
# grid.table(flagstatRes)
# dev.off()



# Distribution of read Lengths (bp)----
extractLengths <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file)
  t(as.matrix(width(fastq)))
}
lengthData <- mapply(extractLengths, row.names(studyDesign))

if (length(studyDesign$filename)>1){
  lengthDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(lengthData)))
  colnames(lengthDataMatrix) <-  row.names(studyDesign)
  lengthMatrixMelt <- reshape2::melt(lengthDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))
  lengthMatrixMelt <- cbind(lengthMatrixMelt, group=studyDesign[match(lengthMatrixMelt$variable,  rownames(studyDesign)), "filename"])
  plot <- ggplot(lengthMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(0, as.numeric(quantile(lengthMatrixMelt$value, probs=c(0.975))))) + xlab("study sample") +  ylab("Distribution of Read Lengths (bp)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read lengths across samples")
  ggsave("Results/QC/read_length.png", width=12, height = 4)
} else{
  lengthDataMatrix <- data.frame(lengthData)
  colnames(lengthDataMatrix) <-  row.names(studyDesign)
  lengthMatrixMelt <- reshape2::melt(lengthDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))
  plot <- ggplot(lengthMatrixMelt, aes(x=variable, y=value)) + geom_violin(fill="cornflowerblue") + scale_y_continuous(limits=c(0, as.numeric(quantile(lengthMatrixMelt$value, probs=c(0.975))))) + xlab("study sample") +  ylab("Distribution of Read Lengths (bp)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution \nof read lengths across samples")
  ggsave("Results/QC/read_length.png", width=6, height = 4)
}

# read quality----
extractQualities <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file)
  t(as.matrix(alphabetScore(fastq) / width(fastq)))
}
qualityData <- mapply(extractQualities, row.names(studyDesign))
if (length(studyDesign$filename)>1){
  qualityDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(qualityData)))
  colnames(qualityDataMatrix) <-  row.names(studyDesign)
  
  qualityMatrixMelt <- reshape2::melt(qualityDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))
  qualityMatrixMelt <- cbind(qualityMatrixMelt, group=studyDesign[match(qualityMatrixMelt$variable,  rownames(studyDesign)), "filename"])
  
  plotQ <- ggplot(qualityMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(min(qualityMatrixMelt$value), max(qualityMatrixMelt$value))) + xlab("study sample") +  ylab("Distribution of Read Qualities (QV)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read qualities across samples")
  ggsave("Results/QC/read_quality.png", width=12, height = 4)
} else {
  #TODO
  qualityDataMatrix <- data.frame(qualityData)
  colnames(qualityDataMatrix) <-  row.names(studyDesign)
  
  qualityMatrixMelt <- reshape2::melt(qualityDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))
  qualityMatrixMelt <- cbind(qualityMatrixMelt, group=studyDesign[match(qualityMatrixMelt$variable,  rownames(studyDesign)), "filename"])
  
  plotQ <- ggplot(qualityMatrixMelt, aes(x=variable, y=value)) + geom_violin(fill="cornflowerblue") + scale_y_continuous(limits=c(min(qualityMatrixMelt$value), max(qualityMatrixMelt$value))) + xlab("study sample") +  ylab("Distribution of Read Qualities (QV)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution \nof read qualities across samples")
  ggsave("Results/QC/read_quality.png", width=6, height = 4)
}

#haplogroup and large-scale deletions
read_output <- function(input) {
  read_delim(input) %>%
    coverage_filter() %>%
    dplyr::filter(coverage_filter == "in",
                  !Pos %in% mitopore_blacklist
                  )
}

write_hsd <- function(input, output, variant_level = 0.05, input_type = "file") {
  if(input_type == "file"){
    df <- read_output(input) %>%
      dplyr::filter(VariantLevel>=variant_level,
                    !Filter == "BLACKLISTED") %>%
      group_by(ID) %>%
      mutate(ID = ID, Range = "1-16569", Haplogroup = "?", Polymorphisms = list(str_c(Pos, Variant))) %>%
      ungroup() %>%
      dplyr::select(ID, Range, Haplogroup, Polymorphisms) %>%
      distinct() %>%
      mutate(Polymorphisms = unlist(lapply(Polymorphisms, paste, collapse = "\t"))) %>%
      write_delim(output, delim = "\t ", quote = "none")
    
  } else if (input_type == "df"){
    df <- df %>%
      dplyr::filter(VariantLevel>=variant_level,
                    !Filter == "BLACKLISTED") %>%
      group_by(ID) %>%
      mutate(ID = ID, Range = "1-16569", Haplogroup = "?", Polymorphisms = list(str_c(Pos, Variant))) %>%
      ungroup()
    dplyr::select(ID, Range, Haplogroup, Polymorphisms) %>%
      distinct() %>%
      mutate(Polymorphisms = unlist(lapply(Polymorphisms, paste, collapse = "\t"))) %>%
      write_delim(output, delim = "\t ", quote = "none")
  }
  
}

write_hsd("Analysis/Results/result1.txt", "haplogrep_hsd.hsd", variant_level = 0.75, input_type = "file")
#run haplogrep 3 and visualize file (test linux version, because windows version is not up to date)

#Detect large deletions -- *.coverage file necessary (make with samtools depth command->get from coverage folder)
extract_consecutive_ranges <- function(vec) {
  split(vec, cumsum(c(TRUE, diff(vec) != 1)))
}
extract_deletions <- function(filename, min = 500) {
  df <- read_delim(filename, col_names = F, show_col_types = F)
  poi <- df %>% dplyr::filter(X3<min)
  pois <- extract_consecutive_ranges(poi %>% pull(X2))
  for(i in 1:length(pois)){
    pois[[i]] <- tibble(POS = pois[[i]], range_length = length(unlist(pois[[i]])), id=as.character(i))
  }
  pois <- pois %>% purrr::reduce(full_join, by = join_by(POS, range_length, id)) %>% dplyr::filter(range_length > 200) 
  if(nrow(pois) == 0){pois <- tibble()} else{
    pois <- pois %>%
      group_by(id) %>%
      summarise(from = min(POS), to = max(POS), len = max(POS)-min(POS))
    pois <- pois %>% dplyr::select(-1)
  }
}
extract_deletions_shift <- function(filename, min = 500) {
  df <- read_delim(filename, col_names = F, show_col_types = F) %>% 
    mutate(X4 = case_when(X2<8285 ~ X2+8285, TRUE ~ X2-8284))
  poi <- df %>% dplyr::filter(X3<min)
  pois <- extract_consecutive_ranges(poi %>% pull(X4))
  for(i in 1:length(pois)){
    pois[[i]] <- tibble(POS = pois[[i]], range_length = length(unlist(pois[[i]])), id=as.character(i))
  }
  pois <- pois %>% purrr::reduce(full_join, by = join_by(POS, range_length, id)) %>% dplyr::filter(range_length > 200) 
  if(nrow(pois) == 0){pois <- tibble()} else{
    pois <- pois %>%
      group_by(id) %>%
      summarise(from = min(POS), to = max(POS), len = max(POS)-min(POS)) %>%
      mutate(from = case_when(
        from > 8286 ~ from - 8285, TRUE ~ from + 8284
      ), 
      to = case_when(
        to > 8286 ~ to - 8285, TRUE ~ to + 8284
      )
      ) %>%
      dplyr::filter(!from > 5000, !to < 15000)
    pois <- pois %>% dplyr::select(-1)
  }
}

#make deletion table (final result)
make_deletion_tibble <- function(filename, min = 500) {
  if(nrow(extract_deletions(filename, min))
     +
     nrow(extract_deletions(filename, min)) == 0){tibble()} else{
       full_join(extract_deletions(filename, min), 
                 extract_deletions_shift(filename, min))
     }
}
deletion_tibbles <- list()
for(files in list.files("Analysis/depths/")) {
  deletion_tibbles[[files]] <- make_deletion_tibble(filename = str_c("Analysis/depths/", 
                                                                     list.files("Analysis/depths/")[str_detect(list.files("Analysis/depths/"), files)]), 
                                                                     300)
}
deletions <- character()
for(i in 1:length(deletion_tibbles)) {
  if(nrow(deletion_tibbles[[i]]) == 0) {
    deletions <- c(deletions, 
                       str_c("No large-scale deletions were detected in the sample ",
                             str_remove(list.files("Analysis/depths")[i], ".txt")))
  } else if(nrow(deletion_tibbles[[i]]) > 0) 
    {deletions <- c(deletions, str_c("Attention! Your data suggests the presence of ", 
                       nrow(deletion_tibbles[[i]]), " large deletion/s in the sample ", 
                       str_remove(list.files("Analysis/depths")[i], ".txt"),
                       "! Please check your sample carefully to verify!"))
}
}



cov_dfs <- str_c("Analysis/depths/", list.files("Analysis/depths/")) %>% 
  set_names() %>% 
  map(.f = read_delim, col_names = F, show_col_types = F)

for(i in 1:length(cov_dfs)) {
  cov_dfs[[i]] <- cov_dfs[[i]] %>% mutate(id = str_remove_all(names(cov_dfs), "Analysis/depths/") %>% str_remove_all(".txt") %>% .[i])
}
cov_dfs <- cov_dfs %>% reduce(full_join, by=c("X1", "X2", "X3", "id"))
pdf("Analysis/Results/deletion_plots.pdf", width = 5, height = 3*length(unique(cov_dfs$id)))
cov_dfs %>%
  ggplot(aes(X2, log10(X3))) +
  geom_point(alpha = 0.5, color = "cornflowerblue") +
  facet_wrap(~id, scales="free", ncol=1)  +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill=scales::alpha("orange", 0.5)),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12, face = "italic"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom") +
  geom_hline(yintercept=log10(300), color="red", linewidth = 1.1, linetype = "dashed") +
  ylab("log 10 Coverage") +
  xlab("Position in mtDNA")
dev.off()


haplo_matrix <- read_delim("Analysis/Results/result1.txt") %>% 
  coverage_filter() %>% 
  dplyr::filter(coverage_filter == "in", VariantLevel > 0.75) %>% 
  dplyr::select(ID, Pos, VariantLevel) %>% 
  pivot_wider(names_from = "ID", values_from = "VariantLevel", values_fn = {mean}, values_fill = 0) %>% 
  dplyr::select(-1) %>% 
  as.matrix() %>% 
  scale

pdf("Analysis/Results/variant_landscape.pdf", width = 5, height = 7)
plot_varscape <- function(){
  if(length(haplo_matrix) > 0) {
    ComplexHeatmap::Heatmap(t(haplo_matrix), row_labels = str_remove_all(rownames(t(haplo_matrix)), ".bam"))
  } else {
    text(x=0.5,y=0.5,"No variants were detected.\nVariance landscape plot cant be displayed!", font=2)
}}
plot_varscape()
dev.off()


#Potentially detrimental mutations
vcf_file %>%
  dplyr::filter(!is.na(ID), is.na(Status.x)) %>%
  group_by(ID) %>%
  dplyr::count()


library(stringr)
library(vcfR)
library(tidyverse)
process_vcf <- function(vcf_file, format = c("GT:AD:AF:DP:F1R2:F2R1:FAD:SB")) {

  vcf_file <- read.vcfR(vcf_file)

  format <- str_split(format, ":", simplify = TRUE)

  cbind(as_tibble(vcf_file@fix), as_tibble(vcf_file@gt)) %>%

    dplyr::select(-c(CHROM, ID, QUAL, FILTER)) %>%

    pivot_longer(cols=c(6:ncol(.)), names_to = "SAMPLE", values_to = "DATA") %>%

    dplyr::select(-INFO) %>%

    separate(DATA, into = c(format), sep = ":", extra = "merge") %>%

    separate(AF, c("GT1", "GT2", "GT3", "GT4", "GT5"), sep = ",", extra = "drop") %>%

    pivot_longer(cols = c(GT1:GT5), names_to = "GType", values_to = "AF") %>%

    mutate(AF = as.double(AF), DP = as.double(DP), POS = as.double(POS)) %>%

    dplyr::select(-FORMAT, -AD, - F1R2, -F2R1, -FAD, -SB, -GT) %>%

    separate(ALT, "ALT", sep = ",", extra = "drop")

}

#filter according to stringent criteria function
#filter according to stringent criteria function

filter_indels <- function(vcf_file,

                          depth = 200,

                          min_mut_rate = 0.2,

                          max_mut_rate = 0.8) {

  vcf_file %>%

    dplyr::filter(DP>depth,

                  AF>min_mut_rate,

                  AF<max_mut_rate,

                  !POS %in% mitopore_blacklist,

                  POS > 580,

                  POS < 16024,

                  !nchar(REF)==2,

                  !nchar(REF) == nchar(ALT),

                  GType == "GT1",

                  !(nchar(REF) - nchar(ALT)) == 1)

}

if (dir.exists("Analysis/INDEL")){
  for (vcf_pro_file in list.files("Analysis/INDEL/")){
    s1<-filter_indels(process_vcf(paste0("Analysis/INDEL/", vcf_pro_file)))
    write.table(s1, paste0("Analysis/INDEL/", vcf_pro_file),sep = "\t", row.names = TRUE, col.names = TRUE)
  }
}

write_delim(as_tibble(deletions), "Analysis/Results/large_deletions.txt", col_names = F)
