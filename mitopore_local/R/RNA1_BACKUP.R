# Loading required packages
suppressMessages(library(digest))
suppressMessages(library(ShortRead))  
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))      
suppressMessages(library(tidyr))
suppressMessages(library(pcaMethods))
suppressMessages(library(kableExtra))  
suppressMessages(library(caTools))     
suppressMessages(library(yaml))
#suppressMessages(library(session))
suppressMessages(library(RColorBrewer))
suppressMessages(library(PoiClaClu))
suppressMessages(library(karyoploteR))
suppressMessages(library(viridis))
suppressMessages(library(gridExtra))
# Custamized functions----
# source("Static/R/common.R")

resultDir <- file.path("Analysis", "Results")
dir.create(resultDir, showWarnings = FALSE, recursive = TRUE) 

# Loading the config file----
config <- yaml::yaml.load_file("config.yaml")

# Create a study design by using the config.yaml-file
studyDesign <- data.frame()
for (i in 1:length(config$Samples)) {
  studyDesign <- rbind(studyDesign,
                       data.frame(filename=unlist(config$Samples[i])))
}



library(vcfR)
library(tidyverse)
if(!interactive()) pdf(NULL)
#vcf_file <- read.vcfR("path_to_vcf")
vcf_file <- read.vcfR("Analysis/Results/result1.vcf")
vcf_file <- if(cbind(vcf_file@fix %>% as_tibble(), vcf_file@gt %>% as_tibble()) %>%
               dplyr::select(-c(CHROM, ID, QUAL, FILTER)) %>%
               ncol() > 6) {cbind(vcf_file@fix %>% as_tibble(), vcf_file@gt %>% as_tibble()) %>%
    dplyr::select(-c(CHROM, ID, QUAL, FILTER)) %>% 
    pivot_longer(cols=c(6:ncol(.)), names_to="SAMPLE", values_to="DATA") %>%
    separate(INFO, into = c("AC", "AN"), sep = ";") %>%
    mutate(AN = str_remove(AN, "..="), AC = str_remove(AC, "..=")) %>%
    separate(DATA, into = c("GT", "AF", "DP"), sep = ":") %>% dplyr::select(-FORMAT) %>%
    mutate(SAMPLE = str_remove(SAMPLE, ".bam")) %>%
    mutate(AF = case_when(
      AF == "." ~ 0,
      TRUE ~ as.numeric(AF)
    )) %>%
    dplyr::filter(!is.na(DP))} else {
      cbind(vcf_file@fix %>% as_tibble(), vcf_file@gt %>% as_tibble()) %>%
        dplyr::select(-c(CHROM, ID, QUAL, FILTER)) %>%
        separate(INFO, into = c("AC", "AN"), sep = ";") %>%
        mutate(AN = str_remove(AN, "..="), AC = str_remove(AC, "..=")) %>%
        mutate(SAMPLE = str_remove(colnames(.)[max(ncol(.))], ".bam")) %>%
        dplyr::select(-FORMAT) %>%
        separate(colnames(.)[max(ncol(.))-1], into = c("GT", "AF", "DP"), sep = ":") %>%
        mutate(AF = case_when(
          AF == "." ~ 0,
          TRUE ~ as.numeric(AF)
        )) %>%
        dplyr::filter(!is.na(DP))
    }
disease_df <- readxl::read_xlsx("/home/ag-rossi/projects/mitopore/mitopore/static/Lists MITOMAP/MutationsCodingControl MITOMAP Foswiki.xlsx", skip = 1) %>% 
  dplyr::select(Position, `Plasmy Reports(Homo/Hetero)`, NucleotideChange, Disease) %>%
  separate(NucleotideChange, into = c("REF", "ALT"), sep = "-") %>%
  separate(`Plasmy Reports(Homo/Hetero)`, into = c("Homoplasmy", "Heteroplasmy"), sep = "/") %>%
  mutate(POS = Position, DiseaseStatus = Disease) %>%
  dplyr::select(-Position, -Disease)

rtRNA_df <- readxl::read_xlsx("/home/ag-rossi/projects/mitopore/mitopore/static/Lists MITOMAP/MutationsRNA MITOMAP Foswiki.xlsx", skip = 1) %>% 
  dplyr::select(Position, Allele, Disease, Status, Homoplasmy, Heteroplasmy) %>%
  separate(Allele, into = c("REF", "ALT"), sep = "\\d+") %>%
  mutate(POS = Position, DiseaseStatus_tr = Disease) %>%
  dplyr::select(-Position, -Disease)

poly_df <- readxl::read_xlsx("/home/ag-rossi/projects/mitopore/mitopore/static/Lists MITOMAP/Polymorphisms MITOMAP Foswiki.xlsx") %>% 
  dplyr::select(pos, ref, alt, aachange) %>%
  mutate(POS = pos, REF = ref, ALT = alt, `Status` = "polymorphism") %>%
  dplyr::select(-pos, -ref, -alt)

orgdb <- EnsDb.Hsapiens.v104::EnsDb.Hsapiens.v104
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

vcf_file <- vcf_file %>% mutate(POS = as.double(POS)) %>%
  full_join(mt_anno) %>%
  mutate(gene_name = case_when(is.na(gene_name) ~ "non-coding", TRUE ~ gene_name),
         gene_biotype = case_when(is.na(gene_name) ~ "non-coding", TRUE ~ gene_name)) %>%
  left_join(poly_df) %>% left_join(disease_df) %>% left_join(rtRNA_df) %>%
  mutate(Status = case_when(
    is.na(Status) ~ "Unknown",
    TRUE ~ Status)) %>%
  mutate(gene_name = case_when(is.na(gene_name) ~ "non-coding",
                               TRUE ~ gene_name),
         gene_biotype = case_when(is.na(gene_biotype) ~ "unknown",
                                  str_detect(gene_biotype, "tRNA") ~ "tRNA",
                                  str_detect(gene_biotype, "rRNA") ~ "rRNA",
                                  TRUE ~ gene_biotype)) %>%
  mutate(Status = case_when(
    !is.na(DiseaseStatus) ~ DiseaseStatus,
    !is.na(DiseaseStatus_tr) ~ DiseaseStatus_tr,
    !is.na(Status) ~ Status,
    TRUE ~ "uda"
  )) %>%
  mutate(snv_type = case_when(
    !Status %in% c("polymorphism", "uda") ~ "pda",
    Status == "polymorphism" ~ "nda",
    TRUE ~ Status)) %>%
  mutate(
    aachangetype = case_when(
      !aachange %in% c("non-coding", "tRNA", "rRNA") ~ ifelse(str_extract(.$aachange, "^[A-Z]") == str_extract(.$aachange, "[A-Z]$"), "synonymous", "non-synonymous"),
      TRUE ~ aachange
    )) %>%
  mutate(aachangetype = case_when(is.na(aachangetype) ~ "unknown", TRUE ~ aachangetype)) %>%
  mutate(type = case_when(
    str_detect(ALT, "[A-Z]{2,}") & !str_detect(REF, "[A-Z]{2,}") ~ "indel",
    TRUE ~ "SNV"
  ))

options(repr.plot.width=16, repr.plot.height=if(vcf_file %>% pull(SAMPLE) %>% unique() %>% length() > 1){18} else{9})
vcf_file %>%
  dplyr::filter(!is.na(SAMPLE)) %>%
  mutate(wrap = case_when(AF>0.75 ~ 1, TRUE ~2)) %>%
  mutate(Status = case_when(!Status == "polymorphism" ~ "no polymorphism", TRUE ~ Status)) %>%
  mutate(trans = factor(str_c(REF, ">", ALT))) %>%
  ggplot(data = ., mapping = aes(POS, AF)) +
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
                             mutate(Status = case_when(!Status == "polymorphism" ~ "no polymorphism", TRUE ~ Status)) %>%
                             dplyr::filter(!Status == "polymorphism", AF >= 0.05, !gene_name == "non-coding") %>%
                             distinct(), aes(label = gene_name), size = 6, fontface="italic") +
  if(vcf_file %>% pull(SAMPLE) %>% unique() %>% length() > 1){facet_wrap(~SAMPLE, ncol=1)}
ggsave('Results/mitomap.png', height = 39, width = 19)


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
    #n50 = ncalc(width(fastq), n=0.5),
    #l50 = lcalc(width(fastq), n=0.5),
    #n90 = ncalc(width(fastq), n=0.9),
    #l90 = lcalc(width(fastq), n=0.9)
  )
}

data <- lapply(row.names(studyDesign), processQCFastq)
qcData <- data.frame(data)
colnames(qcData) <- row.names(studyDesign)
png("Results/QC/sample_summary.png",  width=600, height = 300)       # Export png
grid.table(qcData)
dev.off()

# Distribution of read Lengths (bp)----
extractLengths <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file)
  t(as.matrix(width(fastq)))
}
lengthData <- mapply(extractLengths, row.names(studyDesign))
lengthDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(lengthData)))
colnames(lengthDataMatrix) <-  row.names(studyDesign)

lengthMatrixMelt <- reshape2::melt(lengthDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))
lengthMatrixMelt <- cbind(lengthMatrixMelt, group=studyDesign[match(lengthMatrixMelt$variable,  rownames(studyDesign)), "filename"])

plot <- ggplot(lengthMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(0, as.numeric(quantile(lengthMatrixMelt$value, probs=c(0.975))))) + xlab("study sample") +  ylab("Distribution of Read Lengths (bp)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read lengths across samples")
ggsave("Results/QC/read_length.png", width=12, height = 4)

# read quality----
extractQualities <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file)
  t(as.matrix(alphabetScore(fastq) / width(fastq)))
}
qualityData <- mapply(extractQualities, row.names(studyDesign))
qualityDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(qualityData)))
colnames(qualityDataMatrix) <-  row.names(studyDesign)

qualityMatrixMelt <- reshape2::melt(qualityDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))
qualityMatrixMelt <- cbind(qualityMatrixMelt, group=studyDesign[match(qualityMatrixMelt$variable,  rownames(studyDesign)), "filename"])

plotQ <- ggplot(qualityMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(min(qualityMatrixMelt$value), max(qualityMatrixMelt$value))) + xlab("study sample") +  ylab("Distribution of Read Qualities (QV)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read qualities across samples")
ggsave("Results/QC/read_quality.png", width=12, height = 4)


# Review of cDNA read mapping----
# Read mappings correspond to the number of unique reads. In our worklflow we disabled secondary mappings.
# A Supplementary alignment could correspond to an chimeric read. For instance a result of a read with a barcode located in the middle.
# The informations are obtained from the .txt-Files with the flagstats, created from the .bam-Files.
flagstatTargets <- file.path("Analysis", "flagstat", 
                             paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE),".txt",sep=""))

loadFlagstat <- function(file) {
  x <- read.table(file, header=FALSE, sep=" ", fill=NA)[c(1:5),c(1,3)]
  x[,1]
}

flagstatRes <- data.frame(matrix(unlist(lapply(flagstatTargets, loadFlagstat)), ncol=length(flagstatTargets)), stringsAsFactors = FALSE)
colnames(flagstatRes) <- rownames(studyDesign)
rownames(flagstatRes) <- c("read mappings", "Secondary", "Supplementary", "Duplicates", "Mapped")


flagstatRes[nrow(flagstatRes)+1,] <- as.numeric(gsub(",","",t(qcData)[, "reads"]))
rownames(flagstatRes)[6] <- "nreads"

getVal <- function(word) {
  sum(as.numeric(unlist(strsplit(word, "/"))))
}

zreads <- unlist(lapply(flagstatRes["read mappings", ], getVal)) -
  unlist(lapply(flagstatRes["Secondary", ], getVal)) -
  unlist(lapply(flagstatRes["Supplementary", ], getVal)) - 
  unlist(lapply(flagstatRes["Duplicates", ], getVal)) 

flagstatRes[nrow(flagstatRes)+1,] <- round(as.numeric(flagstatRes["nreads", ]) / as.numeric(flagstatRes["read mappings", ]) * 100, digits = 2)

flagstatRes <- flagstatRes[c(6,1,2,3,4,7),]
rownames(flagstatRes)[6] <- "%mapping"

png("Results/QC/mtDNA_mapping_summary.png",  width=600, height = 300)     # Export PDF
grid.table(flagstatRes)
dev.off()