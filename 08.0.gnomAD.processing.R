library(VariantAnnotation)
library(data.table)
library(dplyr)


process_vcf_with_af <- function(vcf_file) {
  # Read the VCF file
  vcf <- readVcf(vcf_file, "hg38")
  
  # Extract the INFO field
  info <- info(vcf)
  
  # Extract specific AF fields
  af_afr <- info$AF_afr
  af_ami <- info$AF_ami
  af_amr <- info$AF_amr
  af_asj <- info$AF_asj
  af_eas <- info$AF_eas
  af_fin <- info$AF_fin
  af_mid <- info$AF_mid
  af_nfe <- info$AF_nfe
  af_remaining <- info$AF_remaining
  af_sas <- info$AF_sas
  
  # Create a data table with the extracted information
  alt = CharacterList(alt(vcf))
  alt = unstrsplit(alt, sep = ",")
  
  snp_data <- data.frame(
    CHROM = as.character(seqnames(vcf)),
    POS = start(vcf),
    ID = names(vcf),
    REF = ref(vcf),
    ALT = alt,
    AF = as.numeric(info$AF),
    AF_afr = as.numeric(af_afr),
    AF_ami = as.numeric(af_ami),
    AF_amr = as.numeric(af_amr),
    AF_asj = as.numeric(af_asj),
    AF_eas = as.numeric(af_eas),
    AF_fin = as.numeric(af_fin),
    AF_mid = as.numeric(af_mid),
    AF_nfe = as.numeric(af_nfe),
    AF_remaining = as.numeric(af_remaining),
    AF_sas = as.numeric(af_sas)
  )
  
  return(snp_data)
}


# Define the path to the filtered VCF file
# vcf_file <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/gnomAD/combined_random_snps.vcf"
# snp_data <- process_vcf_with_af(vcf_file)


################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/gnomAD/")

ancestry = c("African/African-American", "Amish", "Admixed American", 
             "Ashkenazi Jewish", "East Asian", "Finnish", "Middle Eastern",
             "Non-Finnish European", "Remaining", "South Asian")
names(ancestry) = c("AF_afr", "AF_ami", "AF_amr", "AF_asj", "AF_eas", 
                    "AF_fin", "AF_mid", "AF_nfe", "AF_remaining", "AF_sas")

vcf_list = list.files(path = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/gnomAD/", pattern = "*bp_snps.vcf")
snp_data = data.frame()
for(vcf_file in vcf_list){
  snp_data_tmp <- process_vcf_with_af(vcf_file)
  if(nrow(snp_data_tmp)>0){
    snp_data = rbind(snp_data, snp_data_tmp)
  }
  print("test")
}

snp_data <- snp_data[snp_data$AF>=0.01, ]

snp_significance <- snp_data %>%
  pivot_longer(cols = starts_with("AF_"), names_to = "Genetic.Ancestry.Group", values_to = "AF_Value") %>%
  mutate(SNP = paste0(ID, "-", ALT)) %>%
  group_by(SNP) %>%
  mutate(
    SNP_Mean = mean(AF_Value),
    SD_Frequency = sd(AF_Value),
    Genetic.Ancestry.Group = ancestry[as.character(Genetic.Ancestry.Group)],
  ) %>%
  ungroup() %>%
  mutate(Significance = case_when(
    (AF_Value - SNP_Mean) > 2*SD_Frequency ~ "High", # 0.05
    (AF_Value - SNP_Mean) < -2*SD_Frequency ~ "Low",
    TRUE ~ "Baseline"
  ))

snp_ids <- snp_significance %>%
  filter(Significance %in% c("High", "Low")) %>%
  distinct(SNP)

filtered_data <- snp_significance %>%
  filter(SNP %in% snp_ids$SNP & ID %in% snp_use$SNP)

ggplot(filtered_data, aes(x = Genetic.Ancestry.Group, y = AF_Value, fill = Significance)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(SNP ~ ., scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = unit(c(0,0,0,1), "cm"),
        panel.spacing = unit(0.3, "lines"),
        legend.position = "right",
        axis.text = element_text(colour = 'black'),
        strip.background = element_rect(colour = "black", fill = NA)) +
  labs(x = "Genetic ancestry group", y = "Mean allele frequency", fill="Relative\nfrequency") +
  scale_fill_manual(values = c("High" = "red", "Low" = "blue", "Baseline" = "grey"))


gwas_bp_beta[gwas_bp_beta$SNP %in% filtered_data$ID, ]
unique(filtered_data$SNP)










################################################################################
# Define the number of iterations for the permutation test
num_iterations <- 100

# Define a function to check population-specific SNPs
check_population_specific <- function(snp_sample, threshold = 2) {
  # Calculate mean and SD across populations
  # mean_af <- rowMeans(snp_sample[, .(AF_afr, AF_ami, AF_amr, AF_asj, AF_eas, AF_fin, AF_mid, AF_nfe, AF_remaining, AF_sas)], na.rm = TRUE)
  # sd_af <- apply(snp_sample[, .(AF_afr, AF_ami, AF_amr, AF_asj, AF_eas, AF_fin, AF_mid, AF_nfe, AF_remaining, AF_sas)], 1, sd, na.rm = TRUE)

  mean_af <- rowMeans(snp_sample[, c("AF_afr", "AF_ami", "AF_amr", "AF_asj", "AF_eas", "AF_fin", "AF_mid", "AF_nfe", "AF_remaining", "AF_sas")], na.rm = TRUE)
  sd_af <- apply(snp_sample[, c("AF_afr", "AF_ami", "AF_amr", "AF_asj", "AF_eas", "AF_fin", "AF_mid", "AF_nfe", "AF_remaining", "AF_sas")], 1, sd, na.rm = TRUE)
  
  # Identify population-specific SNPs
  # population_specific <- abs(snp_sample$AF_afr - mean_af) > threshold * sd_af
  # population_specific <- snp_sample$AF_afr - mean_af > threshold * sd_af
  population_specific <- snp_sample$AF_afr - mean_af < -1 * threshold * sd_af
  sum(population_specific, na.rm = TRUE)
}

# Perform permutation test
set.seed(123)  # For reproducibility
results <- replicate(num_iterations, {
  # Sample 25 SNPs
  sampled_snps <- snp_data[sample(nrow(snp_data), 25), ]
  
  # Check for population-specific SNPs
  num_population_specific <- check_population_specific(sampled_snps)
  
  return(num_population_specific)
})

# Calculate p-value
observed_value <- 5  # Replace with your observed value
p_value <- mean(results >= observed_value)

# Print the p-value
print(paste("P-value:", p_value))





