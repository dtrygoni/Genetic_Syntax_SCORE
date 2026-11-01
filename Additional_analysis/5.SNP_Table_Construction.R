library(readxl)
data<-read_excel("./GESS_SNPs_final.xlsx")
library(SNPassoc)
data$Outcome<-ifelse(data$TOTAL.SS>0,1,0)
SNP_COLUMNS <- grep(
  "^rs",
  names(data),
  value = TRUE
)

GENOTYPE_SEPARATOR <- ""

if (length(SNP_COLUMNS) != 228) {
  warning(
    paste(
      length(SNP_COLUMNS),
      "SNP columns detected; 228 were expected."
    )
  )
}
dat<-data

snp_column_numbers <- match(
  SNP_COLUMNS,
  names(dat)
)
snp_data <- SNPassoc::setupSNP(
  data = dat,
  colSNPs = snp_column_numbers,
  sep = GENOTYPE_SEPARATOR
)




res<-tableHWE(snp_data,Outcome)
summary_ss<-summary(snp_data,info=TRUE)


snps_all<-row.names(res)
HWE_p<-round(res[,2],3)
missingness<-summary_ss$missing
overall_MAF<-100-summary_ss$major.allele.freq
alleles<-summary_ss$alleles
minor_effect_alleles<-c()
for(allele in alleles){
  r<-strsplit(allele,split="/")[[1]]
  if(length(r)==2){
    minor_effect_alleles<-c(minor_effect_alleles,r[2])
  }else{
    minor_effect_alleles<-c(minor_effect_alleles,NA)
  }
}




### Call-rates
calculate_call_rates <- function(
    data,
    snp_columns,
    group = NULL
) {
  
  if (!all(snp_columns %in% names(data))) {
    stop("Some SNP columns were not found in the dataframe.")
  }
  
  if (
    !is.null(group) &&
    length(group) != nrow(data)
  ) {
    stop("The group vector must have one value per participant.")
  }
  
  calculate_one <- function(
    snp_name,
    index,
    group_name
  ) {
    
    genotype <- data[[snp_name]][index]
    
    n_total <- length(genotype)
    n_missing <- sum(is.na(genotype))
    n_called <- n_total - n_missing
    
    data.frame(
      SNP = snp_name,
      Group = group_name,
      N = n_total,
      Called_n = n_called,
      Missing_n = n_missing,
      Call_rate = n_called / n_total,
      Call_rate_percent = 100 * n_called / n_total,
      stringsAsFactors = FALSE
    )
  }
  
  results <- list()
  
  # Overall call rates
  for (snp_name in snp_columns) {
    
    results[[length(results) + 1]] <-
      calculate_one(
        snp_name = snp_name,
        index = rep(TRUE, nrow(data)),
        group_name = "Overall"
      )
  }
  
  # Group-specific call rates
  if (!is.null(group)) {
    
    group <- as.character(group)
    group_levels <- unique(group[!is.na(group)])
    
    for (group_value in group_levels) {
      
      group_index <- !is.na(group) &
        group == group_value
      
      for (snp_name in snp_columns) {
        
        results[[length(results) + 1]] <-
          calculate_one(
            snp_name = snp_name,
            index = group_index,
            group_name = group_value
          )
      }
    }
  }
  
  results <- do.call(
    rbind,
    results
  )
  
  rownames(results) <- NULL
  
  results
}
call_rates <- calculate_call_rates(
  data = snp_data,
  snp_columns = SNP_COLUMNS
)

call_rate_vector<-round(call_rates$Call_rate,3)



res_df<-data.frame(SNP=snps_all,Call.rate=call_rate_vector,Minor.Allele=minor_effect_alleles,Overall.MAF=overall_MAF,HWE.p.value=HWE_p,Missing=missingness)







################ MAF per group ###################
calculate_allele_frequency <- function(
    genotype,
    allele,
    participant_index
) {
  
  # Select the relevant participants
  genotype <- genotype[participant_index]
  
  # Convert SNPassoc SNP object to character
  genotype <- as.character(genotype)
  
  # Standardize missing values
  genotype[genotype %in% c(
    "",
    "NA",
    "N/A",
    "./.",
    ".|.",
    "--"
  )] <- NA_character_
  
  # Standardize separators
  genotype <- gsub(
    "\\|",
    "/",
    genotype
  )
  
  genotype <- gsub(
    "\\s+",
    "",
    genotype
  )
  
  parse_one <- function(x) {
    
    if (is.na(x)) {
      return(c(
        NA_character_,
        NA_character_
      ))
    }
    
    # Genotypes such as A/G
    if (grepl("/", x, fixed = TRUE)) {
      
      alleles <- strsplit(
        x,
        split = "/",
        fixed = TRUE
      )[[1]]
      
      # Genotypes such as AG
    } else if (nchar(x) == 2) {
      
      alleles <- strsplit(
        x,
        split = ""
      )[[1]]
      
    } else {
      
      return(c(
        NA_character_,
        NA_character_
      ))
    }
    
    if (length(alleles) != 2) {
      return(c(
        NA_character_,
        NA_character_
      ))
    }
    
    alleles
  }
  
  parsed <- do.call(
    rbind,
    lapply(genotype, parse_one)
  )
  
  valid <- complete.cases(parsed)
  
  if (sum(valid) == 0) {
    return(NA_real_)
  }
  
  observed_alleles <- c(
    parsed[valid, 1],
    parsed[valid, 2]
  )
  
  # Return percentage to match Overall.MAF
  100 * mean(observed_alleles == allele)
}
zero_index <- dat$Outcome == 0
nonzero_index <- dat$Outcome == 1
res_df$Minor_allele_frequency_zero <- mapply(
  FUN = function(snp_name, minor_allele) {
    
    calculate_allele_frequency(
      genotype = snp_data[[snp_name]],
      allele = minor_allele,
      participant_index = zero_index
    )
  },
  snp_name = res_df$SNP,
  minor_allele = res_df$Minor.Allele
)

res_df$Minor_allele_frequency_nonzero <- mapply(
  FUN = function(snp_name, minor_allele) {
    
    calculate_allele_frequency(
      genotype = snp_data[[snp_name]],
      allele = minor_allele,
      participant_index = nonzero_index
    )
  },
  snp_name = res_df$SNP,
  minor_allele = res_df$Minor.Allele
)

res_df$Minor_allele_frequency_zero<-round(res_df$Minor_allele_frequency_zero,1)
res_df$Minor_allele_frequency_nonzero<-round(res_df$Minor_allele_frequency_nonzero,1)
writexl::write_xlsx(res_df,path='.../Supplementary_Table_3.xlsx')
