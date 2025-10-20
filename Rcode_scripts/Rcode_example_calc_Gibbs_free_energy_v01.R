
# Your provided example data frame
df_Ml_sub <- structure(list(speciesabbr = "Mnelei", 
                            FprimNm = "Mnelei_its2_F04", 
                            Fprimseq = "ACGGTCCCTTGAAGTAGAGC", 
                            RprimNm = "Mnelei_its2_R06", 
                            RprimSeq = "TCTGAGAAGGCTTCGGACAT", 
                            ProbeSeq = "GTGCCTCTCGGTGTGGTAGCAATATCT", 
                            Targetfragm = "ACGGTCCCTTGAAGTAGAGCGATCCCGAGTTGCGCGGAAGCCCCTTCGTGGCCCGTGCCTCTCGGTGTGGTAGCAATATCTCACCGAGCGGCCGAGCTTTGCAGGATGTCCGAAGCCTTCTCAGA"
), row.names = 19L, class = "data.frame")

# Nearest-neighbor thermodynamic parameters (SantaLucia 1998)
nn_params <- data.frame(
  dinuc = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"),
  dH     = c(-7.9, -8.4, -7.8, -7.2, -8.5, -8.0, -10.6, -7.8, -8.2, -9.8, -8.0, -8.4, -7.2, -8.2, -8.5, -7.9),
  dS     = c(-22.2, -22.4, -21.0, -20.4, -22.7, -19.9, -27.2, -21.0, -22.2, -24.4, -19.9, -22.4, -21.3, -22.2, -22.7, -22.2)
)

#1. Function to Calculate Tm (Melting Temperature)
# Function to calculate Tm using nearest-neighbor method
calculate_tm <- function(dna_seq, Na_conc = 40, primer_conc = 50e-5) {
  seq <- toupper(dna_seq)
  n <- nchar(seq)
  if (n == 0) return(NA)
  
  # Initialize
  dH <- 0
  dS <- 0
  
  # Calculate ΔH and ΔS
  for (i in 1:(n-1)) {
    dinuc <- paste(substr(seq, i, i), substr(seq, i+1, i+1), sep = "")
    params <- nn_params[nn_params$dinuc == dinuc, ]
    if (nrow(params) > 0) {
      dH <- dH + params$dH[1]  # Take the first match
      dS <- dS + params$dS[1]
    }
  }
  
  # Add initiation and symmetry corrections
  dH <- dH + 0.2  # initiation
  dS <- dS - 5.7  # initiation
  
  # Salt correction
  R <- 1.987  # Gas constant in cal/mol·K
  salt_corr <- 16.6 * log10(Na_conc / 1000)
  
  # Calculate Tm
  tm <- (1000 * dH) / (dS + R * log(primer_conc)) - 273.15 + salt_corr
  
  return(tm)
}
# 2. Function to Calculate ΔG (Gibbs Free Energy) at a Given Temperature
# Function to calculate ΔG (Gibbs free energy) in kcal/mol at a given temperature
calculate_deltaG <- function(dna_seq, Na_conc = 40, primer_conc = 50e-5, temp_C = 60) {
  seq <- toupper(dna_seq)
  n <- nchar(seq)
  if (n == 0) return(NA)
  
  # Initialize
  dH <- 0
  dS <- 0
  
  # Calculate ΔH and ΔS
  for (i in 1:(n-1)) {
    dinuc <- paste(substr(seq, i, i), substr(seq, i+1, i+1), sep = "")
    params <- nn_params[nn_params$dinuc == dinuc, ]
    if (nrow(params) > 0) {
      dH <- dH + params$dH[1]  # Take the first match
      dS <- dS + params$dS[1]
    }
  }
  
  # Add initiation and symmetry corrections
  dH <- dH + 0.2  # initiation
  dS <- dS - 5.7  # initiation
  
  # Convert temperature to Kelvin
  temp_K <- temp_C + 273.15
  
  # Calculate ΔG in kcal/mol
  deltaG <- dH - (temp_K / 1000) * dS
  
  return(deltaG)
}

# Calculate Tm
df_Ml_sub$FprimTm <- calculate_tm(df_Ml_sub$Fprimseq)
df_Ml_sub$RprimTm <- calculate_tm(df_Ml_sub$RprimSeq)
df_Ml_sub$ProbeTm <- calculate_tm(df_Ml_sub$ProbeSeq)

# Calculate ΔG at 60°C
df_Ml_sub$FprimDeltaG <- calculate_deltaG(df_Ml_sub$Fprimseq, temp_C = 60)
df_Ml_sub$RprimDeltaG <- calculate_deltaG(df_Ml_sub$RprimSeq, temp_C = 60)
df_Ml_sub$ProbeDeltaG <- calculate_deltaG(df_Ml_sub$ProbeSeq, temp_C = 60)

# Print results
print(df_Ml_sub[, c("FprimNm", "FprimTm", "FprimDeltaG",
                    "RprimNm", "RprimTm", "RprimDeltaG",
                    "ProbeSeq", "ProbeTm", "ProbeDeltaG")])

