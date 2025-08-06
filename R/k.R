# Dose, eACH, and CADR Calculations for Multiple Susceptibility Constants
# This code assumes you have already run the UV lamp calculation code 
# and have 'results' dataframe with z, E_avg, and V_cone columns

# USER INPUT: Define all susceptibility constants
k_values <- c(4.11, 3, 5.9, 1.8, 1.14, 0.22, 0.41, 5.75, 17.14, 20.08, 14.26)
pathogen_names <- c("Human coronavirus", "Phage P22", "Human coronavirus", "Influenza virus", "Human coronavirus", 
                   "Phage MS2", "Phage MS2", "Phage MS2", "Phage P22", "Phage Phi6", 
                   "Human coronavirus")

unique_pathogen_ids <- paste0(pathogen_names, "_k", k_values)
unique_pathogen_ids

cat("Susceptibility constants (k values):\n")
for(i in 1:length(k_values)) {
  cat(sprintf("  %s: %.2f cm²/mJ\n", pathogen_names[i], k_values[i]))
}
cat("\n")

exposure_time_hours <- 1
exposure_time_seconds <- exposure_time_hours * 3600  # Convert to seconds

if (!exists("results")) {
  stop("Error: 'results' dataframe not found. Please run the UV lamp calculation code first.")
}

target_z <- 2.45
z_index <- which.min(abs(results$z - target_z))
actual_z <- results$z[z_index]

cat(sprintf("Target z: %.2f m, Actual z used: %.3f m (index: %d)\n", target_z, actual_z, z_index))


results$dose_per_slice <- NA          # J/m² (dose per slice)
results$dose_per_slice_mJ_cm2 <- NA   # mJ/cm² (for eACH calculation)


for(i in 1:length(k_values)) {
  pathogen_id <- unique_pathogen_ids[i]
  results[[paste0("eACH_", pathogen_id)]] <- NA
  results[[paste0("CADR_", pathogen_id)]] <- NA
}

cat("Calculating dose, eACH, and CADR for each slice and pathogen...\n")

for (i in 1:nrow(results)) {
  z <- results$z[i]
  E_avg <- results$E_avg[i]
  V_slice <- results$V_cone[i]
  
  # Calculate dose for 1 hour
  # Dose = Irradiance (W/m²) × Time (s) = J/m²
  dose_J_m2 <- E_avg * exposure_time_seconds
  
  # Convert dose to mJ/cm² for eACH calculation
  # Derivation: 1 J = 1000 mJ, 1 m² = 10,000 cm²
  # So 1 J/m² = 1000 mJ / 10,000 cm² = 0.1 mJ/cm²
  dose_mJ_cm2 <- dose_J_m2 * 0.1  
  
  # Store dose results (same for all pathogens)
  results$dose_per_slice[i] <- dose_J_m2
  results$dose_per_slice_mJ_cm2[i] <- dose_mJ_cm2
  
  # Calculate eACH and CADR for each pathogen
  for(j in 1:length(k_values)) {
    pathogen_id <- unique_pathogen_ids[j]
    k_susceptibility <- k_values[j]
    
    # Calculate eACH (equivalent Air Changes per Hour)
    # For 1-hour exposure: eACH = k × dose (this gives air changes per hour directly)
    eACH <- k_susceptibility * dose_mJ_cm2
    
    # Calculate CADR (Clean Air Delivery Rate)
    # CADR = eACH × Volume of slice
    CADR <- eACH * V_slice
    
    results[[paste0("eACH_", pathogen_id)]][i] <- eACH
    results[[paste0("CADR_", pathogen_id)]][i] <- CADR
  }
}

cat("\n=== SUMMARY STATISTICS BY PATHOGEN ===\n")
cat(sprintf("Total number of slices: %d\n", nrow(results)))
cat(sprintf("Z range: %.3f to %.3f m\n", min(results$z), max(results$z)))
cat(sprintf("E_avg range: %.6f to %.6f W/m²\n", min(results$E_avg), max(results$E_avg)))
cat(sprintf("Dose range: %.4f to %.4f J/m² (%.4f to %.4f mJ/cm²)\n", 
            min(results$dose_per_slice), max(results$dose_per_slice),
            min(results$dose_per_slice_mJ_cm2), max(results$dose_per_slice_mJ_cm2)))

pathogen_summary <- data.frame(
  Pathogen = pathogen_names,
  Pathogen_ID = unique_pathogen_ids,
  K_Value = k_values,
  CADR_at_z245 = numeric(length(k_values)),  # CADR at z=2.45m
  Volume_Weighted_eACH = numeric(length(k_values)),
  Min_eACH = numeric(length(k_values)),
  Max_eACH = numeric(length(k_values)),
  Min_CADR = numeric(length(k_values)),
  Max_CADR = numeric(length(k_values)),
  eACH_at_z245 = numeric(length(k_values))   # eACH at z=2.45m
)

total_volume <- sum(results$V_cone)

for(i in 1:length(k_values)) {
  pathogen_id <- unique_pathogen_ids[i]
  
  eACH_col <- paste0("eACH_", pathogen_id)
  CADR_col <- paste0("CADR_", pathogen_id)
  
  pathogen_summary$CADR_at_z245[i] <- results[[CADR_col]][z_index]
  pathogen_summary$eACH_at_z245[i] <- results[[eACH_col]][z_index]
  
  pathogen_summary$Min_eACH[i] <- min(results[[eACH_col]])
  pathogen_summary$Max_eACH[i] <- max(results[[eACH_col]])
  pathogen_summary$Min_CADR[i] <- min(results[[CADR_col]])
  pathogen_summary$Max_CADR[i] <- max(results[[CADR_col]])
  
  cat(sprintf("\n%s (k = %.2f cm²/mJ):\n", pathogen_names[i], k_values[i]))
  cat(sprintf("  CADR at z=%.2fm: %.2f m³/h\n", actual_z, pathogen_summary$CADR_at_z245[i]))
  cat(sprintf("  eACH at z=%.2fm: %.6f h⁻¹\n", actual_z, pathogen_summary$eACH_at_z245[i]))
  cat(sprintf("  eACH range: %.6f to %.6f h⁻¹\n", pathogen_summary$Min_eACH[i], pathogen_summary$Max_eACH[i]))
  cat(sprintf("  CADR range: %.6f to %.6f m³/h\n", pathogen_summary$Min_CADR[i], pathogen_summary$Max_CADR[i]))
}

write.csv(results, "uv_lamp_dose_each_cadr_all_pathogens_results.csv", row.names = FALSE)
cat(sprintf("\nDetailed results saved to 'uv_lamp_dose_each_cadr_all_pathogens_results.csv'\n"))

write.csv(pathogen_summary, "pathogen_cadr_summary.csv", row.names = FALSE)
cat(sprintf("Summary results saved to 'pathogen_cadr_summary.csv'\n"))

par(mfrow = c(2, 2))

# CADR comparison - sorted by k value
sorted_indices <- order(pathogen_summary$K_Value)
barplot(pathogen_summary$CADR_at_z245[sorted_indices], 
        names.arg = pathogen_summary$Pathogen_ID[sorted_indices],
        las = 2, cex.names = 0.6,
        main = sprintf("CADR at z=%.2fm by Pathogen (sorted by k)", actual_z),
        ylab = "CADR (m³/h)",
        col = rainbow(length(k_values)))

# eACH comparison - sorted by k value
barplot(pathogen_summary$eACH_at_z245[sorted_indices], 
        names.arg = pathogen_summary$Pathogen_ID[sorted_indices],
        las = 2, cex.names = 0.6,
        main = sprintf("eACH at z=%.2fm by Pathogen (sorted by k)", actual_z),
        ylab = "eACH (h⁻¹)",
        col = rainbow(length(k_values)))

# CADR vs Distance for selected pathogens (top 4 by k value)
top_4_indices <- order(k_values, decreasing = TRUE)[1:4]
colors <- c("red", "blue", "green", "purple")

plot(results$z, results[[paste0("CADR_", unique_pathogen_ids[top_4_indices[1]])]], 
     type = "l", lwd = 2, col = colors[1],
     xlab = "Z height (m)", ylab = "CADR (m³/h)",
     main = "CADR vs Distance (Top 4 Pathogens)")

abline(v = actual_z, col = "black", lty = 2, lwd = 2)
text(actual_z, max(results[[paste0("CADR_", unique_pathogen_ids[top_4_indices[1]])]])*0.9, 
     sprintf("z=%.2fm", actual_z), pos = 4)

for(i in 2:4) {
  lines(results$z, results[[paste0("CADR_", unique_pathogen_ids[top_4_indices[i]])]], 
        lwd = 2, col = colors[i])
}

legend("topright", legend = unique_pathogen_ids[top_4_indices], 
       col = colors, lwd = 2, cex = 0.7)
grid(TRUE)

# K values comparison - sorted
barplot(pathogen_summary$K_Value[sorted_indices], 
        names.arg = pathogen_summary$Pathogen_ID[sorted_indices],
        las = 2, cex.names = 0.6,
        main = "Susceptibility Constants (k) by Pathogen (sorted)",
        ylab = "k (cm²/mJ)",
        col = rainbow(length(k_values)))

par(mfrow = c(1, 1))

# Function to get CADR for specific pathogen and multiplier
get_pathogen_cadr <- function(pathogen_name, multiplier = 1) {
  if (!pathogen_name %in% pathogen_names) {
    stop(paste("Pathogen not found. Available pathogens:", paste(unique(pathogen_names), collapse = ", ")))
  }
  
  pathogen_index <- which(pathogen_names == pathogen_name)[1]  # Take first match
  base_cadr <- pathogen_summary$CADR_at_z245[pathogen_index]  # CORRECTED: Use z=2.45m value
  adjusted_cadr <- base_cadr * multiplier
  
  cat(sprintf("Pathogen: %s\n", pathogen_name))
  cat(sprintf("Base CADR at z=%.2fm: %.2f m³/h\n", actual_z, base_cadr))
  cat(sprintf("Multiplier: %d\n", multiplier))
  cat(sprintf("Adjusted CADR: %.2f m³/h\n", adjusted_cadr))
  
  return(adjusted_cadr)
}

pathogen_cadr_summary <<- pathogen_summary   
