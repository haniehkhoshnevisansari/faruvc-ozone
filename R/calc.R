# Dose, eACH, and CADR Calculations
# This code assumes you have already run the far-UVC lamp calculation code 
# and have 'results' dataframe with z, E_avg, and V_cone columns

# USER INPUT: Set susceptibility constant k
k_susceptibility <- 4.11  # cm²/mJ 

# Parameters 
exposure_time_hours <- 1
exposure_time_seconds <- exposure_time_hours * 3600  # Convert to seconds

if (!exists("results")) {
  stop("Error: 'results' dataframe not found. Please run the UV lamp calculation code first.")
}

results$dose_per_slice <- NA          # J/m² (dose per slice)
results$dose_per_slice_mJ_cm2 <- NA   # mJ/cm² (for eACH calculation)
results$eACH <- NA                    # equivalent Air Changes per Hour
results$CADR <- NA                    # Clean Air Delivery Rate (m³/h)

cat("Calculating dose, eACH, and CADR for each slice...\n")
cat("Slice\tZ (m)\tE_avg (W/m²)\tDose (J/m²)\tDose (mJ/cm²)\teACH\t\tCADR (m³/h)\tV_slice (m³)\n")
cat("----\t-----\t-----------\t----------\t------------\t----\t\t----------\t------------\n")

for (i in 1:nrow(results)) {
  z <- results$z[i]
  E_avg <- results$E_avg[i]
  V_slice <- results$V_cone[i]
  
  # Dose = Irradiance (W/m²) × Time (s) = J/m²
  dose_J_m2 <- E_avg * exposure_time_seconds
  
  # Convert dose to mJ/cm² for eACH calculation
  # Derivation: 1 J = 1000 mJ, 1 m² = 10,000 cm²
  # So 1 J/m² = 1000 mJ / 10,000 cm² = 0.1 mJ/cm²
  dose_mJ_cm2 <- dose_J_m2 * 0.1  
  
  # Calculate eACH (equivalent Air Changes per Hour)
  # For 1-hour exposure: eACH = k × dose (this gives air changes per hour directly)
  # Since k is in cm²/mJ and dose is in mJ/cm², the result is dimensionless per hour
  eACH <- k_susceptibility * dose_mJ_cm2
  
  # Calculate CADR (Clean Air Delivery Rate)
  # CADR = eACH × Volume of slice
  # eACH is per hour, so CADR will be in m³/h
  CADR <- eACH * V_slice
  
  results$dose_per_slice[i] <- dose_J_m2
  results$dose_per_slice_mJ_cm2[i] <- dose_mJ_cm2
  results$eACH[i] <- eACH
  results$CADR[i] <- CADR
  
  if (i %% 250 == 0 || i <= 10 || i > (nrow(results) - 10)) {
    cat(sprintf("%d\t%.3f\t%.6f\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.6f\n", 
                i, z, E_avg, dose_J_m2, dose_mJ_cm2, eACH, CADR, V_slice))
  }
}

cat("\n=== SUMMARY STATISTICS ===\n")
cat(sprintf("Total number of slices: %d\n", nrow(results)))
cat(sprintf("Z range: %.3f to %.3f m\n", min(results$z), max(results$z)))
cat(sprintf("E_avg range: %.6f to %.6f W/m²\n", min(results$E_avg), max(results$E_avg)))
cat(sprintf("Dose range: %.4f to %.4f J/m² (%.4f to %.4f mJ/cm²)\n", 
            min(results$dose_per_slice), max(results$dose_per_slice),
            min(results$dose_per_slice_mJ_cm2), max(results$dose_per_slice_mJ_cm2)))
cat(sprintf("eACH range: %.6f to %.6f h⁻¹\n", min(results$eACH), max(results$eACH)))
cat(sprintf("CADR range: %.6f to %.6f m³/h\n", min(results$CADR), max(results$CADR)))

par(mfrow = c(2, 2))

# Dose vs Z
plot(results$z, results$dose_per_slice_mJ_cm2, type = "l", lwd = 2,
     xlab = "Z height (m)", ylab = "Dose (mJ/cm²)",
     main = sprintf("UV Dose vs Distance (1h exposure, k = %.2f cm²/mJ)", k_susceptibility),
     col = "blue")
grid(TRUE)

# eACH vs Z
plot(results$z, results$eACH, type = "l", lwd = 2,
     xlab = "Z height (m)", ylab = "eACH (h⁻¹)",
     main = "Equivalent Air Changes per Hour vs Distance",
     col = "red")
grid(TRUE)

# CADR vs Z
plot(results$z, results$CADR, type = "l", lwd = 2,
     xlab = "Z height (m)", ylab = "CADR (m³/h)",
     main = "Clean Air Delivery Rate vs Distance",
     col = "green")
grid(TRUE)

par(mfrow = c(1, 1))

write.csv(results, "uv_lamp_dose_each_cadr_results.csv", row.names = FALSE)
cat(sprintf("\nResults saved to 'uv_lamp_dose_each_cadr_results.csv'\n"))

# Function to recalculate with different k value - For 1 hour exposure
recalculate_with_new_k <- function(new_k) {
  cat(sprintf("\n=== RECALCULATING WITH k = %.2f cm²/mJ ===\n", new_k))
  
  # Recalculate eACH and CADR with new k - For 1 hour exposure
  results$eACH_new <- new_k * results$dose_per_slice_mJ_cm2  # Direct calculation for 1 hour
  results$CADR_new <- results$eACH_new * results$V_cone
  
  total_CADR_new <- sum(results$CADR_new)
  volume_weighted_eACH_new <- sum(results$eACH_new * results$V_cone) / total_volume
  
  cat(sprintf("Total CADR with new k: %.2f m³/h\n", total_CADR_new))
  cat(sprintf("Volume-weighted average eACH with new k: %.6f h⁻¹\n", volume_weighted_eACH_new))
  
  return(list(
    total_CADR = total_CADR_new,
    volume_weighted_eACH = volume_weighted_eACH_new,
    results_with_new_k = results[, c("z", "E_avg", "V_cone", "dose_per_slice_mJ_cm2", "eACH_new", "CADR_new")]
  ))
}

cat("• To recalculate with a different k value, use: recalculate_with_new_k(new_k_value)\n")
