# Photon Flux and Ozone Production Rate Calculations
# Based on Link et al. 2023 methodology for GUV222 lamps
# This code extends the UV lamp calculation results with ozone production modeling
# Ozone generation reporting in mg/h

if (!exists("results")) {
  stop("Error: 'results' dataframe not found. Please run the UV lamp calculation code first.")
}

# Physical constants and parameters
h <- 6.62607015e-34  # Planck constant (J·s)
c <- 299792458       # Speed of light (m/s)
lambda_222 <- 222e-9 # Wavelength 222 nm in meters
N_A <- 6.02214076e23 # Avogadro's number (mol⁻¹)

# Ozone molecular parameters
M_O3 <- 47.998       # Molecular weight of O3 (g/mol)

# O2 photolysis parameters 
sigma_O2 <- 4.30e-24    # O2 absorption cross section (cm²)
phi_O2 <- 1.0           # Quantum yield for O2 photolysis

# Rate constants (for 20°C)
k1 <- 7.96e-15         # O + O3 -> 2O2 rate constant (cm³ molecule⁻¹ s⁻¹)
k2 <- 6.10e-34         # O + O2 + M -> O3 + M rate constant (cm⁶ molecule⁻² s⁻¹)

# Air composition and density (at 20°C, 1 atm)
M_air <- 2.46e19       # Air number density (molecules cm⁻³)
O2_fraction <- 0.21    # Fraction of O2 in air
O2_density <- M_air * O2_fraction  # O2 number density (molecules cm⁻³)

results$photon_flux <- NA          # Photon flux (photons cm⁻² s⁻¹)
results$j_O2 <- NA                 # O2 photolysis rate constant (s⁻¹)
results$O3_production_rate <- NA   # O3 production rate (molecules cm⁻³ s⁻¹)
results$O3_production_ppbv_h <- NA # O3 production rate (ppbv h⁻¹)
results$O3_production_mg_h <- NA   # O3 production rate (mg h⁻¹)

cat("\nCalculating photon flux and ozone production for each slice...\n")
cat("Slice\tZ (m)\tE_avg (W/m²)\tPhoton Flux\t\tj_O2 (s⁻¹)\tO3 Prod Rate\tO3 Prod (ppbv/h)\tO3 Prod (mg/h)\n")
cat("\t\t\t\t(photons/cm²/s)\t\t\t(molec/cm³/s)\n")
cat("----\t-----\t-----------\t--------------\t----------\t-----------\t---------------\t-------------\n")

# Photon energy at 222 nm
E_photon <- h * c / lambda_222  # Energy per photon (J)

for (i in 1:nrow(results)) {
  z <- results$z[i]
  E_avg_W_m2 <- results$E_avg[i]  # Irradiance in W/m²
  V_cone_cm3 <- results$V_cone[i] * 1e6  # Convert m³ to cm³
  
  # Convert irradiance to photon flux
  # E_avg (W/m²) = E_avg (J/s/m²)
  # Photon flux = E_avg / E_photon (photons/m²/s)
  # Convert to photons/cm²/s: divide by 10,000
  photon_flux_cm2_s <- (E_avg_W_m2 / E_photon) / 10000
  
  # Calculate O2 photolysis rate constant (j_O2)
  # j_O2 = σ_O2 × Φ_O2 × F 
  # where F is photon flux in photons cm⁻² s⁻¹
  j_O2 <- sigma_O2 * phi_O2 * photon_flux_cm2_s
  
  # Calculate O3 production rate (molecules cm⁻³ s⁻¹)
  # From the steady-state assumption and Equation 2:
  # d[O3]/dt = 2 × j_O2 × [O2] - k1 × [O] × [O3] - k_loss × [O3]
  # At steady state for O atoms: [O] = (j_O2 × [O2]) / (k1 × [O3] + k2 × [O2] × [M])
  # For simplicity and following the paper's approach, the net O3 production rate is:
  # ≈ 2 × j_O2 × [O2] (the factor of 2 comes from each O2 molecule producing 2 O atoms)
  O3_production_rate <- 2 * j_O2 * O2_density
  
  # Convert O3 production rate to ppbv/h
  # molecules cm⁻³ s⁻¹ to ppbv h⁻¹:
  # ppbv = (molecules cm⁻³) / (total molecules cm⁻³) × 10⁹
  # ppbv h⁻¹ = (molecules cm⁻³ s⁻¹) / M_air × 10⁹ × 3600 s/h
  O3_production_ppbv_h <- (O3_production_rate / M_air) * 1e9 * 3600
  
  # Calculate O3 production rate in mg/h 
  # molecules cm⁻³ s⁻¹ → molecules s⁻¹ → mol s⁻¹ → mg s⁻¹ → mg h⁻¹
  O3_molecules_per_s <- O3_production_rate * V_cone_cm3  # molecules/s 
  O3_mol_per_s <- O3_molecules_per_s / N_A               # mol/s
  O3_g_per_s <- O3_mol_per_s * M_O3                      # g/s
  O3_mg_per_h <- O3_g_per_s * 1000 * 3600               # mg/h
  
  results$photon_flux[i] <- photon_flux_cm2_s
  results$j_O2[i] <- j_O2
  results$O3_production_rate[i] <- O3_production_rate
  results$O3_production_ppbv_h[i] <- O3_production_ppbv_h
  results$O3_production_mg_h[i] <- O3_mg_per_h
  
  if (i %% 250 == 0 || i <= 10 || i > (nrow(results) - 10)) {
    cat(sprintf("%d\t%.3f\t%.6f\t\t%.2e\t\t%.2e\t%.2e\t%.4f\t\t%.6f\n", 
                i, z, E_avg_W_m2, photon_flux_cm2_s, j_O2, O3_production_rate, O3_production_ppbv_h, O3_mg_per_h))
  }
}

cat(sprintf("Total number of slices: %d\n", nrow(results)))
cat(sprintf("Photon flux range: %.2e to %.2e photons cm⁻² s⁻¹\n", 
            min(results$photon_flux), max(results$photon_flux)))
cat(sprintf("j_O2 range: %.2e to %.2e s⁻¹\n", 
            min(results$j_O2), max(results$j_O2)))
cat(sprintf("O3 production rate range: %.2e to %.2e molecules cm⁻³ s⁻¹\n", 
            min(results$O3_production_rate), max(results$O3_production_rate)))
cat(sprintf("O3 production rate range: %.4f to %.4f ppbv h⁻¹\n", 
            min(results$O3_production_ppbv_h), max(results$O3_production_ppbv_h)))
cat(sprintf("O3 production rate range: %.6f to %.6f mg h⁻¹\n", 
            min(results$O3_production_mg_h), max(results$O3_production_mg_h)))


par(mfrow = c(2, 4))

# Photon flux vs Z
plot(results$z, results$photon_flux, type = "l", lwd = 2,
     xlab = "Z height (m)", ylab = "Photon flux (photons cm⁻² s⁻¹)",
     main = "Photon Flux vs Distance",
     col = "blue")
grid(TRUE)

# j_O2 vs Z
plot(results$z, results$j_O2, type = "l", lwd = 2,
     xlab = "Z height (m)", ylab = "j_O2 (s⁻¹)",
     main = "O2 Photolysis Rate vs Distance",
     col = "red")
grid(TRUE)

# O3 production rate vs Z
plot(results$z, results$O3_production_rate, type = "l", lwd = 2,
     xlab = "Z height (m)", ylab = "O3 production (molecules cm⁻³ s⁻¹)",
     main = "O3 Production Rate vs Distance",
     col = "green")
grid(TRUE)

# O3 production in ppbv/h vs Z
plot(results$z, results$O3_production_ppbv_h, type = "l", lwd = 2,
     xlab = "Z height (m)", ylab = "O3 production (ppbv h⁻¹)",
     main = "O3 Production Rate vs Distance",
     col = "orange")
grid(TRUE)

# O3 production in mg/h vs Z
plot(results$z, results$O3_production_mg_h, type = "l", lwd = 2,
     xlab = "Z height (m)", ylab = "O3 production (mg h⁻¹)",
     main = "O3 Production Rate vs Distance",
     col = "darkgreen")
grid(TRUE)

par(mfrow = c(1, 1))

write.csv(results, "uv_lamp_ozone_production_results.csv", row.names = FALSE)
cat(sprintf("\nOzone production results saved to 'uv_lamp_ozone_production_results.csv'\n"))

# Enhanced function to recalculate with different O2 absorption cross section
recalculate_ozone_with_new_sigma <- function(new_sigma_O2) {
  cat(sprintf("\n=== RECALCULATING OZONE WITH σ_O2 = %.2e cm² ===\n", new_sigma_O2))
  
  results$j_O2_new <- new_sigma_O2 * phi_O2 * results$photon_flux
  results$O3_production_rate_new <- 2 * results$j_O2_new * O2_density
  results$O3_production_ppbv_h_new <- (results$O3_production_rate_new / M_air) * 1e9 * 3600
  
  # Calculate new mg/h values
  V_cone_cm3 <- results$V_cone * 1e6
  O3_molecules_per_s_new <- results$O3_production_rate_new * V_cone_cm3
  O3_mol_per_s_new <- O3_molecules_per_s_new / N_A
  O3_g_per_s_new <- O3_mol_per_s_new * M_O3
  results$O3_production_mg_h_new <- O3_g_per_s_new * 1000 * 3600
  
  total_O3_new_ppbv <- sum(results$O3_production_ppbv_h_new * results$V_cone) / sum(results$V_cone)
  total_O3_new_mg_h <- sum(results$O3_production_mg_h_new)
  
  cat(sprintf("New volume-weighted average O3 production rate: %.4f ppbv h⁻¹\n", total_O3_new_ppbv))
  cat(sprintf("New total O3 production rate: %.3f mg h⁻¹\n", total_O3_new_mg_h))
  
  return(list(
    total_O3_production_ppbv = total_O3_new_ppbv,
    total_O3_production_mg_h = total_O3_new_mg_h,
    results_with_new_sigma = results[, c("z", "photon_flux", "j_O2_new", "O3_production_ppbv_h_new", "O3_production_mg_h_new")]
  ))
}

cat("• To recalculate with different σ_O2, use: recalculate_ozone_with_new_sigma(new_value)\n")
