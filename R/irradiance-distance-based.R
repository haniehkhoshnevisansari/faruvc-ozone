# UV Lamp E_avg Calculation
# Based on equations from the supplemental material of Link et al. (2023)

library(pracma) 

E0 <- 0.286  # W/m^(2-a) 
a <- 1.52   

# Angular distribution parameters
A <- 0.00479
B <- 0.9809
b <- 0.7180
c <- 0.1196

V_chamber <- 31.5  # m^3 
theta_max <- 55 * pi/180  
h_max <- 2.45 # m


z_min <- 0.001  

irradiance_func <- function(z, rho) {
  if (z < z_min) {
    return(105)  # W/m^2 - maximum irradiance at lamp surface
  }
  
  r <- sqrt(z^2 + rho^2)
  
  E_mag <- E0 / (r^a)
  
  if (z <= 0) {
    return(105)  # Safety check
  }
  
  angle_ratio <- rho / z
  angular_factor <- A + (B - A) / (1 + exp((angle_ratio - b) / c))
  
  result <- E_mag * angular_factor
  
  if (result < 0 || is.na(result) || !is.finite(result)) {
    cat(sprintf("Warning: Invalid result at z=%.6f, rho=%.6f: E_mag=%.6f, angular_factor=%.6f, result=%.6f\n", 
                z, rho, E_mag, angular_factor, result))
    return(105) 
  }
  
  return(result)
}

# Function to calculate E_avg for a given z value (height of cone)
calculate_E_avg <- function(z_height) {
  rho_max <- z_height * tan(theta_max)
  
  V_cone <- (1/3) * pi * rho_max^2 * z_height
  
  integrand <- function(params) {
    z <- params[1]
    rho <- params[2]
    phi <- params[3]
    
    # The Jacobian for cylindrical coordinates is rho
    return(irradiance_func(z, rho) * rho)
  }
  
  # Perform numerical integration using adaptIntegral for 3D
  # We need to integrate over z from z_min to z_height,
  # rho from 0 to rho(z) = z*tan(theta_max), and phi from 0 to 2*pi
  # Since the function is symmetric in phi, we can multiply by 2*pi
  # and integrate only over z and rho
  
  integrand_2d <- function(z_val, rho_val) {
    max_rho_at_z <- z_val * tan(theta_max)
    if (rho_val > max_rho_at_z) {
      return(0)
    }
    return(irradiance_func(z_val, rho_val) * rho_val)
  }
  
  total_flux <- 0
  
  z_points <- seq(z_min, z_height, length.out = 50)
  
  for (i in 1:(length(z_points)-1)) {
    z_low <- z_points[i]
    z_high <- z_points[i+1]
    z_mid <- (z_low + z_high) / 2
    
    rho_max_at_z <- z_mid * tan(theta_max)
    
    if (rho_max_at_z > 0) {
      # Integrate over rho at this z level
      rho_points <- seq(0, rho_max_at_z, length.out = 30)
      
      for (j in 1:(length(rho_points)-1)) {
        rho_low <- rho_points[j]
        rho_high <- rho_points[j+1]
        rho_mid <- (rho_low + rho_high) / 2
        
        dz <- z_high - z_low
        drho <- rho_high - rho_low
        dphi <- 2 * pi
        
        dV <- rho_mid * drho * dphi * dz
        flux_contribution <- irradiance_func(z_mid, rho_mid) * dV
        total_flux <- total_flux + flux_contribution
      }
    }
  }
  
  E_avg <- total_flux / V_cone
  
  return(list(
    z_height = z_height,
    E_avg = E_avg,
    V_cone = V_cone,
    rho_max = rho_max,
    total_flux = total_flux
  ))
}

z_values <- seq(0.001, 2.5, by = 0.001)
results <- data.frame(
  z = numeric(),
  E_avg = numeric(),
  V_cone = numeric(),
  rho_max = numeric(),
  total_flux = numeric()
)

total_points <- length(z_values)
progress_interval <- max(1, floor(total_points / 100))

for (i in seq_along(z_values)) {
  z <- z_values[i]
  result <- calculate_E_avg(z)
  
  results <- rbind(results, data.frame(
    z = result$z_height,
    E_avg = result$E_avg,
    V_cone = result$V_cone,
    rho_max = result$rho_max,
    total_flux = result$total_flux
  ))
  
  if (i %% progress_interval == 0 || i == total_points) {
    cat(sprintf("Progress: %d/%d (%.1f%%) - z=%.3f, E_avg=%.6f\n", 
                i, total_points, 100*i/total_points, result$z_height, result$E_avg))
  }
}

# Find z that gives chamber volume
# V = (1/3) * pi * (z*tan(theta_max))^2 * z = 31.5
# Solving: (1/3) * pi * z^3 * tan^2(theta_max) = 31.5
z_chamber <- (31.5 * 3 / (pi * tan(theta_max)^2))^(1/3)

cat(sprintf("Z for chamber volume: %.3f m\n", z_chamber))

chamber_result <- calculate_E_avg(z_chamber)
cat(sprintf("E_avg for chamber: %.6f W/m²\n", chamber_result$E_avg))

# Convert to photon flux (assuming 222 nm wavelength)
h <- 6.626e-34  # Planck constant (J⋅s)
c <- 3e8        # Speed of light (m/s)
lambda <- 222e-9 # Wavelength (m)
photon_energy <- h * c / lambda  # Energy per photon (J)

photon_flux <- chamber_result$E_avg / photon_energy  # photons/(m²⋅s)
photon_flux_cm2 <- photon_flux / 1e4  # Convert to photons/(cm²⋅s)

cat(sprintf("Photon flux: %.2e photons/(cm²⋅s)\n", photon_flux_cm2))

plot(results$z, results$E_avg, type = "l", lwd = 2,
     xlab = "Z height (m)", ylab = "E_avg (W/m²)",
     main = "Average Irradiance vs Cone Height (Fine Resolution)",
     col = "blue")
grid(TRUE)

points(z_chamber, chamber_result$E_avg, col = "red", pch = 19, cex = 1.5)
text(z_chamber, chamber_result$E_avg, 
     sprintf("Chamber\n(%.4f W/m²)", chamber_result$E_avg),
     pos = 3, col = "red")


cat(sprintf("E_avg for chamber volume (%.1f m³): %.6f W/m²\n", 
            V_chamber, chamber_result$E_avg))
cat(sprintf("Corresponding photon flux: %.2e photons/(cm²⋅s)\n", photon_flux_cm2))
cat(sprintf("Document reported E_avg: 0.032700 W/m²\n"))
cat(sprintf("Document reported photon flux: 3.65e12 photons/(cm²⋅s)\n"))
