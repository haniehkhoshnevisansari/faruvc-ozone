# Far-UVC Lamp E_avg Calculation for Variable Theta and Z Values
# Extension to analyze E_avg as a function of both cone angle (theta) and height (z)

library(pracma) 

E0 <- 0.286  # W/m^(2-a) from equation S3
a <- 1.52    
A <- 0.00479
B <- 0.9809
b <- 0.7180
c <- 0.1196

z_min <- 0.001

irradiance_func <- function(z, rho) {
  if (z < z_min) {
    return(105)  
  }
  
  r <- sqrt(z^2 + rho^2)
  
  E_mag <- E0 / (r^a)
  
  if (z <= 0) {
    return(105)  
  }
  
  angle_ratio <- rho / z
  angular_factor <- A + (B - A) / (1 + exp((angle_ratio - b) / c))
  
  result <- E_mag * angular_factor
  
  if (result < 0 || is.na(result) || !is.finite(result)) {
    return(105)  
  }
  
  return(result)
}

# Function to calculate E_avg for a given z height and theta
calculate_E_avg_theta_z <- function(z_height, theta_degrees) {
  theta_rad <- theta_degrees * pi / 180  # Convert to radians
  
  rho_max <- z_height * tan(theta_rad)
  
  V_cone <- (1/3) * pi * rho_max^2 * z_height
  
  total_flux <- 0

  z_points <- seq(z_min, z_height, length.out = 50)
  
  for (i in 1:(length(z_points)-1)) {
    z_low <- z_points[i]
    z_high <- z_points[i+1]
    z_mid <- (z_low + z_high) / 2
    
    rho_max_at_z <- z_mid * tan(theta_rad)
    
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
    theta_degrees = theta_degrees,
    theta_radians = theta_rad,
    E_avg = E_avg,
    V_cone = V_cone,
    rho_max = rho_max,
    total_flux = total_flux
  ))
}


z_values <- seq(0.001, 2.5, by = 0.01)  # Same range as first code but with larger steps for efficiency
theta_values <- seq(10, 80, by = 1)     # Every 1 degree from 10° to 80°

results_matrix <- data.frame(
  z = numeric(),
  theta_degrees = numeric(),
  theta_radians = numeric(),
  E_avg = numeric(),
  V_cone = numeric(),
  rho_max = numeric()
)

cat(sprintf("Total calculations: %d x %d = %d\n", 
            length(z_values), length(theta_values), 
            length(z_values) * length(theta_values)))

total_calculations <- length(z_values) * length(theta_values)
current_calc <- 0
progress_interval <- max(1, floor(total_calculations / 100))  # Show progress every 1%

# Main calculation loop
start_time <- Sys.time()

for (theta in theta_values) {
  
  for (z in z_values) {
    current_calc <- current_calc + 1
    
    result <- calculate_E_avg_theta_z(z, theta)
    
    results_matrix <- rbind(results_matrix, data.frame(
      z = result$z_height,
      theta_degrees = result$theta_degrees,
      theta_radians = result$theta_radians,
      E_avg = result$E_avg,
      V_cone = result$V_cone,
      rho_max = result$rho_max
    ))
    
    if (current_calc %% progress_interval == 0 || current_calc == total_calculations) {
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      estimated_total <- elapsed_time * total_calculations / current_calc
      remaining_time <- estimated_total - elapsed_time
      
      cat(sprintf("Progress: %d/%d (%.1f%%) - Elapsed: %.1f min, Est. remaining: %.1f min\n", 
                  current_calc, total_calculations, 100*current_calc/total_calculations,
                  elapsed_time, remaining_time))
    }
  }
}

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

theta_summary <- data.frame(
  theta_degrees = numeric(),
  min_E_avg = numeric(),
  max_E_avg = numeric(),
  mean_E_avg = numeric(),
  E_avg_at_z_1m = numeric(),
  E_avg_at_z_2m = numeric(),
  E_avg_at_z_2_45m = numeric()
)

for (theta in theta_values) {
  theta_data <- results_matrix[results_matrix$theta_degrees == theta, ]
  
  E_avg_1m <- theta_data$E_avg[which.min(abs(theta_data$z - 1.0))]
  E_avg_2m <- theta_data$E_avg[which.min(abs(theta_data$z - 2.0))]
  E_avg_2_45m <- theta_data$E_avg[which.min(abs(theta_data$z - 2.45))]
  
  if (length(E_avg_1m) == 0) E_avg_1m <- NA
  if (length(E_avg_2m) == 0) E_avg_2m <- NA
  if (length(E_avg_2_45m) == 0) E_avg_2_45m <- NA
  
  theta_summary <- rbind(theta_summary, data.frame(
    theta_degrees = theta,
    min_E_avg = min(theta_data$E_avg, na.rm = TRUE),
    max_E_avg = max(theta_data$E_avg, na.rm = TRUE),
    mean_E_avg = mean(theta_data$E_avg, na.rm = TRUE),
    E_avg_at_z_1m = E_avg_1m,
    E_avg_at_z_2m = E_avg_2m,
    E_avg_at_z_2_45m = E_avg_2_45m
  ))
  
  cat(sprintf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
              theta, min(theta_data$E_avg), max(theta_data$E_avg), 
              mean(theta_data$E_avg), E_avg_1m, E_avg_2m, E_avg_2_45m))
}

par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

# E_avg vs Z for different theta values (show selected theta values)
selected_thetas <- c(15, 25, 35, 45, 55, 65, 75)
colors <- rainbow(length(selected_thetas))

plot(0, 0, type = "n", xlim = c(0, 2.5), ylim = c(0, max(results_matrix$E_avg)),
     xlab = "Z height (m)", ylab = "E_avg (W/m²)",
     main = "E_avg vs Z for Different Theta Values")

for (i in seq_along(selected_thetas)) {
  theta <- selected_thetas[i]
  theta_data <- results_matrix[results_matrix$theta_degrees == theta, ]
  lines(theta_data$z, theta_data$E_avg, col = colors[i], lwd = 2)
}

legend("topright", legend = paste0(selected_thetas, "°"), 
       col = colors, lwd = 2, cex = 0.8, title = "Theta")
grid(TRUE)

theta_55_data <- results_matrix[results_matrix$theta_degrees == 55, ]
lines(theta_55_data$z, theta_55_data$E_avg, col = "red", lwd = 3)

# E_avg at z=2.45m vs Theta
plot(theta_summary$theta_degrees, theta_summary$E_avg_at_z_2_45m, 
     type = "b", lwd = 2, pch = 16, col = "blue",
     xlab = "Theta (degrees)", ylab = "E_avg at z=2.45m (W/m²)",
     main = "E_avg at z=2.45m vs Theta")
grid(TRUE)

max_theta_245_idx <- which.max(theta_summary$E_avg_at_z_2_45m)
points(theta_summary$theta_degrees[max_theta_245_idx], 
       theta_summary$E_avg_at_z_2_45m[max_theta_245_idx], 
       col = "red", pch = 19, cex = 1.5)
text(theta_summary$theta_degrees[max_theta_245_idx], 
     theta_summary$E_avg_at_z_2_45m[max_theta_245_idx],
     sprintf("Max\n(%.0f°)", theta_summary$theta_degrees[max_theta_245_idx]),
     pos = 3, col = "red")


cat(sprintf("Total data points calculated: %d\n", nrow(results_matrix)))
cat(sprintf("Z range: %.3f to %.3f m\n", min(results_matrix$z), max(results_matrix$z)))
cat(sprintf("Theta range: %.0f° to %.0f°\n", 
            min(results_matrix$theta_degrees), max(results_matrix$theta_degrees)))
cat(sprintf("E_avg range: %.6f to %.6f W/m²\n", 
            min(results_matrix$E_avg), max(results_matrix$E_avg)))

global_max_idx <- which.max(results_matrix$E_avg)
cat(sprintf("Global maximum E_avg: %.6f W/m² at z=%.3f m, theta=%.0f°\n",
            results_matrix$E_avg[global_max_idx],
            results_matrix$z[global_max_idx],
            results_matrix$theta_degrees[global_max_idx]))

original_idx <- which(abs(results_matrix$z - 2.45) < 0.02 & 
                     results_matrix$theta_degrees == 55)
if (length(original_idx) > 0) {
  cat(sprintf("Original case (55°, z=2.45m): %.6f W/m²\n",
              results_matrix$E_avg[original_idx[1]]))
}

write.csv(results_matrix, "eavg_theta_z_matrix.csv", row.names = FALSE)
write.csv(theta_summary, "eavg_theta_summary.csv", row.names = FALSE)

par(mfrow = c(1, 1))
