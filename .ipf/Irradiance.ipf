#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Far-UVC Lamp E_avg Calculation 
// based on the method described by Link et al. (2023)

// Global constants
constant E0 = 0.286		// W/m^(2-a) from equation S3
constant a_exp = 1.52		// exponent from equation S3

// Angular distribution parameters from equation S4
constant A_param = 0.00479
constant B_param = 0.9809
constant beta_param = 0.7180
constant gamma_param = 0.1196

// Chamber parameters
constant V_chamber = 31.5		// m^3 - total chamber volume
constant theta_max_deg = 55		// Maximum half-angle in degrees
constant theta_max = 0.9599		// Maximum half-angle in radians (55 * pi / 180)
constant h_max = 2.45			// Maximum height of cone (m)

// Function to calculate irradiance at given z and rho coordinates
Function irradiance_func_calculated(z, rho)
	Variable z, rho
	
	if (z <= 0)
		return NaN
	endif
	
	// Calculate distance from lamp
	Variable r = sqrt(z^2 + rho^2)
	
	// Calculate irradiance magnitude
	Variable E_mag = E0 / (r^a_exp)
	
	// Calculate angular factor
	Variable angle_ratio = rho / z
	Variable angular_factor = A_param + (B_param - A_param) / (1 + exp((angle_ratio - beta_param) / gamma_param))
	
	Variable result = E_mag * angular_factor
	
	if (result < 0 || numtype(result) != 0)
		return NaN
	endif
	
	return result
End

// Function to calculate E_avg for a given z value using only calculated values
Function calculate_E_avg_no_limit(z_height, E_avg_out, V_cone_out, rho_max_out, total_flux_out)
	Variable z_height
	Variable &E_avg_out, &V_cone_out, &rho_max_out, &total_flux_out
	
	if (z_height <= 0)
		E_avg_out = NaN
		V_cone_out = NaN
		rho_max_out = NaN
		total_flux_out = NaN
		return -1
	endif
	
	rho_max_out = z_height * tan(theta_max)
	
	// Calculate volume of cone up to height z
	V_cone_out = (1/3) * pi * rho_max_out^2 * z_height
	
	Variable total_flux = 0
	Variable valid_integration = 1
	
	Variable z_start = 1e-6
	
	// Create z points for integration
	Variable num_z_points = 50
	Variable i, j
	Variable z_low, z_high, z_mid
	Variable rho_low, rho_high, rho_mid
	Variable rho_max_at_z
	Variable dz, drho, dphi, dV
	Variable irradiance_val, flux_contribution
	
	for (i = 0; i < num_z_points - 1; i += 1)
		z_low = z_start + i * (z_height - z_start) / (num_z_points - 1)
		z_high = z_start + (i + 1) * (z_height - z_start) / (num_z_points - 1)
		z_mid = (z_low + z_high) / 2
		
		rho_max_at_z = z_mid * tan(theta_max)
		
		if (rho_max_at_z > 0)
			Variable num_rho_points = 30
			
			for (j = 0; j < num_rho_points - 1; j += 1)
				rho_low = j * rho_max_at_z / (num_rho_points - 1)
				rho_high = (j + 1) * rho_max_at_z / (num_rho_points - 1)
				rho_mid = (rho_low + rho_high) / 2
				
				dz = z_high - z_low
				drho = rho_high - rho_low
				dphi = 2 * pi
				
				dV = rho_mid * drho * dphi * dz
				irradiance_val = irradiance_func_calculated(z_mid, rho_mid)
				
				if (numtype(irradiance_val) != 0)
					valid_integration = 0
					break
				endif
				
				flux_contribution = irradiance_val * dV
				total_flux += flux_contribution
			endfor
			
			if (valid_integration == 0)
				break
			endif
		endif
	endfor
	
	if (valid_integration == 0 || V_cone_out <= 0)
		E_avg_out = NaN
		total_flux_out = NaN
		return -1
	endif
	
	// Calculate average irradiance
	E_avg_out = total_flux / V_cone_out
	total_flux_out = total_flux
	
	return 0
End

// Function to find z value that gives E_avg = 105 W/m²
Function find_z_min(target_E_avg)
	Variable target_E_avg
	
	print "Searching for appropriate range to find z_min..."
	
	// Test different z values to find appropriate range
	Make/O/N=8 test_z_values = {0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2}
	Make/O/N=8 test_E_avg_values
	Make/O/N=8 test_valid
	
	Variable i
	Variable E_avg_temp, V_cone_temp, rho_max_temp, total_flux_temp
	
	for (i = 0; i < numpnts(test_z_values); i += 1)
		Variable result = calculate_E_avg_no_limit(test_z_values[i], E_avg_temp, V_cone_temp, rho_max_temp, total_flux_temp)
		
		if (result == 0)
			test_E_avg_values[i] = E_avg_temp
			test_valid[i] = 1
			printf "z = %.4f m: E_avg = %.2f W/m²\r", test_z_values[i], E_avg_temp
		else
			test_E_avg_values[i] = NaN
			test_valid[i] = 0
		endif
	endfor
	
	// Find bounds where E_avg crosses target value
	Variable z_lower, z_upper
	Variable found_bounds = 0
	
	// Find crossing point
	for (i = 0; i < numpnts(test_z_values) - 1; i += 1)
		if (test_valid[i] & test_valid[i+1])
			Variable crosses_up = (test_E_avg_values[i] < target_E_avg) & (test_E_avg_values[i+1] > target_E_avg)
			Variable crosses_down = (test_E_avg_values[i] > target_E_avg) & (test_E_avg_values[i+1] < target_E_avg)
			if (crosses_up | crosses_down)
				z_lower = min(test_z_values[i], test_z_values[i+1])
				z_upper = max(test_z_values[i], test_z_values[i+1])
				found_bounds = 1
				break
			endif
		endif
	endfor
	
	if (found_bounds == 0)
		// Use default bounds
		z_lower = 0.001
		z_upper = 0.1
		print "Using default search range"
	endif
	
	printf "Search range: [%.6f, %.6f] m\r", z_lower, z_upper
	
	// Use binary search to find z_min
	Variable tolerance = 1e-8
	Variable z_mid, E_avg_mid
	Variable iterations = 0
	Variable max_iterations = 100
	
	do
		z_mid = (z_lower + z_upper) / 2
		calculate_E_avg_no_limit(z_mid, E_avg_mid, V_cone_temp, rho_max_temp, total_flux_temp)
		
		if (abs(E_avg_mid - target_E_avg) < tolerance)
			break
		endif
		
		if (E_avg_mid > target_E_avg)
			z_lower = z_mid
		else
			z_upper = z_mid
		endif
		
		iterations += 1
	while (iterations < max_iterations && (z_upper - z_lower) > tolerance)
	
	return z_mid
End

// Piecewise irradiance function
Function irradiance_func_piecewise(z, rho, z_min_val)
	Variable z, rho, z_min_val
	
	if (z < z_min_val)
		return 105  // W/m^2 - as defined in equation S5
	else
		return irradiance_func_calculated(z, rho)
	endif
End

// Function to calculate E_avg using piecewise function
Function calculate_E_avg_piecewise(z_height, z_min_val, E_avg_out, V_cone_out, rho_max_out, total_flux_out, has_forced_region_out)
	Variable z_height, z_min_val
	Variable &E_avg_out, &V_cone_out, &rho_max_out, &total_flux_out, &has_forced_region_out
	
	if (z_height <= 0)
		E_avg_out = NaN
		V_cone_out = NaN
		rho_max_out = NaN
		total_flux_out = NaN
		has_forced_region_out = 0
		return -1
	endif
	
	rho_max_out = z_height * tan(theta_max)
	V_cone_out = (1/3) * pi * rho_max_out^2 * z_height
	has_forced_region_out = z_height > z_min_val
	
	Variable total_flux = 0
	Variable z_start = 1e-6
	
	// Integration parameters
	Variable num_z_points = 50
	Variable num_rho_points = 30
	Variable i, j
	Variable z_low, z_high, z_mid
	Variable rho_low, rho_high, rho_mid
	Variable rho_max_at_z
	Variable dz, drho, dphi, dV
	Variable irradiance_val, flux_contribution
	
	for (i = 0; i < num_z_points - 1; i += 1)
		z_low = z_start + i * (z_height - z_start) / (num_z_points - 1)
		z_high = z_start + (i + 1) * (z_height - z_start) / (num_z_points - 1)
		z_mid = (z_low + z_high) / 2
		
		rho_max_at_z = z_mid * tan(theta_max)
		
		if (rho_max_at_z > 0)
			for (j = 0; j < num_rho_points - 1; j += 1)
				rho_low = j * rho_max_at_z / (num_rho_points - 1)
				rho_high = (j + 1) * rho_max_at_z / (num_rho_points - 1)
				rho_mid = (rho_low + rho_high) / 2
				
				dz = z_high - z_low
				drho = rho_high - rho_low
				dphi = 2 * pi
				
				dV = rho_mid * drho * dphi * dz
				irradiance_val = irradiance_func_piecewise(z_mid, rho_mid, z_min_val)
				
				flux_contribution = irradiance_val * dV
				total_flux += flux_contribution
			endfor
		endif
	endfor
	
	E_avg_out = total_flux / V_cone_out
	total_flux_out = total_flux
	
	return 0
End

// Main function to run the calculation
Function RunUVLampCalculation()
	
	print "=== UV Lamp E_avg Calculation ==="
	print "Finding z_min where E_avg = 105 W/m²..."
	
	// Find z_min
	Variable z_min = find_z_min(105)
	Variable E_avg_at_z_min, V_cone_at_z_min, rho_max_at_z_min, total_flux_at_z_min
	calculate_E_avg_no_limit(z_min, E_avg_at_z_min, V_cone_at_z_min, rho_max_at_z_min, total_flux_at_z_min)
	
	printf "\r=== FOUND z_min ===\r"
	printf "z_min = %.6f m\r", z_min
	printf "E_avg at z_min = %.6f W/m²\r", E_avg_at_z_min
	printf "Volume at z_min = %.6f m³\r", V_cone_at_z_min
	
	// Calculate E_avg for different z values using piecewise function
	Variable z_start = 0.001
	Variable z_end = 2.5
	Variable z_step = 0.01
	Variable num_points = floor((z_end - z_start) / z_step) + 1
	
	Make/O/N=(num_points) z_values, E_avg_values, V_cone_values, rho_max_values, total_flux_values, has_forced_region
	
	printf "\rCalculating E_avg using piecewise function with z_min = %.6f m...\r", z_min
	
	Variable i
	Variable E_avg_temp, V_cone_temp, rho_max_temp, total_flux_temp, has_forced_temp
	
	for (i = 0; i < num_points; i += 1)
		Variable z_val = z_start + i * z_step
		z_values[i] = z_val
		
		calculate_E_avg_piecewise(z_val, z_min, E_avg_temp, V_cone_temp, rho_max_temp, total_flux_temp, has_forced_temp)
		
		E_avg_values[i] = E_avg_temp
		V_cone_values[i] = V_cone_temp
		rho_max_values[i] = rho_max_temp
		total_flux_values[i] = total_flux_temp
		has_forced_region[i] = has_forced_temp
		
		if (mod(i, 50) == 0)
			printf "Progress: %d/%d - z=%.3f, E_avg=%.3f\r", i+1, num_points, z_val, E_avg_temp
		endif
	endfor
	
	// Calculate for chamber volume
	print "\rCalculating E_avg for full chamber volume (31.5 m³)..."
	
	Variable z_chamber = (V_chamber * 3 / (pi * tan(theta_max)^2))^(1/3)
	printf "Z for chamber volume: %.3f m\r", z_chamber
	
	Variable chamber_E_avg, chamber_V_cone, chamber_rho_max, chamber_total_flux, chamber_forced
	calculate_E_avg_piecewise(z_chamber, z_min, chamber_E_avg, chamber_V_cone, chamber_rho_max, chamber_total_flux, chamber_forced)
	printf "E_avg for chamber: %.6f W/m²\r", chamber_E_avg
	
	// Convert to photon flux (assuming 222 nm wavelength)
	Variable h_planck = 6.626e-34		// Planck constant (J⋅s)
	Variable c_light = 3e8			// Speed of light (m/s)
	Variable lambda_nm = 222e-9		// Wavelength (m)
	Variable photon_energy = h_planck * c_light / lambda_nm	// Energy per photon (J)
	
	Variable photon_flux = chamber_E_avg / photon_energy		// photons/(m²⋅s)
	Variable photon_flux_cm2 = photon_flux / 1e4			// Convert to photons/(cm²⋅s)
	
	printf "Photon flux: %.2e photons/(cm²⋅s)\r", photon_flux_cm2
	
	// Create plots
	Display E_avg_values vs z_values
	ModifyGraph mode=0, lsize=2, rgb=(0,0,65535)
	Label left "E_avg (W/m²)"
	Label bottom "Z height (m)"
	SetAxis left 0,*
	SetAxis bottom 0,2.5
	
	// Add vertical line at z_min
	Make/O/N=2 z_min_line = z_min
	Make/O/N=2 y_line_for_z_min = {0, 200}
	AppendToGraph y_line_for_z_min vs z_min_line
	ModifyGraph mode[1]=0, lsize[1]=2, rgb[1]=(65535,0,0), lstyle[1]=2
	
	// Add point at z_min
	Make/O/N=1 z_min_point = z_min
	Make/O/N=1 E_avg_point = 105
	AppendToGraph E_avg_point vs z_min_point
	ModifyGraph mode[2]=3, marker[2]=19, msize[2]=4, rgb[2]=(65535,0,0)
	
	// Add chamber point
	Make/O/N=1 z_chamber_point = z_chamber
	Make/O/N=1 chamber_E_avg_point = chamber_E_avg
	AppendToGraph chamber_E_avg_point vs z_chamber_point
	ModifyGraph mode[3]=3, marker[3]=19, msize[3]=4, rgb[3]=(0,0,0)
	
	TextBox/C/N=text0/F=0/A=RT "z_min = " + num2str(z_min) + " m\\rE_avg = 105 W/m²"
	
	// Volume plot
	Display V_cone_values vs z_values
	ModifyGraph mode=0, lsize=2, rgb=(32768,0,32768)
	Label left "Volume (m³)"
	Label bottom "Z height (m)"
	
	// Add horizontal line at chamber volume
	Make/O/N=2 x_line_for_volume = {0, 2.5}
	Make/O/N=2 chamber_volume_line = V_chamber
	AppendToGraph chamber_volume_line vs x_line_for_volume
	ModifyGraph mode[1]=0, lsize[1]=2, rgb[1]=(65535,0,0), lstyle[1]=2
	
	// Add vertical line at z_chamber
	Make/O/N=2 z_chamber_line = z_chamber
	Make/O/N=2 y_line_for_volume = {0, 40}
	AppendToGraph y_line_for_volume vs z_chamber_line
	ModifyGraph mode[2]=0, lsize[2]=2, rgb[2]=(65535,0,0), lstyle[2]=2
	
	// Add chamber point
	Make/O/N=1 z_chamber_vol_point = z_chamber
	Make/O/N=1 V_chamber_point = V_chamber
	AppendToGraph V_chamber_point vs z_chamber_vol_point
	ModifyGraph mode[3]=3, marker[3]=19, msize[3]=4, rgb[3]=(65535,0,0)
	
	// Print summary
	print "\r=== SUMMARY ==="
	printf "Found z_min = %.6f m (where E_avg = 105 W/m²)\r", z_min
	printf "Volume at z_min = %.6f m³\r", V_cone_at_z_min
	printf "Chamber z value: %.6f m\r", z_chamber
	printf "E_avg for chamber volume (%.1f m³): %.6f W/m²\r", V_chamber, chamber_E_avg
	printf "Corresponding photon flux: %.2e photons/(cm²⋅s)\r", photon_flux_cm2
	printf "Paper reported E_avg: 0.032700 W/m²\r"
	printf "Paper reported photon flux: 3.65e12 photons/(cm²⋅s)\r"
	printf "Ratio (calculated/paper): %.2f\r", chamber_E_avg / 0.0327
	
	// Store key results in global variables
	Variable/G g_z_min = z_min
	Variable/G g_chamber_E_avg = chamber_E_avg
	Variable/G g_photon_flux_cm2 = photon_flux_cm2
	Variable/G g_z_chamber = z_chamber
	
	print "Calculation complete. Key results stored in global variables."
	print "Run 'RunUVLampCalculation()' to execute this analysis."
End
