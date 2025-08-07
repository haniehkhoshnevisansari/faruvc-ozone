"""
Far-UVC Lamp Average Irradiance Calculator

This module calculates the average irradiance (E_avg) for a far-uvc lamp system
Based on equations from the supplemental material of Link et al. (2023)

Author: Hanieh Khoshnevis Ansari
Date: 2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Tuple, Optional
import warnings
from dataclasses import dataclass
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class PhysicalConstants:
    """Physical constants for photon flux calculations."""
    PLANCK_CONSTANT: float = 6.626e-34  # J⋅s
    SPEED_OF_LIGHT: float = 3.0e8  # m/s
    UV_WAVELENGTH: float = 222e-9  # m (222 nm)

    @property
    def photon_energy(self) -> float:
        """Calculate energy per photon."""
        return self.PLANCK_CONSTANT * self.SPEED_OF_LIGHT / self.UV_WAVELENGTH


@dataclass
class LampParameters:
    """UV lamp model parameters from equations S3 and S4."""
    E0: float = 0.286  # W/m^(2-a) - base irradiance
    exponent_a: float = 1.52  # distance decay exponent

    # Angular distribution parameters
    angular_A: float = 0.00479
    angular_B: float = 0.9809
    angular_b: float = 0.7180
    angular_c: float = 0.1196


@dataclass
class ChamberGeometry:
    """Chamber geometric parameters."""
    volume: float = 31.5  # m³
    max_half_angle_deg: float = 55.0  # degrees
    max_height: float = 2.45  # m
    min_z: float = 0.001  # m (to avoid singularities)

    @property
    def max_half_angle_rad(self) -> float:
        """Convert max half angle to radians."""
        return np.deg2rad(self.max_half_angle_deg)


class UVLampIrradianceModel:
    """
    Model for calculating far-uvc lamp irradiance distribution and average values.

    This class implements the mathematical models from equations S3-S5 to
    calculate spatial irradiance distribution.
    """

    def __init__(self, lamp_params: LampParameters, chamber_geom: ChamberGeometry):
        self.lamp = lamp_params
        self.chamber = chamber_geom
        self.constants = PhysicalConstants()

    def calculate_irradiance(self, z: float, rho: float) -> float:
        """
        Calculate irradiance at position (z, rho) using equation S5.

        Args:
            z: Axial distance from lamp (m)
            rho: Radial distance from axis (m)

        Returns:
            Irradiance in W/m²
        """
        # Handle edge cases
        if z < self.chamber.min_z:
            return 105.0  # Maximum irradiance at lamp surface

        try:
            # Distance from lamp center
            r = np.sqrt(z ** 2 + rho ** 2)

            # Base irradiance magnitude (equation S3)
            E_magnitude = self.lamp.E0 / (r ** self.lamp.exponent_a)

            # Angular correction factor (equation S4)
            angle_ratio = rho / z
            angular_factor = (self.lamp.angular_A +
                              (self.lamp.angular_B - self.lamp.angular_A) /
                              (1 + np.exp((angle_ratio - self.lamp.angular_b) /
                                          self.lamp.angular_c)))

            result = E_magnitude * angular_factor

            # Validate result
            if not (np.isfinite(result) and result >= 0):
                logger.warning(f"Invalid irradiance at z={z:.6f}, rho={rho:.6f}: {result}")
                return 105.0

            return float(result)

        except Exception as e:
            logger.error(f"Error calculating irradiance at z={z}, rho={rho}: {e}")
            return 105.0

    def calculate_cone_volume_average(self, z_height: float,
                                      integration_points: Tuple[int, int] = (50, 30)) -> Dict[str, float]:
        """

        Args:
            z_height: Height of the cone (m)
            integration_points: Tuple of (z_points, rho_points) for numerical integration

        Returns:
            Dictionary containing calculation results
        """
        z_points, rho_points = integration_points

        # Cone geometry
        rho_max = z_height * np.tan(self.chamber.max_half_angle_rad)
        cone_volume = (np.pi / 3) * rho_max ** 2 * z_height

        # Numerical integration using Riemann sums
        total_flux = 0.0
        z_array = np.linspace(self.chamber.min_z, z_height, z_points)

        for i in range(len(z_array) - 1):
            z_mid = (z_array[i] + z_array[i + 1]) / 2
            dz = z_array[i + 1] - z_array[i]

            rho_max_at_z = z_mid * np.tan(self.chamber.max_half_angle_rad)

            if rho_max_at_z > 0:
                rho_array = np.linspace(0, rho_max_at_z, rho_points)

                for j in range(len(rho_array) - 1):
                    rho_mid = (rho_array[j] + rho_array[j + 1]) / 2
                    drho = rho_array[j + 1] - rho_array[j]

                    # Volume element in cylindrical coordinates
                    dV = rho_mid * drho * (2 * np.pi) * dz

                    # Add flux contribution
                    irradiance = self.calculate_irradiance(z_mid, rho_mid)
                    total_flux += irradiance * dV

        # Calculate average
        E_avg = total_flux / cone_volume if cone_volume > 0 else 0.0

        return {
            'z_height': z_height,
            'E_avg': E_avg,
            'cone_volume': cone_volume,
            'rho_max': rho_max,
            'total_flux': total_flux
        }

    def find_chamber_height(self) -> float:
        """
        Calculate the height that produces the target chamber volume.

        Returns:
            Height in meters
        """
        # Solve: V = (π/3) * (z*tan(θ))² * z = target_volume
        tan_theta = np.tan(self.chamber.max_half_angle_rad)
        z_chamber = ((3 * self.chamber.volume) / (np.pi * tan_theta ** 2)) ** (1 / 3)
        return z_chamber

    def convert_to_photon_flux(self, irradiance: float) -> float:
        """
        Convert irradiance to photon flux density.

        Args:
            irradiance: Irradiance in W/m²

        Returns:
            Photon flux in photons/(cm²⋅s)
        """
        photon_flux_m2 = irradiance / self.constants.photon_energy
        return photon_flux_m2 / 1e4  # Convert to cm⁻²


class UVLampAnalyzer:
    """
    High-level analyzer for UV lamp performance calculations.
    """

    def __init__(self, lamp_params: Optional[LampParameters] = None,
                 chamber_geom: Optional[ChamberGeometry] = None):
        self.lamp_params = lamp_params or LampParameters()
        self.chamber_geom = chamber_geom or ChamberGeometry()
        self.model = UVLampIrradianceModel(self.lamp_params, self.chamber_geom)

    def run_full_analysis(self, z_step: float = 0.001, z_max: float = 2.5) -> Tuple[pd.DataFrame, Dict]:
        """
        Perform complete analysis over range of heights.

        Args:
            z_step: Step size for height sweep (m)
            z_max: Maximum height to analyze (m)

        Returns:
            Tuple of (results_dataframe, chamber_result_dict)
        """
        logger.info("Starting UV lamp irradiance analysis")
        logger.info(f"Height range: {self.chamber_geom.min_z:.3f} to {z_max:.3f} m (step: {z_step:.3f})")

        z_values = np.arange(self.chamber_geom.min_z, z_max + z_step, z_step)
        results = []

        total_points = len(z_values)
        progress_interval = max(1, total_points // 20)  # 20 progress updates

        for i, z in enumerate(z_values):
            result = self.model.calculate_cone_volume_average(z)
            results.append(result)

            if i % progress_interval == 0 or i == total_points - 1:
                progress = 100 * (i + 1) / total_points
                logger.info(f"Progress: {i + 1:4d}/{total_points} ({progress:5.1f}%) - "
                            f"z={z:.3f}m, E_avg={result['E_avg']:.6f} W/m²")

        df_results = pd.DataFrame(results)

        logger.info(f"Calculating for target chamber volume: {self.chamber_geom.volume} m³")
        z_chamber = self.model.find_chamber_height()
        chamber_result = self.model.calculate_cone_volume_average(z_chamber)

        logger.info(f"Chamber height: {z_chamber:.3f} m")
        logger.info(f"Chamber E_avg: {chamber_result['E_avg']:.6f} W/m²")

        return df_results, chamber_result

    def create_visualization(self, df_results: pd.DataFrame, chamber_result: Dict) -> plt.Figure:
        """

        Args:
            df_results: DataFrame with height sweep results
            chamber_result: Dictionary with chamber-specific results

        Returns:
            Matplotlib figure object
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

        # Main irradiance plot
        ax1.plot(df_results['z_height'], df_results['E_avg'], 'b-', linewidth=2,
                 label='Average Irradiance')
        ax1.scatter(chamber_result['z_height'], chamber_result['E_avg'],
                    color='red', s=100, zorder=5, label='Target Chamber')

        ax1.set_xlabel('Cone Height (m)')
        ax1.set_ylabel('Average Irradiance (W/m²)')
        ax1.set_title('UV Lamp Average Irradiance vs Cone Height')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Volume plot
        ax2.plot(df_results['z_height'], df_results['cone_volume'], 'g-', linewidth=2,
                 label='Cone Volume')
        ax2.axhline(y=self.chamber_geom.volume, color='red', linestyle='--',
                    label=f'Target Volume ({self.chamber_geom.volume} m³)')

        ax2.set_xlabel('Cone Height (m)')
        ax2.set_ylabel('Volume (m³)')
        ax2.set_title('Cone Volume vs Height')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        plt.tight_layout()
        return fig

    def generate_report(self, chamber_result: Dict) -> str:
        """
        Generate a formatted analysis report.

        Args:
            chamber_result: Dictionary with chamber calculation results

        Returns:
            Formatted report string
        """
        photon_flux = self.model.convert_to_photon_flux(chamber_result['E_avg'])

        report = f"""
UV LAMP IRRADIANCE ANALYSIS REPORT
{'=' * 50}

CHAMBER PARAMETERS:
  Volume: {self.chamber_geom.volume:.1f} m³
  Max Half-Angle: {self.chamber_geom.max_half_angle_deg:.1f}°
  Required Height: {chamber_result['z_height']:.3f} m
  Max Radius: {chamber_result['rho_max']:.3f} m

LAMP PARAMETERS:
  Base Irradiance (E₀): {self.lamp_params.E0:.3f} W/m^(2-a)
  Distance Exponent (a): {self.lamp_params.exponent_a:.2f}
  Angular Parameters: A={self.lamp_params.angular_A:.5f}, B={self.lamp_params.angular_B:.4f}

RESULTS:
  Average Irradiance: {chamber_result['E_avg']:.6f} W/m²
  Total Flux: {chamber_result['total_flux']:.2e} W
  Photon Flux Density: {photon_flux:.2e} photons/(cm²⋅s)

REFERENCE VALUES:
  Document E_avg: 0.032700 W/m²
  Document Photon Flux: 3.65e12 photons/(cm²⋅s)

CALCULATION NOTES:
  - Used cylindrical coordinate numerical integration
  - Minimum z-value: {self.chamber_geom.min_z:.3f} m (singularity avoidance)
  - Integration points: 50×30 (z×ρ)
"""
        return report


def main():
    """Main execution function."""

    # Initialize analyzer with default parameters
    analyzer = UVLampAnalyzer()

    # Run full analysis
    results_df, chamber_result = analyzer.run_full_analysis(z_step=0.001, z_max=2.5)

    # Create visualization
    fig = analyzer.create_visualization(results_df, chamber_result)
    plt.show()

    # Generate and display report
    report = analyzer.generate_report(chamber_result)
    print(report)

    # Save results
    output_file = "uv_lamp_analysis_results.csv"
    results_df.to_csv(output_file, index=False)
    logger.info(f"Results saved to '{output_file}'")

    return results_df, chamber_result


if __name__ == "__main__":
    results_df, chamber_result = main()
