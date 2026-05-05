# Py6S Overview
> Credit to truegeometry.com for this concise article that has subsequently been removed from the internet.

**1. Underlying Principles: Radiative Transfer & Multiple Scattering**

**Radiative Transfer Equation (RTE):** At its core, 6S is based on the RTE. This equation describes how light (electromagnetic radiation) propagates through a medium (in this case, the atmosphere). It accounts for absorption, scattering, and emission of light as it travels. The RTE is complex, involving integrals over all directions and wavelengths.
Scattering is Key: Atmospheric correction is primarily about dealing with scattering. Sunlight interacts with atmospheric particles (molecules, aerosols, clouds) and is scattered in various directions. This scattered light reaches the sensor, contaminating the signal that originated from the surface being observed.
**Multiple Scattering:** A crucial aspect is multiple scattering. A photon might be scattered once, then scattered again, and again, before reaching the sensor. Ignoring multiple scattering leads to significant errors in atmospheric correction. This is where 6S's "six-stream" approach comes in.

**2. Components of the 6S Model**

6S simplifies the RTE by using a stream approach. Here's a breakdown of the key components:

**Atmospheric Layers:** The atmosphere is divided into a series of discrete layers. The number of layers is a user-defined parameter (typically 10-30 layers). Each layer has its own properties.
Streams: The "6S" name comes from the fact that the model uses six streams of radiation. These streams represent different directions of propagation:
**Two Upward Streams:** Representing radiation traveling upwards from the surface.
**Two Downward Streams:** Representing radiation traveling downwards towards the sensor.
**Two Horizontal Streams:** Representing radiation traveling horizontally. These are important for accounting for scattering between layers.
**Optical Properties:** Each layer is characterized by its optical properties:
**Absorption Coefficient (α):** How much light is absorbed by the gases and particles in the layer. This depends on the wavelength and the composition of the atmosphere (e.g., ozone, water vapor, carbon dioxide, aerosols).
**Scattering Coefficient (σ):** How much light is scattered by the particles in the layer. This is wavelength-dependent and depends on the size, shape, and refractive index of the particles.
**Asymmetry Factor (g):** Describes the directionality of scattering. g = 0 means isotropic scattering (equal in all directions). g = 1 means forward scattering. g = -1 means backward scattering.
**Single Scattering Albedo (ω):** The ratio of scattering to extinction (scattering + absorption). ω = σ / (α + σ). It indicates how much light is scattered relative to how much is absorbed.
**Surface Reflectance (ρs):** The reflectance of the surface being observed. This is often an unknown that the model tries to estimate.
**Sensor Geometry:** The solar zenith angle (θs), viewing zenith angle (θv), and relative azimuth angle (φ) are crucial inputs. These define the angles between the sun, the surface, and the sensor.

**3. The Calculation Process (Simplified)**

The 6S model iteratively solves a simplified version of the RTE using the stream approach. Here's a high-level overview:

**Initialization:** The model starts with initial guesses for the stream intensities (the amount of radiation in each stream). Often, the streams are initialized with the direct beam radiation.
**Iteration:** The model then iterates through the atmospheric layers, performing the following steps for each layer and each stream:
**Calculate Extinction:** Determine how much radiation is lost due to absorption and scattering within the layer.
**Calculate Scattering:** Calculate how much radiation is scattered into the current stream from other streams (including the streams in adjacent layers). This is the most computationally intensive part, as it involves multiple scattering calculations. The asymmetry factor (g) is critical here.
**Update Stream Intensity:** Update the intensity of each stream based on the extinction and scattering calculations.
**Surface Interaction:** When the streams reach the surface, the surface reflectance (ρs) is applied. The streams are split into upward and downward components based on the viewing geometry.
**Sensor Measurement:** The downward streams that reach the sensor are combined (weighted by their respective areas on the sensor) to simulate the sensor's measurement.
**Convergence:** The iteration continues until the stream intensities converge (i.e., they change very little between iterations). This indicates that the model has reached a stable solution.
**Atmospheric Transmittance & Path Radiance:** From the converged stream intensities, the model calculates:
**Atmospheric Transmittance (T):** The fraction of direct solar radiation that reaches the surface without being scattered or absorbed.
**Path Radiance (L):** The radiance of the sky (scattered light) that reaches the sensor.

**4. Atmospheric Correction using 6S**

Once the atmospheric transmittance (T) and path radiance (L) are calculated, they can be used to correct the sensor data:

**Simple Correction (Radiance-based):**

ρc = (ρs * T) / (1 - L/Es)
ρc: Corrected surface reflectance
ρs: Measured surface reflectance (from the sensor data)
T: Atmospheric transmittance
L: Path radiance
Es: Extraterrestrial solar irradiance (the amount of solar energy reaching the top of the atmosphere)

**More Complex Approaches:** 6S can also be used in more sophisticated atmospheric correction schemes that involve iterative optimization to estimate surface reflectance.

**5. Limitations of 6S**

**Cloud Handling:** 6S is not designed to handle clouds directly. Clouds are treated as opaque areas, and the model typically requires a cloud mask to be applied to the data before atmospheric correction. There are extensions to 6S that attempt to account for thin clouds, but they are not as robust as dedicated cloud correction algorithms.
**Aerosol Characterization:** Accurate aerosol characterization (size, shape, composition) is crucial for 6S. The model often relies on pre-defined aerosol models or requires external aerosol data (e.g., from AERONET). Errors in aerosol characterization can lead to significant errors in atmospheric correction.
**Surface BRDF:** 6S assumes a relatively simple surface reflectance model. It doesn't explicitly account for the Bidirectional Reflectance Distribution Function (BRDF) of the surface, which describes how reflectance varies with viewing and illumination geometry. This can introduce errors, especially for surfaces with complex BRDFs (e.g., vegetation).
**Rayleigh Scattering:** While 6S handles Rayleigh scattering (scattering by air molecules) well, it may not be as accurate for very high wavelengths (e.g., near-infrared) where the assumptions about Rayleigh scattering may not hold.

Resources for Further Learning:

6S Documentation: https://www.spectral.co.uk/6S (This is the official website and has detailed documentation.)
Numerous research papers: Search for "6S atmospheric correction" on Google Scholar or Web of Science.
