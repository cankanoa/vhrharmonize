# preprocess-high-resolution-satellite-imagery-pipeline

Very high-resolution satellite imagery greatly enhances conservation and environmental analyses by offering broad coverage and enabling useful methodologies. However, preprocessing is complex and time-consuming because radiometric and geometric distortions pose mosaicking and harmonizing challenges. This work presents an open-source library to processes raw Worldview satellite imagery into analysis-ready, mosaicked imagery. This result is spectrally matched, allowing direct comparison between images, enhancing land cover model training, and improving analyses that rely on spectral signatures. It is also aligned with reference data which allows for highly accurate co-registration and absolute geolocation accuracy.

The pipeline begins with atmospheric correction which mitigates scattering, absorption, and adjacency effects, converting digital numbers into reflectance. Next, orthorectification removes terrain distortions, corrects sensor geometry, aligns imagery with reference data, and ensures uniform scale. Next, pansharpening fuses multispectral and panchromatic data to enhance detail. Next, global histogram matching standardizes brightness and contrast across images, assuming overlapping areas have the same spectral signature. Next, local histogram matching refines these adjustments at a block level, producing smoother transitions. Lastly, seamlines are generated between images to minimize visual artifacts and the resulting images are merged into a seamless, tiled, analysis-ready mosaic.

The mosaic is valuable for assessing land cover and land use (LCLU), detecting forest dynamics, assessing native species ranges, land managers, and other analyses of large spatial scales. The pipeline was developed and is being used to create seamless imagery of the Big Island to perform species-level forest mapping and estimate forest carbon.



# Setup
- Clone repo
- Install dependancies
- Mount drives if working on external drive


# Requirements
## Atmospheric Correction
 - ENVI >= 5.7
