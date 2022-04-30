# DORMOUSE: Detection Of Reflected Microscopic Optical UltraSound Emission
**by Kenton Kwok, Yuji Go**

part of Year 3 BSc Project, Department of Physics, Imperial College London

This is the GitHub repository of our BSc project under the supervision of Dr Chris Rowlands (Bioengineering, Imperial College London). We use a computer simulation approach to improve the understanding of a problem in an optics laboratory. 

## Projection Mechanism
Mainly worked on by Kenton

This is a k-Wave simulation to investigate the feasibility of a particular setup used for all-optical ultrasound generation. In particular we are concerned with a phased-array implementation of ultrasound scanning.


Two k-Wave simulations exist and are in MATLAB code. Running them requires the k-Wave simulation files which can be downloaded [here](http://www.k-wave.org/).

* `simulations/twod_scan.m` which looks at the degradation of resolution
* `simulations/twod_phantom.m` which looks at the image performance when the waves are directed into a scattering phantom

They are based on the same example code by Treeby et. al. We could not use the default transducer MATLAB object provided as it did not give us the flexibility to modify the time delays. The result was to calculate the input signals and delay them externally. 

The reset of the files are to do with analysis and plotting graphs
* `analysis/basic.ipynb` concerns simple plots for the project report
* `analysis/image.ipynb` concerns making plots of the phantoms
*  `analysis/analysis.ipynb` is the main file and concerns the investigation of resolution and plotting figures

The analysis files will refer to raw `.mat` data files obtained using MATLAB but they have not been uploaded due to large file size.

The `phantoms` folder contains 
* the ultrasound images used 
* ImageJ ROI `.zip` files for annotation and analysis
* Excel files for analysis of the contrast 

Two gifs are also uploaded in the `videos` folder for illustration

## Detection Mechanism
Mainly worked on by Yuji
