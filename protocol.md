    Authors: Matthew D. Lew,  Alexander R. S. von Diezmann & W. E. Moerner


### Abstract

Automated processing of double-helix (DH) microscope images of single molecules (SMs) streamlines the protocol required to obtain super-resolved three-dimensional (3D) reconstructions of ultrastructures in biological samples by single-molecule active control microscopy. Here, we present a [suite of MATLAB subroutines](https://code.google.com/p/easy-dhpsf/), bundled with an easy-to-use [graphical user interface](http://www.nature.com/protocolexchange/protocols/2622#f0), that facilitates 3D localization of single emitters (e.g. SMs, fluorescent beads, or quantum dots) with precisions of tens of nanometers in multi-frame movies acquired using a wide-field DH epifluorescence microscope. The algorithmic approach is based upon template matching for SM recognition and least-squares fitting for 3D position measurement, both of which are computationally expedient and precise. Overlapping images of SMs are ignored, and the precision of least-squares fitting is not as high as maximum likelihood-based methods. However, once calibrated, the algorithm can fit 15-30 molecules per second on a 3 GHz Intel Core 2 Duo workstation, thereby producing a 3D super-resolution reconstruction of 100,000 molecules over a 20×20×2 μm field of view (processing 128×128 pixels × 20000 frames) in [75 min](http://www.nature.com/protocolexchange/protocols/2622#time_taken).

### Introduction
 
Single-molecule (SM)-based super-resolution microscopy can resolve neighboring objects separated closer than the diffraction limit of λ/(2 NA) ≈ 250 nm (1, 2), where λ is the wavelength of light and NA is the numerical aperture of the microscope’s objective lens. This is achieved through the fusion of several key ideas (3): 1) dense labeling of the structures of interest with SMs, 2) active control of the emission states of these SMs such that only a sparse, non-overlapping subset is emitting at any given point in time, and 3) precise measurement of the position of each emitting SM. Dense labeling is required because these techniques sample the underlying structure and reconstruct it in a pointillist manner; in order to avoid aliasing, the labeling density must satisfy the Nyquist-Shannon criterion (4, 5). Active control of SM emission states can be accomplished via a variety of mechanisms, as exemplified by the terms (fluorescence) Photoactivation Localization Microscopy [(f)PALM] (6, 7), STochastic Optical Reconstruction Microscopy (STORM) (8), fluorescent protein photocontrol/blinking (9) and Points Accumulation for Imaging in Nanoscale Topography (PAINT) (10), but we commonly refer to the family of SM super-resolution methods by the mechanism-independent moniker Single-Molecule Active Control Microscopy (SMACM) (2, 11).

Here, we focus on the third aforementioned key idea, namely the precise measurement of SM positions via widefield imaging using the double-helix (DH) microscope (12). The DH point spread function (DH-PSF) is an engineered rotating response added to a conventional widefield epifluorescence microscope that enables precise measurement of the three-dimensional (3D) position and orientation of SMs (13, 14). Three-dimensional position determination works as follows. A SM near focus in a DH microscope appears as two spots on a detector. The positions of these two spots can be measured via fitting of the image to the sum of two Gaussian functions in a least-squares manner (15), hereafter referred to as a double-Gaussian fit. The midpoint between these two spots yields the lateral (x, y) position of the SM, while the angle of the line connecting the two spots (relative to a fixed direction, such as the horizontal x axis) yields the axial (z) position of the SM from a separate calibration measurement. Thus, the 3D position of multiple, non-overlapping SMs can be measured with a single snapshot; no mechanical or optical scanning is required. If a molecule moves in the axial direction, the midpoint of the two spots remains constant while the two spots revolve around one another, hence tracing out a double-helix in 3D space. Axial position estimation is usually calibrated by scanning a bright fluorescent particle, such as a fluorescent bead, in known axial steps and measuring the rotation of the DH-PSF (see Figure 3).

The computational challenge in 3D SMACM with the DH microscope is to recognize the spot pairs that correspond to images of SMs accurately and precisely measure their locations to extract 3D position. Template matching is a process in which we try to find small, known patches t(x, y) in an image s(x, y) by selecting locations (p, q) that minimize the mean squared error quantity

![Eq 1](http://i.imgur.com/8BbK2Of.gif "Eq 1")

where FT{} denotes the two-dimensional Fourier transform, FT-1{} denotes the two-dimensional inverse Fourier transform, and the prefilter h is a Gaussian lowpass filter with a width σ ≈ 1.5 pixels. Phase correlation emphasizes high frequency components in both the input image and the template (16), and these spatial frequencies typically have a lower signal to noise ratio than lower frequencies in the image. The Gaussian lowpass filter h emphasizes lower spatial frequencies present in both the input image and the template and thus leads to more reliable estimates.

Local maxima (i.e. peaks) above a user-defined threshold (calibrated in Procedure, step 2) are detected in each image rk, and the locations and magnitudes of each peak are stored in memory. The width of the Gaussian prefilter and the threshold values are usually chosen empirically to minimize false positive matches (matches to images that are not the DH-PSF) and false negative matches (failure to match to an image of the DH-PSF) simultaneously. In general, the optimal value of these parameters can be derived if the additive noise in the original image is modeled with a simple noise distribution. However, typical live-cell SM data is complex, and modeling the noise distribution in the test data is beyond the scope of this work.

Since the peak threshold is usually set low enough to detect weak single molecules, extraneous matches occur for stronger fluorescent signals. If the candidate DH-PSF is bright enough, the template will match even if only one of the spots in the template overlaps with only one of the spots in the input. These “ghost” matches are filtered out by enforcing a minimum distance between all of the template matches in a given frame of the input data. This minimum distance is chosen to be the typical diameter of the DH-PSF image (~1.2 μm). If multiple template matches occur too close to one another, the strongest match is chosen, and the others are discarded. The template matches that satisfy this filter are validated template matches. These validated matches are then fit to the aforementioned double-Gaussian model of the DH-PSF (see Procedure, step 4).

The MATLAB optimization function lsqnonlin is used to fit each candidate SM image identified by the template-matching algorithm by minimizing the mean-square error between an eight-parameter double-Gaussian model of the DH-PSF U(x, y) and the raw SM image. This model is given by

![Eq 2](http://i.imgur.com/9XrHQv7.gif "Eq 2")

where the four parameters used for each circularly-symmetric Gaussian function are amplitude (A1 and A2), x-center location (μx1 and μx2), y-center location (μy1 and μy2), and width (σ1 and σ2). While this double-Gaussian model is only an approximation of the true shape of the DH-PSF (15), it represents a good compromise between fitting accuracy and computational complexity. These fits are then filtered to ensure that they match the data sufficiently accurately and produce reasonable reconstructions of the DH-PSF. This filtering step rejects a variety of unphysical situations as described below in Troubleshooting (step 4). The precision of these localizations can be estimated from measuring the number of photons contained within each image of the DH-PSF (13).

Sample drift can be corrected via imaging bright fluorescent beads that serve as fiduciary markers. These markers are tracked separately from the SMs throughout a fluorescence movie (see Procedure, step 3) to measure the drift. Their movements can then be subtracted from the three-dimensional localizations of SMs in order to remove thermal and mechanical motion artifacts from the super-resolution images.

Below, we demonstrate the [open-source Easy-DHPSF software](https://code.google.com/p/easy-dhpsf/) by analyzing a [SM dataset of Alexa647-immunolabeled microtubules in BSC-1 cells](https://code.google.com/p/easy-dhpsf/#Sample_single-molecule_DH-PSF_data) (17). We have included an easy-to-use [graphical user interface](http://www.nature.com/protocolexchange/protocols/2622#f0) to streamline the processing steps and minimize user error. This code has been [validated in several different imaging conditions](http://www.nature.com/protocolexchange/protocols/2622#related-articles) (17-19) and is provided free-of-charge under the New BSD License at https://code.google.com/p/easy-dhpsf/. However, this code is [provided ‘as is’ without guarantees of any kind](http://opensource.org/licenses/BSD-3-Clause).


### Equipment

**Software**

The Easy-DHPSF package consists of MATLAB functions and subroutines and requires the image processing, optimization, statistics, and wavelet toolboxes. Easy-DHPSF has been tested on workstations running MATLAB 2011b and 2012b. There are no known issues with earlier versions of MATLAB, but some built-in functions (e.g. fft2, imread) may have reduced performance in earlier versions of MATLAB.

**Imaging**

To create the double-helical PSF, it is necessary to use a spatial light modulator (12, 13) or a custom transmissive phase mask (19, 20). The Easy-DHPSF code anticipates that raw SM images are detected with an electron-multiplying (EM) CCD. If no EM gain is used, the ‘EM Gain’ fields in the program may simply be set to 1. It is necessary to measure the conversion gain (camera analog-to-digital converter counts per incident photon) of the EMCCD (21) and size of the camera pixels in object space before running this software. The combination of the EM gain and the conversion gain is used to compute the total number of detected photons.

Images should be acquired as .tif stacks readable by MATLAB’s imread function. Fluorescent beads or other bright, stationary, long-lasting emitters must be available in order to calibrate the z response of the DH-PSF; these may also be used to correct for drift if desired. A SM imaging dataset may be split over multiple sequential .tif files.

**Data file considerations**

Easy-DHPSF requires several .tif stacks (readable by MATLAB’s imread function) to run.

The first is the raw data file(s) with SMs to be localized. If desired, the field of view may contain one or more fiducials to correct for stage drift, ideally located near the edge of the field of view in order to not overlap with the region containing single molecules to be localized.

The second image stack is a calibration scan. This should comprise a series of frames of a fiducial at different known z positions generated using an axial nanopositioner in the microscope. The sequence of steps used to generate the scan is fed to the program in the form of a ‘sequence log file’, which should be an ASCII text file with a .dat extension and can be generated using the included ‘writeCalDatFile.m’ script. As shown in the example data, this scan should consist of approximately 50 steps sampling the entire z range of the DH-PSF, with ~10 acquisitions at each step to improve the precision of the measurement and provide localization precision information at each position.

The last image stacks are ‘dark offsets’ for the calibration and the data files. The frames in the dark offset .tif stacks are averaged and subtracted before data processing to account for the offset and dark counts from the EMCCD. Acquire a series of a few dozen frames with the camera shutters closed, and use the same acquisition parameters (e.g. integration time, EM gain, shift speed) as those used for the calibration and data acquisition. If these parameters are the same between the calibration and raw data images, only one dark offset file is necessary.

Easy-DHPSF works by processing these images with a series of modules that are run independently. Each of these saves its output to a separate directory, along with additional diagnostic figures. These modules and outputs are managed with a central GUI, which stores project data (e.g. detector information, file locations) as a MATLAB structure. Projects can be saved and loaded between sessions. Each module creates subfolders within the same folder as each raw .tif file that is processed, labeled with the name of the .tif file and the current time. For example, if it is noon on March 1, 2013, your calibration image is named ‘Calibration 1.tif’, and it is contained in C:\YourDataFolder, then the figures and data from the calibration will be saved in ‘C:\YourDataFolder\Calibration 1\calibration 20130301 1200.’ Ensure that you have read and write access to these folders.

### Procedure

**Step 0. Image acquisition and program setup**

Acquire the raw data, calibration, and dark offset images described in Equipment, and generate a .dat sequence log file.

Add the folder containing all Easy-DHPSF .m files to the MATLAB path using File>Set Path, and set the current folder to the folder containing your data. To initialize the GUI, run easy_dhpsf.m. This window (see Figure 1) controls all the modules needed to process the data.

Before running any of the modules, set the conversion gain and pixel size of the detector (in object/sample space) using the fields at the top of the main window. These parameters should be the same for all .tif files within a project. Setting the pixel size incorrectly may cause artifacts in the recognition and localization algorithms. Changing these parameters after processing the data does not retroactively change the calculated results; use caution!

![Fig 1](http://i.imgur.com/c3NOMha.png "Fig 1")

**Figure 1. The Easy-DHPSF command window: (left) before running any modules and (right) after localizing single molecules**.

**Step 1. DH-PSF Calibration**

The first module calibrates the relationship between the angle of the two lobes of the double helix and the axial position of single fiducial emitters, usually fluorescent beads. It also generates template images of the DH-PSF which are used for template matching later in the program.

- 1.1: Click the ‘run’ button in the ‘Calibrate DHPSF’ panel.
- 1.2: Set the EM gain used during the acquisition of the calibration images, and open the calibration images file, the dark offset image file matching the calibration, and the sequence log file for the calibration. The program will prompt for these data via dialog boxes.
- 1.3: Select a region of interest (ROI) containing at least one bright fluorescent object by adjusting the box and double-clicking inside it (see Figure 2). Even though the calibration bead may be isolated (bright object toward lower-left of excitation spot in Figure 2), it is helpful to select an entire widefield excitation region so that the background fluorescence can be estimated accurately.
- 1.4: Select fluorescent objects by clicking in the center of each DH-PSF. Each frame of raw data, as well as the accompanying double Gaussian fit, is displayed as it is processed (see Figure 3). Note that by mistake you might select a bright single molecule instead of a bead, but you will be able to detect and ignore this in subsequent steps.
- 1.5: After processing, an output of calibration information will be generated (see Figure 4). Inspect this for each bead to ensure that the angle vs. z position curve is roughly linear, and that the variation of angle measurements in each step is not unreasonably high (error bars in upper-left plot). Additional information from the calibration, such as the observed shift in x and y position with axial position, may be useful when aligning phase masks.
- 1.6: Select the most reliable calibration bead in the Easy-DHPSF main window using the drop-down menu. The calibration information and templates generated from this bead will be used in all later processing steps. 

![Fig 2](http://i.imgur.com/bV830Tj.png "Fig 2")

**Figure 2. Choose region of interest for calibrating the DH-PSF**.

![Fig 3](http://i.imgur.com/Fgqraq8.png "Fig 3")

**Figure 3. A single raw movie frame being analyzed by the DH-PSF calibration module**. Left: raw data after dark offset subtraction. Right: the reconstructed image of two beads in the raw data that were selected for fitting.

![Fig 4](http://i.imgur.com/pjWvifQ.png "Fig 4")

**Figure 4**. DH-PSF calibration statistics. For a good calibration standard, the calibration curve (upper-left) should be roughly linear, the xy position (upper-right) should deviate less than ±50 nm, and the localization precision (lower-left, lower-middle) in x, y, and z should be ≤15 nm over the course of the scan.

**Step 2. Single-Molecule Detection Calibration**

In this module, the templates generated from the calibration subroutine are used to generate a large array of matches to the raw image data. This step serves to assess what value of phase correlation is typical for a good match to a single molecule. After this step, the user will define a threshold for each template such that only DH-PSF images of reasonable quality will be analyzed with the double-Gaussian fitting algorithm.

- 2.1: Click the ‘run’ button in the ‘Calibrate SM identification’ panel.
- 2.2: Set the EM gain and select the templates to be used. By default, templates are chosen by the program such that the angle of the line connecting the two lobes is given by {-60°, -30°, 0°, 30°, 60°, 90°}, and this selection should generally be appropriate (see Figure 5).
- 2.3: Select a region containing the single molecules of interest. This exact ROI will be used for the single-molecule fitting module and ideally should not contain fluorescent objects that are much (≥10×) brighter than the SMs you wish to analyze. Template matches are indicated as circles drawn over the raw data, and stronger matches are drawn as larger circles (see Figure 6).
- 2.4: After the template-matching module is complete, open the ’threshold [DateAndTime]’ folder. This contains a selection of .png files that represent potential template matches. The filenames describe which template matched the data and the value of the phase correlation for the match (e.g. ‘template 3 threshold 272.png’).
- 2.5: Select appropriate thresholds for each template, and enter these into the ‘threshold’ field in the GUI (see Figure 1). These thresholds should be chosen such that the two lobes of the DH-PSF are faintly visible and that there is a low rate of false matches for higher correlation values, as in Figure 7. You should not be concerned if there are a few incorrect matches, say to one lobe of the DH-PSF from a very bright molecule, because these will be rejected by the subsequent double-Gaussian fit module. 

![Fig 5](http://i.imgur.com/4cNCk6l.png "Fig 5")

**Figure 5. Typical DH-PSF templates**. Each template is chosen from the aforementioned calibration scan, where the rotation of the double helix between templates is evenly spaced at ~30°.

![Fig 6](http://i.imgur.com/Se0MZ8S.png "Fig 6")

**Figure 6. A single raw movie frame being analyzed by the SM-detection calibration module**. Left: phase correlation of the raw image with the DH templates. Right: raw data with template matches indicated by colored circles. Larger circles indicate stronger matches, and will be saved with a higher threshold value.

![Fig 7](http://i.imgur.com/gbXgVw6.png "Fig 7")

**Figure 7. Potential DH-PSF SMs identified by template matching**. The image with an appropriate threshold is highlighted. Note that while thresholds 88, 92, 94, 102, 106, and 107 represent poor matches to the DH-PSF, the other matches for correlation values ≥84 are satisfactory.

**Step 3. Fiducial Tracking (Optional)**

This module tracks the movement of one or more fluorescent beads or other stationary markers. Stage drift during a SM experiment can be thus be removed by subtracting the movement of the fiduciary marker from the SM localizations.

- 3.1: Click the ‘run’ button in the ‘Track fiduciaries’ panel.
- 3.2: Select a region of interest containing at least one bright fiducial by adjusting the box and double-clicking inside it.
- .3: Select one or more fiducials by clicking in the center of each DH-PSF.
- 3.4: After fitting the fiducials, select whether or not to apply the averaged correction to the output data by clicking the ‘use fiducials’ checkbox. This may be changed at any time before outputting the data.

**Step 4. Single-Molecule Localization**

Using the results of the previous processing steps, this module identifies single molecules that score above the template matching threshold identified in step 2.5, then fits them using nonlinear least squares minimization to a double-Gaussian function.

- 4.1: Click the ‘run’ button in the ‘Localize DHPSF SMs’ panel.
- 4.2: Choose the frames from the raw data to process: by default, the entire image stack is processed.
- 4.3: Whenever starting a new experiment, monitor this process as it proceeds to ensure that the threshold used includes the single molecules visible in the raw data, and to ensure that the fits are performed successfully. Only matches that score above the thresholds are fit, and good fits are plotted in the frame-by-frame reconstruction (see Figure 8). Here the definition of a good double-Gaussian fit is a fit that meets several straightforward criteria, hard-coded into the module. The tests performed to define a good fit are listed below under Troubleshooting (step 4).

![Fig 8](http://i.imgur.com/jabXGHC.png "Fig 8")

**Figure 8. A single raw movie frame being analyzed by the SM-fitting module**. Left: phase correlation of the raw image with the DH templates. Center: raw data with template matches indicated by circles. Right: the image reconstructed from good double-Gaussian fits.

**Step 5. Viewing/Exporting the Processed Data**

Easy-DHPSF gives three options for output of the processed data after all modules have been run, which are discussed in greater detail in the Anticipated Results.

- 5.1: Click ‘Export to csv’ to generate a list of localizations in .csv format for direct manipulation. The columns of this file are labeled. Note that this listed z position does NOT correct for index of refraction mismatch between the sample and the objective! The user should supply that correction to the data before producing final reconstructions.
- 5.2: Click ‘3D scatterplot’ to generate a scatterplot of localizations at a particular ROI. Set the index of refraction for the immersion lens and the sample to apply a simple correction for index mismatch; by default, oil immersion and an aqueous sample are assumed.
- 5.3: Click ‘2D histogram’ to generate a histogram of localizations with median z position for each bin coded by color. As with the scatterplot, this module will correct for index mismatch.

### Timing 

*Benchmarked on a 3 GHz Intel Core 2 Duo workstation running 64-bit Windows 7 Professional with 8 GB of RAM*

- Step 1. DH-PSF Calibration
  - 2 min (2 beads, 256×256 pixels × 500 frames [50 z-steps, 10 acquisitions per step])
- Step 2. Single-Molecule Detection Calibration
  - 36 min (128×128 pixels × 20,000 frames)
- Step 3. Fiducial Tracking
  - 13 min (1 bead, 64×64 pixels × 9,000 frames)
- Step 4. Single-Molecule Localization
  - 75 min (100,000 molecules in 128×128 pixels × 20,000 frames)

### Troubleshooting

**General troubleshooting**

If the phase mask used to generate the DH-PSF is grossly misaligned, then the fitting may not be feasible. Problems with the alignment should be noticeable in the output from the calibration module.

If the processed data are moved after saving a project file, the GUI will not be able to locate it. Please avoid moving processed data once it is created. Otherwise, manually edit the file paths saved in the project .MAT file to restore operation.

If the user does not have read and write access to the directory where the raw datafiles reside, file access errors will occur. Move the raw data to a folder for which you have read and write access.

**Step 2. Single-Molecule Detection Calibration**

If the ROI chosen for the template calibration is not circularly periodic (i.e., the edges of the ROI do not have approximately the same signal level from left to right, top to bottom), it is possible that erroneous template matches will appear. If this occurs, change the ROI such that its edges are approximately continuous.

**Step 4. Single-Molecule Localization**

If the DH-PSF fitting algorithm has a high rate of false negatives (not fitting clearly visible images of the DH-PSF present in the raw data), there are several possible causes:

First, if the threshold has been set too high, then good PSFs identified by the template-matching algorithm may be thrown out and never fit. If this occurs, the central panel of the fitting module (see Figure 8) will not show circles around DH-PSFs. To fix this, lower the thresholds.

If the localization fails during the double-Gaussian fitting step, it is likely that the DH-PSF has failed one or more filtering criteria. These are flagged, and can be viewed by hitting the ‘debug’ button after the fitting step is complete. Doing so generates a .csv file containing all data for good and bad fits, as well as a diagnostic showing how many fits failed, and why. This will indicate possible reasons for the errors in fitting. For example, if the pixel size is set improperly, many bounds used in assuring that the fits are valid will be distorted. This might lead to many localizations being flagged as having values of σ for the double-Gaussian fit that are out of bounds. The possible flags are:

- -1007: Error in the linear least squares fit too high
- -1006: Ratio of the amplitudes of the two Gaussian lobes too high
- -1005: Distance between DH-PSF lobes out of bounds
- -1004: Ratio of σ values out of bounds
- -1003: Absolute σ values out of bounds
- -1002: One or more peaks in the fit too far away from template match
- -1001: Amplitude of one or more of the peaks <0
- -1000: Initial guess of location outside of selected ROI
 
### Anticipated Results
 
Upon successful fitting of the SM data, the user may use the ‘Output DHPSF SM localizations’ panel to export the SM 3D positions for further post-processing or visualization of the data. Clicking on the ‘Export to csv’ button will prompt the user for a filename and then subsequently save the valid 3D localizations to a comma-separated text file. Note that these localizations have not been corrected for focal shift resulting from imaging into mismatched media.

The program can also plot the 3D localizations as solid circles in a 3D scatterplot; note that this type of plot is quick to render but conveys localization density poorly when localizations are dense, and the validity of the use of a crude scatterplot for publication is the responsibility of the user. Clicking the ‘3D scatterplot’ button will first prompt the user to verify settings used to plot the SM data. For example, the immersion lens index of refraction and sample media index of refraction can be entered to correct the SMs’ z positions for focal shift due to imaging into mismatched media (note that this correction of ![Eq 3](http://i.imgur.com/pYSmvyu.gif "Eq 3") is an approximation). Furthermore, the user can specify if they want the scatter points to be color-coded in time (useful for discerning sample drift, for example). The user then selects a region of interest (ROI) to plot in 3D. Once displayed, the 3D scatterplot can be manipulated using the zoom, pan, and rotate tools in the figure toolbar.

![Fig 9](http://i.imgur.com/4ITN9oe.png "Fig 9")

**Figure 9. Select region of interest for 3D scatterplot**.

![Fig 10](http://i.imgur.com/Gyctg5Q.png "Fig 10")

**Figure 10. Sample 3D scatterplot of a section of microtubules**.

For continuous extended structures, the ‘2D histogram’ function also available. Again, the validity of the use of any particular reconstruction for publication is the responsibility of the user. This function subdivides the ROI into a series of square 2D xy bins of a user-specified size. It then calculates the number and median z position of the localizations within each bin. Similar to the scatterplot function, clicking on the ‘2D histogram’ button also prompts the user for histogram plotting parameters and an ROI. Then two plots appear, one with the median z position color-coded for each bin (no localization number/density information) and the other with both z position color-coded for each bin and number of localizations proportional to the brightness of each bin.

![Fig 11](http://i.imgur.com/lNP7wmQ.png "Fig 11")

**Figure 11. Sample histograms of labeled microtubules conveying 3D localization information**. Left: Localization bins color-coded by median z position. Right: Localization bins color-coded by both number of localizations and median z position to convey density information.

### References

1. Hell, S. W. Microscopy and its focal switch. *Nat. Methods* 6, 24-32 (2009).
- Moerner, W. E. Microscopy beyond the diffraction limit using actively controlled single molecules. *J. Microsc*. 246, 213-220 (2012).
- Thompson, M. A., Lew, M. D. & Moerner, W. E. Extending Microscopic Resolution with Single-Molecule Imaging and Active Control. Annu. Rev. Biophys. 41, 321-342 (2012).
- Nyquist, H. Certain Topics in Telegraph Transmission Theory. *Trans. AIEE* 47, 617-644 (1928).
- Shannon, C. E. Communication in the Presence of Noise. Proc. IRE 37, 10-21 (1949).
- Hess, S. T., Girirajan, T. P. K. & Mason, M. D. Ultra-high resolution imaging by fluorescence photoactivation localization microscopy. *Biophys. J*. 91, 4258-4272 (2006).
- Betzig, E. et al. Imaging Intracellular Fluorescent Proteins at Nanometer Resolution. *Science* 313, 1642-1645 (2006).
- Rust, M. J., Bates, M. & Zhuang, X. Sub-diffraction-limit imaging by stochastic optical reconstruction microscopy (STORM). *Nat. Methods* 3, 793-796 (2006).
- Biteen, J. S. et al. Super-resolution imaging in live Caulobacter crescentus cells using photoswitchable EYFP. *Nat. Methods* 5, 947-949 (2008).
- Sharonov, A. & Hochstrasser, R. M. Wide-field subdiffraction imaging by accumulated binding of diffusing probes. *Proc. Natl. Acad. Sci*. USA 103, 18911-18916 (2006).
- Biteen, J. S., Thompson, M. A., Tselentis, N. K., Shapiro, L. & Moerner, W. E. Superresolution Imaging in Live Caulobacter crescentus Cells Using Photoswitchable Enhanced Yellow Fluorescent Protein. *Proc. SPIE* 7185, 71850I (2009).
- Pavani, S. R. P. et al. Three-dimensional, single-molecule fluorescence imaging beyond the diffraction limit by using a double-helix point spread function. *Proc. Natl. Acad. Sci. USA* 106, 2995-2999 (2009).
- Thompson, M. A., Lew, M. D., Badieirostami, M. & Moerner, W. E. Localizing and Tracking Single Nanoscale Emitters in Three Dimensions with High Spatiotemporal Resolution Using a Double-Helix Point Spread Function. *Nano Lett*. 10, 211-218 (2010).
- Backlund, M. P. et al. Simultaneous, accurate measurement of the 3D position and orientation of single molecules. *Proc. Natl. Acad. Sci. USA* 109, 19087-19092 (2012).
- Lew, M. D., Thompson, M. A., Badieirostami, M. & Moerner, W. E. In vivo three-dimensional superresolution fluorescence tracking using a double-helix point spread function. *Proc. SPIE* 7571, 75710Z (2010).
- Brunelli, *R. in Template matching techniques in computer vision: theory and practice* (Wiley, Chichester, West Sussex, U.K., 2009).
- Lee, H. D., Sahl, S. J., Lew, M. D. & Moerner, W. E. The double-helix microscope super-resolves extended biological structures by localizing single blinking molecules in three dimensions with nanoscale precision. *App. Phys. Lett*. 100, 153701 (2012).
- Lew, M. D. et al. Three-dimensional superresolution colocalization of intracellular protein superstructures and the cell surface in live Caulobacter crescentus. *Proc. Natl. Acad. Sci. USA* 108, E1102-E1110 (2011).
- Gahlmann, A. et al. Quantitative multicolor subdiffraction imaging of bacterial protein ultrastructures in 3D. *Nano Lett*. (2013).
- Grover, G., Quirin, S., Fiedler, C. & Piestun, R. Photon efficient double-helix PSF microscopy with application to 3D photo-activation localization imaging. *Biomed Opt Expr* 2, 3010-3020 (2011).
- http://www.mirametrics.com/tech_note_ccdgain.htm.

### Acknowledgements
 
We thank and recognize everyone who has contributed to writing, debugging, and generally improving this code over the years: Andreas Gahlmann, Mikael P. Backlund, Steven F. Lee, Steffen J. Sahl, Alex Chang, and Scott S. Hsieh. We also acknowledge R. Piestun for providing the DH-PSF phase mask design. This work was supported by a National Science Foundation Graduate Research Fellowship and 3Com Corporation Stanford Graduate Fellowship (to M.D.L.) and a National Institute of General Medical Sciences Grant R01GM085437 (to W.E.M.).

### Associated Publications
 
1. **PNAS Plus: Three-dimensional superresolution colocalization of intracellular protein superstructures and the cell surface in live Caulobacter crescentus**. M. D. Lew, S. F. Lee, J. L. Ptacin, M. K. Lee, R. J. Twieg, L. Shapiro, and W. E. Moerner. *Proceedings of the National Academy of Sciences* 108 (46) E1102 - E1110 15/11/2011 [doi:10.1073/pnas.1114444108](http://dx.doi.org/10.1073/pnas.1114444108)
- **The double-helix microscope super-resolves extended biological structures by localizing single blinking molecules in three dimensions with nanoscale precision**. Hsiao-lu D. Lee, Steffen J. Sahl, Matthew D. Lew, and W. E. Moerner. *Applied Physics Letters* 100 (15) [doi:10.1063/1.3700446](http://dx.doi.org/10.1063/1.3700446)
- **Quantitative multicolor subdiffraction imaging of bacterial protein ultrastructures in 3D**. Andreas Gahlmann, Jerod Louis Ptacin, Ginni Grover, Sean Quirin, Alexander R.S. von Diezmann, Marissa K. Lee, Mikael Paul Backlund, Lucy Shapiro, Rafael Piestun, and W.E. Moerner. *Nano Letters* 15/02/2013 [doi:10.1021/nl304071h](http://dx.doi.org/10.1021/nl304071h)

### Author information
 
**Matthew D. Lew, Alexander R. S. von Diezmann & W. E. Moerner**, Moerner Laboratory, Department of Chemistry, Stanford University

 Correspondence to: W. E. Moerner (wmoerner@stanford.edu)






*Source: [Protocol Exchange](http://www.nature.com/protocolexchange/protocols/2622) (2013) doi:10.1038/protex.2013.026. Originally published online 25 February 2013*. 