DHM background to zero<br/>
Autor: tonreg, team UMI, CNP-CHUV Lausanne<br/>
Version: 20260108

Programme to flatten phase images obtained with Digital Holographic Microscopes. Find image background. Polynomial fit to background. Subtract fitted background from image.

Open an image sequence from a LynceeTec DHM recording. Control background search and fitting parameters. Flatten all images of a bnr sequence.<br/>
Accepted file format: "LynceeTec format" ("bin", "bnr"). The frames of the sequence (bin or bnr) are 2D arrays of float.<br/>
* "bin" - Select any bin file in the folder. Bin files are single frames of the sequence.<br/>
* "bnr" - Binary file containing the whole sequence.<br/>

A) Make a "noisy mask" to exclude pixels from further image treatment:<br/>

* Excludes pixels with very low phase values (supposedly wrong due to phase jumps).
* Excludes noisy pixels by noise detection by phase treshold. Pixels which vary between two frames more that the treshold will be added to the exclusion mask.
* Mask from excluded pixels will be enhanced by dilation and hole filling

low\_treshold\_para = 115; Unit: Degree; Parameter to exclude pixel with very low values (supposedly wrong due to phase jumps).<br/>
noise\_treshold = 15; Unit: Degree; Pixels which vary between two frames more that this value, will be added to the exclusion mask.<br/>
extend\_para = 2; Number of pixels by which the final mask will be extendet.

B) Make a "global mask" with iterative background estimation:<br/>
Detects the background of a sequence of DHM images. Uses Interative\_BG\_detection\* to estimate the backgrounds of a number of frames.<br/>
Combines the masks to obtain the final background, i.e., makes a combination of all (logical AND).

\*Interative\_BG\_detection: Estimates image background by iterative background detection and polynomial fit to the background.

polynom\_order = 3; Polynomial order of the background fit.<br/>
noise\_delta = 10; Unit: Degree; Some tolerance to noise in the background. E.g. diffraction fringes from dust particel in the optical path.<br/>
convergence\_treshold = 80%; Parameter to check convergence; Stop iteration, when given percentage of image surface is foreground.<br/>
mask\_var\_trsh = 0.01; Mask variantion treshold; Parameter for adaptive tresholding; Stop when the change of mask becomes negligible.<br/>
bg\_frame\_num = 10; Number of frames on which the background is extracted. Final background is combination of all (logical AND). E.g.: 1 = middle frame, 2 = first and last frame, 3 = first, middle, and last frame,...

C) Flattening all images of a bnr sequence:<br/>
On every frame of the sequencs: Applies the background mask, fits the background, substracts the resulting fit from the image.<br/>

polynom\_order = 3; Polynomial order of the background fit.

