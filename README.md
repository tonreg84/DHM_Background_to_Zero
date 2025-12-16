DHM background to zero
Autor: tonreg, team UMI, CNP-CHUV Lausanne

Open a "bnr" file from a a LynceeTec DHM recording.
Control background search and fitting parameters.
Flatten all images of a bnr sequence.
    
A) Make a "noisy mask" to exclude pixels from further image treatment
    - Excludes pixels with very low phase values (supposedly wrong due to phase jumps).
    - Excludes noisy pixels by noise detection by phase treshold. Pixels which vary between two frames more that the treshold will be added to the exclusion mask.
    - Mask from excluded pixels will be enhanced by dilation and hole filling

low_treshold_para = 115 in degrree, parameter to exclude pixel with very low values (supposedly wrong due to phase jumps)
noise_treshold = 15 in degree, pixels which vary between two frames more that this value, will be added to the exclusion mask.
extend_para = 2 # Number of pixels by which the final mask will be extendet

B) Make a "global mask" with iterative background estimation
Detect the background of a sequence of DHM images. Uses Interative_BG_detection* to estimate the backgrounds of a number of frames.
Combines the masks to obtain the final background, i.e., makes a combination of all (logical AND).
E.g.: 1=middle frame, 2=first and last frame, 3=first, middle, and last frame,...

*Interative_BG_detection: Estimates image background by iterative background detection and polynomial fit to the background.

polynom_order = 3 # Polynomial order for the background fit.
noise_delta = 10 in degree, Some tolerance to noise in the background. E.g. diffraction fringes from dust particel in the optical path.
convergence_treshold = 80%, Parameter to check convergence; Stop iteration, when given percentage of image surface is foreground.
mask_var_trsh = 0.01 # Mask variantion treshold; Parameter for adaptive tresholding; stop when the change of mask becomes negligible.
bg_frame_num = 10 # Number of frames on which the background extracted. Final background is combination of all (logical AND). E.g.: 1=middle frame, 2=first and last frame, 3=first, middle, and last frame,...

C) Flattening all images of a bnr sequence.
On every frame of the sequencs: Applies the background mask, fits the background, substracts the resulting fit from the image.

polynom_order = 3 # Polynomial order for the background fit.
noise_delta = 10 in degree, Some tolerance to noise in the background. E.g. diffraction fringes from dust particel in the optical path.
