# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 10:04:20 2025

@author: ge1582

"phase image iterative background-to-zero"

Tools to post-process cell culture images recorded with a Digital Holographic Microscope.            
"""
import os
import numpy as np
from scipy.linalg import lstsq
import matplotlib.pyplot as plt
import binkoala
from scipy.ndimage import binary_fill_holes, binary_dilation
from tkinter import messagebox
from scipy.ndimage import map_coordinates


def diagonal_profile(array):
    N, M = array.shape
    # Number of points along the diagonal
    num_points = max(N, M)  # take enough points to cover diagonal proportionally
    
    # Generate coordinates along the diagonal
    rows = np.linspace(0, N-1, num_points)
    cols = np.linspace(0, M-1, num_points)
    
    # Sample values along the diagonal
    diag = map_coordinates(array, [rows, cols])
    return diag


def Make_noisy_mask(path,low_treshold_para,noise_treshold,extend_para):
    '''
    Make a mask to exclude pixels from further image treatment
    - Excludes pixels with very low phase values (supposedly wrong due to phase jumps).
    - Excludes noisy pixels by noise detection by phase treshold.
    - Mask from excluded pixels will be enhanced by dilation and hole filling
    
    Args:
        - path (string): File path to the phase image sequence in "LynceeTec format" ("bin", "bnr")
                        - If "bin", select any bin file in the folder. Bin files are single frames of the sequence.
                        - "bnr" - binary file containing the whole sequence.
                        The frames of the sequence (bin or bnr) are 2D arrays of float.
        - low_treshold_para (float): Value to be subtracted from average phase of image to set lower treshold.
                                    All pixels with values below lower treshold will be added to the mask
        - noise_treshold (float): Treshold of phase jump; Pixels which vary between two frames more that the treshold,
                                 will be added to the mask.
        - extend_para (int): Number of pixels by which the final mask will be extendet
    Returns:
        final_mask (numpy.ndarray of boolean): Exclusion mask (pixels to exclude = True)
    '''
    print("Make noisy mask")
    file_name, file_type = os.path.splitext(path) # get file suffix

    image_ref = None
    image_src = None

    # Depending on file suffix, read sequence lenght and file header
    go_on = True
    if file_type == ".bin":
        
        # get the list of bin files to process
        binfolder=os.path.dirname(path)
        bin_files = []
        bin_files = sorted([f for f in os.listdir(binfolder) if f.endswith(('.bin'))])
        nImages = len(bin_files) # sequence length
        
        # Get header info and 1st image from sequence
        (image_ref,in_file_header)=binkoala.read_mat_bin(binfolder +os.sep+ bin_files[0])
        hv=str(in_file_header['version'][0])
        end=str(in_file_header['endian'][0])
        hz=str(in_file_header['head_size'][0])
        W=str(in_file_header['width'][0])
        H=str(in_file_header['height'][0])
        pz=in_file_header['px_size'][0]
        hconv=str(in_file_header['hconv'][0])
        uc=str(in_file_header['unit_code'][0])
        
    elif file_type == ".bnr":
        # Get header info from bnr file
        fileID = open(path, 'rb')
        nImages = np.fromfile(fileID, dtype="i4", count=1)
        nImages = nImages[0]  # sequence length
        w = np.fromfile(fileID, dtype="i4", count=1)
        W=w[0]
        h = np.fromfile(fileID, dtype="i4", count=1)
        H=h[0]
        pz = np.fromfile(fileID, dtype="f4", count=1)
        pz=pz[0]
        wavelength = np.fromfile(fileID, dtype="f4", count=1)
        wavelength=wavelength[0]
        n_1 = np.fromfile(fileID, dtype="f4", count=1)
        n_2 = np.fromfile(fileID, dtype="f4", count=1)
        timestamps = [0] * nImages
        for k in range(0,nImages):
            timestamps[k] = np.fromfile(fileID, dtype="f4", count=1)
            
        # get 1st image from sequence
        image_ref = np.zeros((H,W))
        image_src = np.zeros((H,W))
        for k in range(H):
            image_ref[k,:] = np.fromfile(fileID, dtype="f4", count=W)
        
    else:
        print("ERROR: Wrong file type")
        go_on = False
        return None
        
    if go_on:
        # # Plotting (Optional)
        # plt.imshow(image_ref)
        # plt.title("phase image")
        # plt.show()

        # Initialize the mask
        mask = np.zeros_like(image_ref, dtype=bool)
    
        # For the first frame, exclude pixels with very low phase values (supposedly wrong due to phase jumps).
        phaseavg=np.mean(image_ref)
        low_teshold = phaseavg - low_treshold_para
        mask[image_ref < low_teshold] = True
    
        # Loop through frames
        for i in range(1,nImages):
            
            print(f"Processing frame {i} of {nImages}...")
    
            # get next image from sequence
            if file_type == ".bin":
                (image_src,in_file_header)=binkoala.read_mat_bin(binfolder +os.sep+ bin_files[i])
            elif file_type == ".bnr":
                for k in range(H):
                    image_src[k,:] = np.fromfile(fileID, dtype="f4", count=W)

            # Exclude pixels with very low phase values (supposedly wrong due to phase jumps).
            phaseavg=np.mean(image_src)
            low_teshold = phaseavg - low_treshold_para
            mask[image_src < low_teshold] = True
            
            # exlude noisy pixels
            mask[abs(image_src-image_ref)>=noise_treshold] = True
                
            image_ref = image_src.copy()
        
        if file_type == ".bnr":
            fileID.close
    
        # # Plotting (Optional)
        # plt.imshow(mask)
        # plt.title("noise mask")
        # plt.show()
    
        # Extend/dilate mask by N pixels
        extend_para = 2  # number of pixels to extend
        structure = np.ones((2*extend_para+1, 2*extend_para+1), dtype=bool)  # square structuring element
        dilated_mask = binary_dilation(mask, structure=structure)
        
        # # Plotting (Optional)
        # plt.imshow(dilated_mask)
        # plt.title("dilated mask")
        # plt.show()
    
        # Fill holes in mask again
        final_mask = binary_fill_holes(dilated_mask)
    
        # # Plotting (Optional)
        # plt.imshow(final_mask)
        # plt.title("filled 2 mask")
        # plt.show()
        
        print("Noisy mask created")
        return final_mask


def fit_2D(image, Poly_Ord):
    """
    Fit a 2D polynomial to the input data.

    Args:
        image (numpy.ndarray): input image array.
        Poly_Ord (int): Polynomial order of the surface fit
        pix (float): Pixel size.

    Returns:
        tuple: Tuple containing the 2D fit of shape image.shape and polynomial coefficients.
    """
    # Generate grid coordinates
    x = np.arange(0, image.shape[1]) # height 
    y = np.arange(0, image.shape[0]) # width
    X, Y = np.meshgrid(x, y)
    
    # Flatten arrays
    z = image.reshape(image.shape[0] * image.shape[1])
    xx = X.reshape(image.shape[0] * image.shape[1])
    yy = Y.reshape(image.shape[0] * image.shape[1])
    
    # Fit polynomial based on the order
    if Poly_Ord == 1:
        AMx = np.column_stack([xx, yy, np.ones(len(xx))])
    elif Poly_Ord == 2:
        AMx = np.column_stack([xx, yy, xx**2, yy**2, xx*yy, np.ones(len(xx))])
    elif Poly_Ord == 3:
        AMx = np.column_stack([xx, yy, xx**2, yy**2, xx*yy, xx**3, yy**3, xx**2*yy, xx*yy**2, np.ones(len(xx))])
    elif Poly_Ord == 4:
        AMx = np.column_stack([xx, yy, xx**2, yy**2, xx*yy, xx**3, yy**3, xx**2*yy, xx*yy**2, xx**4, yy**4,
                               xx**3*yy, xx*yy**3, xx**2*yy**2, np.ones(len(xx))])
    elif Poly_Ord == 5:
        AMx = np.column_stack([xx, yy, xx**2, yy**2, xx*yy, xx**3, yy**3, xx**2*yy, xx*yy**2, xx**4, yy**4,
                               xx**3*yy, xx*yy**3, xx**2*yy**2, xx**5, yy**5, xx**4*yy, xx*yy**4, xx**3*yy**2,
                               xx*yy**3, np.ones(len(xx))])
    elif Poly_Ord == 6:
        AMx = np.column_stack([xx, yy, xx**2, yy**2, xx*yy, xx**3, yy**3, xx**2*yy, xx*yy**2, xx**4, yy**4,
                               xx**3*yy, xx*yy**3, xx**2*yy**2, xx**5, yy**5, xx**4*yy, xx*yy**4, xx**3*yy**2,
                               xx*yy**3, xx**6, yy**6, xx**5*yy, xx*yy**5, xx**4*yy**2, xx**2*yy**4,
                               xx**3*yy**3, xx*yy**6, xx**6*yy, np.ones(len(xx))])
    else:
        raise ValueError("Invalid polynomial order. Supported orders are from 1 to 6.")
    
    # Solve least squares problem
    P = lstsq(AMx, z)[0]
    
    # Compute fitted 2D polynomial
    fit2Dim = np.dot(AMx, P).reshape(image.shape)

    return fit2Dim, P


def masked_fit_2D(image, Poly_Ord, mask):
    """
    Fit a 2D polynomial to the input data.

    Args:
        image (numpy.ndarray): input image array.
        Poly_Ord (int): Polynomial order of the surface fit
        mask (numpy.ndarray of boolean): Mask to exclude pixel during linear fit (pixels to exclude = True)
    Returns:
        tuple: Tuple containing the 2D fit of shape image.shape and polynomial coefficients.
    """
    # Generate grid coordinates
    x = np.arange(0, image.shape[1]) # height 
    y = np.arange(0, image.shape[0]) # width
    X, Y = np.meshgrid(x, y)
    
    # Flatten arrays
    z = image.reshape(image.shape[0] * image.shape[1])
    mmm = mask.reshape(mask.shape[0] * mask.shape[1])
    xx = X.reshape(image.shape[0] * image.shape[1])
    yy = Y.reshape(image.shape[0] * image.shape[1])
    
    xm = xx[~mmm]
    ym = yy[~mmm]
    zm  = z[~mmm]
    
    # Fit polynomial based on the order
    if Poly_Ord == 1:
        AMxM = np.column_stack([xm, ym, np.ones(len(xx))])
        AMx = np.column_stack([xx, yy, np.ones(len(xx))])
    elif Poly_Ord == 2:
        AMxM = np.column_stack([xm, ym, xm**2, ym**2, xm*ym, np.ones(len(xm))])
        AMx = np.column_stack([xx, yy, xx**2, yy**2, xx*yy, np.ones(len(xx))])
    elif Poly_Ord == 3:
        AMxM = np.column_stack([xm, ym, xm**2, ym**2, xm*ym, xm**3, ym**3, xm**2*ym, xm*ym**2, np.ones(len(xm))])
        AMx = np.column_stack([xx, yy, xx**2, yy**2, xx*yy, xx**3, yy**3, xx**2*yy, xx*yy**2, np.ones(len(xx))])
    elif Poly_Ord == 4:
        AMxM = np.column_stack([xm, ym, xm**2, ym**2, xm*ym, xm**3, ym**3, xm**2*ym, xm*ym**2, xm**4, ym**4,
                               xm**3*ym, xm*ym**3, xm**2*ym**2, np.ones(len(xm))])
        AMx = np.column_stack([xx, yy, xx**2, yy**2, xx*yy, xx**3, yy**3, xx**2*yy, xx*yy**2, xx**4, yy**4,
                               xx**3*yy, xx*yy**3, xx**2*yy**2, np.ones(len(xx))])
    elif Poly_Ord == 5:
        AMxM = np.column_stack([xm, ym, xm**2, ym**2, xm*ym, xm**3, ym**3, xm**2*ym, xm*ym**2, xm**4, ym**4,
                               xm**3*ym, xm*ym**3, xm**2*ym**2, xm**5, ym**5, xm**4*ym, xm*ym**4, xm**3*ym**2,
                               xm*ym**3, np.ones(len(xm))])
        AMx = np.column_stack([xx, yy, xx**2, yy**2, xx*yy, xx**3, yy**3, xx**2*yy, xx*yy**2, xx**4, yy**4,
                               xx**3*yy, xx*yy**3, xx**2*yy**2, xx**5, yy**5, xx**4*yy, xx*yy**4, xx**3*yy**2,
                               xx*yy**3, np.ones(len(xx))])
    elif Poly_Ord == 6:
        AMxM = np.column_stack([xm, ym, xm**2, ym**2, xm*ym, xm**3, ym**3, xm**2*ym, xm*ym**2, xm**4, ym**4,
                               xm**3*ym, xm*ym**3, xm**2*ym**2, xm**5, ym**5, xm**4*ym, xm*ym**4, xm**3*ym**2,
                               xm*ym**3, xm**6, ym**6, xm**5*ym, xm*ym**5, xm**4*ym**2, xm**2*ym**4,
                               xm**3*ym**3, xm*ym**6, xm**6*ym, np.ones(len(xm))])
        AMx = np.column_stack([xx, yy, xx**2, yy**2, xx*yy, xx**3, yy**3, xx**2*yy, xx*yy**2, xx**4, yy**4,
                               xx**3*yy, xx*yy**3, xx**2*yy**2, xx**5, yy**5, xx**4*yy, xx*yy**4, xx**3*yy**2,
                               xx*yy**3, xx**6, yy**6, xx**5*yy, xx*yy**5, xx**4*yy**2, xx**2*yy**4,
                               xx**3*yy**3, xx*yy**6, xx**6*yy, np.ones(len(xx))])
    else:
        raise ValueError("Invalid polynomial order. Supported orders are from 1 to 6.")
    
    # print(type(AMx))
    # print(AMx.shape)
    # print(type(AMxM))
    # print(AMxM.shape)
    
    # Solve least squares problem
    P = lstsq(AMxM, zm)[0]
    
    # Compute fitted 2D polynomial
    fit2Dim = np.dot(AMx, P).reshape(image.shape)

    return fit2Dim, P


def CheckForLess(list1, val):
    return(all(x < val for x in list1)) 


def Interative_BG_detection(phase, poly_order, pix, method, starting_mask, noise_delta, convergence_treshold, mask_var_trsh):
    '''
    starting from function baseLine2DSurFit from https://github.com/lrnp/DHM_Reconstruction/blob/main/Fct_PDHM_Rec.py
    (version 21.07.2025)
    added "starting mask", "evolving mask", more plotting, adaptive thresholding
    
    Estimate image background by iterative background detection and polynomial fit to the background.

    Args:
       phase (numpy.ndarray): DHM phase image
       Poly_Ord (int): Polynomial order for the background fit.
       pix (float): Pixel size.
       method (str): Method of background selection. #TODO
                       - "clamp" : Clamp test image (y0) to not exceed background (bk) -> y0 = np.minimum(y0, bk)
                                   This ensures that the evolving test image never exceeds the current baseline.
                       - "mask" ***under construction*** : Create binary mask for evolving foregound/background, to exlude foreground from surface fit.
       starting_mask (numpy.ndarray of boolean): Initial mask to exclude parts from linear fit (pixels to exclude = True),
                                                 If "None" or "False": start with empty mask
       noise_delta (float): Constant to be added to background fit to account for noise in the background
       convergence_treshold (float): Parameter to check convergence; stop when given percentage of image surface is foreground
       mask_var_trsh (float): # Parameter for adaptive tresholding; stop when the change of mask becomes negligible
    Returns:
       bk (numpy.ndarray of float): Estimated background, same shape as phase
       bg_mask (numpy.ndarray of boolean): Final mask of foregound/background (foreground = True)
    '''
    print("Estimate background...")
    # Get the shape of the phase image
    N, M = phase.shape
    y0 = phase.copy() # Initialize test image
    
    # Initialize mask (foreground = True)
    if isinstance(starting_mask, np.ndarray):
        bg_mask = starting_mask.copy()
    else:
        bg_mask = np.zeros_like(phase, dtype=bool)
    
    # # Fit a 2D polynomial to the raw phase image
    # bk, _ = fit_2D(y0, poly_order, pix)
    # /
    # Fit a 2D polynomial to the raw phase image, exluding masked pixels
    bk, _ = masked_fit_2D(y0, poly_order, bg_mask)#TODO
    
    # # Plotting (Optional)
    # xl = np.arange(N)
    # plt.figure()
    # plt.plot(xl, diagonal_profile(y0), "g", label="Raw")
    # plt.plot(xl, diagonal_profile(bk), "r", label="Fit") # *10-(-2.5*9)
    # plt.xlabel("Diagonal (pixel)")
    # plt.ylabel("Diagonal profile (rad)")
    # plt.title("Diagonal profile, test image and + 1st fit")
    # plt.show()
    
    # Mark foreground from the first fit into the mask
    bg_mask[phase > bk + noise_delta] = True
            
    # Modify test image to exlude foreground pixels from initial mask
    y0[bg_mask] = bk[bg_mask]
    
    # # Plotting (Optional)
    # plt.figure()
    # plt.plot(xl, diagonal_profile(y0), "g", label="Raw")
    # plt.plot(xl, diagonal_profile(bk), "r", label="Fit")
    # plt.xlabel("Diagonal (pixel)")
    # plt.ylabel("Diagonal profile (rad)")
    # plt.title("Diagonal profile, test image 1st update + 1st fit")
    # plt.show()
    
    # plt.imshow(bk - y0)
    # plt.title("BG - y0")
    # plt.show()
    
    # plt.imshow(bg_mask, cmap='gray')
    # plt.title("Foreground Mask")
    # plt.show()
    
    # # 3D Plotting (Optional)
    # # Create grid for X and Y axes
    # X, Y = np.meshgrid(np.arange(M), np.arange(N))
    # # Create the plot
    # fig = plt.figure(figsize=(10, 6))
    # ax = fig.add_subplot(111, projection='3d')
    # # Surface plot
    # surf = ax.plot_surface(X, Y, y0, cmap='viridis', edgecolor='none')
    # # Optional: add color bar and labels
    # fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
    # ax.set_title("False 3D View of Phase Image")
    # ax.set_xlabel("X pixels")
    # ax.set_ylabel("Y pixels")
    # ax.set_zlabel("Phase (rad)")
    # plt.tight_layout()
    # plt.show()
    
    # plt.close()
    
    c = 0
    # Iteratively refine the background estimate
    while not CheckForLess(abs(y0 - bk).reshape(N * M), 0.001):
        c += 1
        
        # bk_old = bk.copy()
        bg_mask_old = bg_mask.copy()
        
        # Fit a 2D polynomial
        bk, _ = masked_fit_2D(y0, poly_order, bg_mask)
        
        # Update mask: mark new foreground (foreground = True)
        bg_mask[y0 > bk + noise_delta] = True
        
        # Modify test image to exlude foreground pixels from initial mask
        y0[bg_mask] = bk[bg_mask]
        
        # # Ensure y0 is not above the baseline
        # y0 = y0.reshape(N * M)
        # bk = bk.reshape(N * M)
        # y0 = np.minimum(y0, bk + noise_delta)
        # y0 = y0.reshape(N, M)
        # bk = bk.reshape(N, M)
    
        # # Refit a 2D polynomial
        # bk_old = bk.copy()
        # bk, _ = fit_2D(y0, poly_order)
        
        # # Plotting (Optional)
        # plt.figure()
        # plt.plot(xl, diagonal_profile(y0), "g", label="Raw")
        # plt.plot(xl, diagonal_profile(bk), "r", label="Fit")
        # plt.xlabel("Diagonal (pixel)")
        # plt.ylabel("Diagonal profile (rad)")
        # plt.title("Diagonal profile, test image update + fit #"+str(c))
        # plt.show()
        
        # # Plotting (Optional)
        # plt.imshow(bk - y0)
        # plt.title("BG - y0")
        # plt.show()
        
        # # Plotting (Optional)
        # plt.imshow(bg_mask, cmap='gray')
        # plt.title("Foreground Mask")
        # plt.show()
        
        # # 3D Plotting (Optional)
        # # Create grid for X and Y axes
        # X, Y = np.meshgrid(np.arange(M), np.arange(N))
        # # Create the plot
        # fig = plt.figure(figsize=(10, 6))
        # ax = fig.add_subplot(111, projection='3d')
        # # Surface plot
        # surf = ax.plot_surface(X, Y, y0, cmap='viridis', edgecolor='none')
        # # Optional: add color bar and labels
        # fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
        # ax.set_title("False 3D View of Phase Image")
        # ax.set_xlabel("X pixels")
        # ax.set_ylabel("Y pixels")
        # ax.set_zlabel("Phase (rad)")
        # plt.tight_layout()
        # plt.show()
        
        # plt.close()

        ### Check conditions
        # tsh = len(np.where(bk_old - y0 == 0)[0]) * 100 / (N * M)
        mask_change = (np.sum(bg_mask) - np.sum(bg_mask_old)) / np.sum(bg_mask_old)
        mask_ratio = np.sum(bg_mask) / (N * M) * 100
        
        # Check convergence, stop when given percentage of image surface is foreground
        if mask_ratio >= convergence_treshold: #tsh >= convergence_treshold:
            print("Surface percentage reached:",round(mask_ratio*100)/100)
            print("Last mask variation:",round(mask_change*10000)/100,"%")
            print("Number of iterations:",c)
            break
        
        # Stop when the change of mask becomes negligible
        if mask_change < mask_var_trsh:
            print("Mask variantion treshold reached:",round(mask_change*10000)/100,"%")
            print("Convergence:",round(mask_ratio*100)/100)
            print("Number of iterations:",c)
            break
    
    # # Plotting (Optional)
    # plt.imshow(bk - y0)
    # plt.title("BG - y0")
    # plt.show()
    
    print("Background estimated.")
    return bk, bg_mask

def Sequence_BG_detection(path, polynom_order, method, starting_mask, noise_delta, convergence_treshold, mask_var_trsh,bg_frame_num):
    '''
    Detect the "global" background of a sequence of DHM images.
    Uses Interative_BG_detection to estimate the backgrounds of a number of frames.
    Combines them to obtain the final background, i.e., makes an intersection of all (logical AND).
    E.g.: 1=middle frame, 2=first and last frame, 3=first, middle, and last frame,...

    Estimates image background by iterative background detection and polynomial fit to the background.

    Args:
       path (string): File path to the phase image sequence in "LynceeTec format" ("bin", "bnr")
                       - If "bin", select any bin file in the folder. Bin files are single frames of the sequence.
                       - "bnr" - binary file containing the whole sequence.
                       The frames of the sequence (bin or bnr) are 2D arrays of float.
       Poly_Ord (int): Polynomial order for the background fit.
       pix (float): Pixel size.
       method (str): Method of background selection.
                       - "clamp" : Clamp test image (y0) to not exceed background (bk) -> y0 = np.minimum(y0, bk)
                                   This ensures that the evolving test image never exceeds the current baseline.
                       - "mask" ***under construction*** : Create binary mask for evolving foregound/background, to exlude foreground from surface fit.
       starting_mask (numpy.ndarray of boolean): Initial mask to exclude parts from linear fit (pixels to exclude = True),
                                                 If "None" or "False": start with empty mask
       noise_delta (float): Constant to be added to background fit to account for noise in the background
       convergence_treshold (float): Parameter to check convergence; stop when given percentage of image surface is foreground
       mask_var_trsh (float): # Parameter for adaptive tresholding; stop when the change of mask becomes negligible
       bg_frame_num (int): (Needs to be positive integer and smaller or equal to sequence length) Number of frames on which the background extracted. E.g.: 1->middle frame, 2->first and last frame, 3->first, middle, and last frame,...
    Returns:
        global_mask (numpy.ndarray of boolean): Final mask of foregound/background (foreground = True)

    '''
    print("Starting background generation:")
    print(path)
    file_name, file_type = os.path.splitext(path) # get file suffix

    # Depending on file suffix, read sequence lenght and file header
    go_on = True
    if file_type == ".bin":
        print("bin-file detected")
        # get the list of bin files to process
        binfolder=os.path.dirname(path)
        bin_files = []
        bin_files = sorted([f for f in os.listdir(binfolder) if f.endswith(('.bin'))])
        nImages = len(bin_files) # sequence length
        # Get header info and 1st image from sequence
        (image,in_file_header)=binkoala.read_mat_bin(binfolder +os.sep+ bin_files[0])
        hv=str(in_file_header['version'][0])
        end=str(in_file_header['endian'][0])
        hz=str(in_file_header['head_size'][0])
        W=str(in_file_header['width'][0])
        H=str(in_file_header['height'][0])
        pz=in_file_header['px_size'][0]
        hconv=str(in_file_header['hconv'][0])
        uc=str(in_file_header['unit_code'][0])
        
        new_folder = os.path.dirname(binfolder) +os.sep+ os.path.basename(binfolder)+"_BGtZ"
        if not os.path.isdir(new_folder):
            os.mkdir(new_folder)
        if len(os.listdir(new_folder)) != 0:
            answ = messagebox.askquestion('Output folder is not empty!', 'Output folder is not empty.\nDo you want to proceed?')
            print("Output folder is not empty:\n",binfolder)
            if answ == "no":
                print("Stack registration cancelled.")
                go_on = False
        
    elif file_type == ".bnr":
        print("bnr-file detected")
        # Get header info from bnr file
        fileID = open(path, 'rb')
        nImages = np.fromfile(fileID, dtype="i4", count=1)
        nImages = nImages[0]  # sequence length
        w = np.fromfile(fileID, dtype="i4", count=1)
        W=w[0]
        h = np.fromfile(fileID, dtype="i4", count=1)
        H=h[0]
        pz = np.fromfile(fileID, dtype="f4", count=1)
        pz=pz[0]
        wavelength = np.fromfile(fileID, dtype="f4", count=1)
        wavelength=wavelength[0]
        n_1 = np.fromfile(fileID, dtype="f4", count=1)
        n_1 = n_1[0]
        n_2 = np.fromfile(fileID, dtype="f4", count=1)
        n_2 = n_2[0]
        timestamps = [0] * nImages
        for k in range(0,nImages):
            TTT = np.fromfile(fileID, dtype="f4", count=1)
            timestamps[k] = TTT[0]
        # initialize image container
        image = np.zeros((H,W))
        
        frame_size = H * W * 4 # Number of bytes per frame
        data_start = fileID.tell() # Cursor position in the binary file
                    
    else:
        print("ERROR: Wrong file type")
        go_on = False
        return None
        
    if go_on:
        
        # make list of frames to analyse
        bg_frame_list = []
        if bg_frame_num == 1:
            bg_frame_list.append(round(nImages/2)-1) # middle frame
        elif bg_frame_num == 2:
            bg_frame_list.append(0) # first frame
            bg_frame_list.append(nImages-1) # last frame
        elif bg_frame_num == 3:
            bg_frame_list.append(0) # first frame
            bg_frame_list.append(round(nImages/2)-1) # middle frame
            bg_frame_list.append(nImages-1) # last frame
        elif bg_frame_num > 3 and bg_frame_num <= nImages:
            for i in range(bg_frame_num-1):
                
                step = round(nImages/(bg_frame_num-1)*1000)/1000

                bg_frame_list.append(round(i*step))
            
            # if not nImages % bg_frame_num == 0:
            bg_frame_list.append(nImages-1)
                
        else: print("Wrong input for bg_frame_num")
        
        # Initialize mask (foreground = True)
        if isinstance(starting_mask, np.ndarray):
            global_mask = starting_mask.copy()
        else:
            global_mask = np.np.zeros((H, W))
        
        print(bg_frame_list)
        print("Datastart",data_start)

        # Loop through bg_frame_list
        cc = 1
        for frame in bg_frame_list:
            print(f"Processing {cc} of {bg_frame_num}... (frame {frame+1} of {nImages})")
            cc+=1
            
            if file_type == ".bin":
                
                (image,in_file_header)=binkoala.read_mat_bin(binfolder +os.sep+ bin_files[i]) # load frame i
 
            elif file_type == ".bnr":
                print(frame)
                google = int(data_start) + int(i) * int(frame_size)
                print(google)
                fileID.seek(google, 0)
                image = np.fromfile(fileID, dtype="f4", count=H*W).reshape((H, W))
            
            background, bg_mask = Interative_BG_detection(image, polynom_order, pz, method, starting_mask, noise_delta, convergence_treshold, mask_var_trsh)
            global_mask[bg_mask] = True # Make intersection of backgrounds (logical OR on foreground)
            
            # # Plotting (Optional)
            # plt.imshow(image)
            # plt.title("Background of frame "+str(i)+" of "+str(nImages))
            # plt.show()

        if file_type == ".bnr":
            fileID.close
    print("Globel background mask found. Convergence:",round(np.sum(global_mask) / (global_mask.shape[0]*global_mask.shape[1])*100*100)/100)
    return global_mask
    

def Sequence_BGtoZ(path, Poly_Ord, BG_mask):
    '''
    To apply background-to-zero on sequence, with a given background mask.

    Args:
       path (string): File path to the phase image sequence in "LynceeTec format" ("bin", "bnr")
                       - If "bin", select any bin file in the folder. Bin files are single frames of the sequence.
                       - "bnr" - binary file containing the whole sequence.
                       The frames of the sequence (bin or bnr) are 2D arrays of float.
       Poly_Ord (int): Polynomial order for the background fit.
       
       BG_mask (numpy.ndarray of boolean): Initial mask to exclude parts from linear fit (pixels to exclude = True),
                                                 If "None" or "False": start with empty mask
    Returns:
        None

    '''
    print("Starting background correction on sequence:")
    print(path)
    file_name, file_type = os.path.splitext(path) # get file suffix

    # Depending on file suffix, read sequence lenght and file header
    go_on = True
    if file_type == ".bin":
        print("bin-file detected")
        # get the list of bin files to process
        binfolder=os.path.dirname(path)
        bin_files = []
        bin_files = sorted([f for f in os.listdir(binfolder) if f.endswith(('.bin'))])
        nImages = len(bin_files) # sequence length
        # Get header info and 1st image from sequence
        (image,in_file_header)=binkoala.read_mat_bin(binfolder +os.sep+ bin_files[0])
        hv=str(in_file_header['version'][0])
        end=str(in_file_header['endian'][0])
        hz=str(in_file_header['head_size'][0])
        W=str(in_file_header['width'][0])
        H=str(in_file_header['height'][0])
        pz=in_file_header['px_size'][0]
        hconv=str(in_file_header['hconv'][0])
        uc=str(in_file_header['unit_code'][0])
        
        new_folder = os.path.dirname(binfolder) +os.sep+ os.path.name(binfolder)+"_BGtZ"
        if not os.path.isdir(new_folder):
            os.mkdir(new_folder)
        if len(os.listdir(new_folder)) != 0:
            answ = messagebox.askquestion('Output folder is not empty!', 'Output folder is not empty.\nDo you want to proceed?')
            print("Output folder is not empty:\n",binfolder)
            if answ == "no":
                print("Stack registration cancelled.")
                go_on = False
        
    elif file_type == ".bnr":
        print("bnr-file detected")
        # Get header info from bnr file
        fileID = open(path, 'rb')
        nImages = np.fromfile(fileID, dtype="i4", count=1)
        nImages = nImages[0]  # sequence length
        w = np.fromfile(fileID, dtype="i4", count=1)
        W=w[0]
        h = np.fromfile(fileID, dtype="i4", count=1)
        H=h[0]
        pz = np.fromfile(fileID, dtype="f4", count=1)
        pz=pz[0]
        wavelength = np.fromfile(fileID, dtype="f4", count=1)
        wavelength=wavelength[0]
        n_1 = np.fromfile(fileID, dtype="f4", count=1)
        n_1 = n_1[0]
        n_2 = np.fromfile(fileID, dtype="f4", count=1)
        n_2 = n_2[0]
        timestamps = [0] * nImages
        for k in range(0,nImages):
            TTT = np.fromfile(fileID, dtype="f4", count=1)
            timestamps[k] = TTT[0]
        # initialize image container
        image = np.zeros((H,W))
        
        # initialize output file
        new_file = os.path.splitext(path)[0] + "_BGtZ.bnr"
        if os.path.isfile(new_file):
            answ = messagebox.askquestion('Output file exist already!', 'Output file exist already.\nDo you want to overwrite?')
            print("Output file exist already:",new_file)
            if answ == "no":
                print("Stack registration cancelled.")
                go_on = False
        #write meta data to bnr file
        outfileID=open(new_file,'w')
        np.array(nImages, dtype=np.int32).tofile(outfileID)
        np.array(W, dtype=np.int32).tofile(outfileID)
        np.array(H, dtype=np.int32).tofile(outfileID)
        np.array(pz, dtype=np.float32).tofile(outfileID)
        np.array(wavelength, dtype=np.float32).tofile(outfileID)
        np.array(n_1, dtype=np.float32).tofile(outfileID)
        np.array(n_2, dtype=np.float32).tofile(outfileID)
        for k in range(0,nImages):
            np.array(timestamps[k], dtype=np.float32).tofile(outfileID)
            
    else:
        print("ERROR: Wrong file type")
        go_on = False
        return None
        
    if go_on:
        # Loop through frames
        for i in range(nImages):
            print("Processing frame",i+1,"of",nImages,"...")
            # get frame i from sequence
            if file_type == ".bin":
                (image,in_file_header)=binkoala.read_mat_bin(binfolder +os.sep+ bin_files[i])
            elif file_type == ".bnr":
                for k in range(H):
                    image[k,:] = np.fromfile(fileID, dtype="f4", count=W)

            BG, _= masked_fit_2D(image, Poly_Ord, BG_mask)
            
            phase_cor = image - BG
            
            # Write corrected frame to file
            #save new bin-file #k
            if file_type == ".bin":
                outfile=new_folder+'/'+bin_files[i]
                binkoala.write_mat_bin(outfile, phase_cor, W, H, pz, hconv, unit_code=1)
            elif file_type == ".bnr":
                phase_cor.astype(np.float32).tofile(outfileID)
        
        if file_type == ".bnr":
            fileID.close
            outfileID.close
    print("Background correction done.")