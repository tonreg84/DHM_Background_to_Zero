# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 13:52:10 2025

@author: ge1582
"""
import Background_tools as BGt
import tkinter as tk
from tkinter import messagebox, filedialog
from os import path
from numpy import uint8
from PIL import Image

# Set parameters
# A) Excludes pixels / make noisy mask
low_treshold_para = 2 # parameter to exclude pixel with very low values (supposedly wrong due to phase jumps)
noise_treshold_deg = 15 # Pixels which vary between two frames more that this value, will be added to the exclusion mask.
noise_treshold = noise_treshold_deg/360*2*3.14
extend_para = 2 # Number of pixels by which the final mask will be extendet
# B) Iterative BG estimation single frame
polynom_order = 3 # Polynomial order for the background fit.
method = "clamp"
noise_delta = 10 /360*2*3.14 # Some tolerance to noise in the background. E.g. diffraction fringes from dust particel in the optical path.
convergence_treshold = 80
mask_var_trsh = 0.01 # Mask variantion treshold
# C) Make global BG mask 
bg_frame_num = 10 # Number of frames on which the background extracted. Final background is intersection of all (logical AND). E.g.: 1=middle frame, 2=first and last frame, 3=first, middle, and last frame,...

# Open DHM sequence
'''
Chose file path of the DHM sequence in "LynceeTec format" ("bin", "bnr")
- If "bin", select any bin file in the folder. Bin files are single frames of the sequence.
- "bnr" - binary file containing the whole sequence.
The frames of the sequence (bin or bnr) are 2D arrays of float.
'''
file_path = filedialog.askopenfilename(title="Select a DHM sequence")

# Proceed
# A) make noisy mask
noisy_mask = BGt.Make_noisy_mask(file_path,low_treshold_para,noise_treshold,extend_para) # make mask to excludes pixels with very low phase values and noisy pixels by noise detection by phase treshold.

# Save noisy mask to png file
file_name, file_extension = path.splitext(file_path)
mask_file = file_name + "_NoisyMask_" + str(noise_treshold_deg) + ".png"
Image.fromarray(noisy_mask.astype(uint8)*255, mode="L").save(mask_file)

# B)
# ...

# C)
# # Read initial mask from png file
# test_mask = Image.open(mask_file).convert("L")
# arr = array(test_mask) # Convert to numpy
# noisy_mask = arr > 0 # Make boolean array by threshold: >0 → True, else False
global_BG = BGt.Sequence_BG_detection(file_path, polynom_order, method, noisy_mask, noise_delta, convergence_treshold, mask_var_trsh, bg_frame_num)
# Save global BG mask to png file
glob_mask_file = file_name + "_globBGmask.png"
Image.fromarray(global_BG.astype(uint8)*255, mode="L").save(glob_mask_file)

# D) Apply background-to-zero on sequence by fitting background frame by frame with a given constant background mask.
BGt.Sequence_BGtoZ(file_path, polynom_order, global_BG)