# -*- coding: utf-8 -*-
"""
DHM Background to Zero
Autor: tonreg, team UMI, CNP-CHUV Lausanne
Version: 20260108

Programme to flatten phase images obtained with Digital Holographic Microscopes. Find image background. Polynomial fit to background. Subtract fitted background from image.

Set parameters:
    
A) Excludes pixels / make noisy mask
low_treshold_para = 115 in degrree, parameter to exclude pixel with very low values (supposedly wrong due to phase jumps)
noise_treshold = 15 in degree, pixels which vary between two frames more that this value, will be added to the exclusion mask.
extend_para = 2 # Number of pixels by which the final mask will be extendet

B) Iterative BG estimation single frame
polynom_order = 3 # Polynomial order for the background fit.
noise_delta = 10 in degree, Some tolerance to noise in the background. E.g. diffraction fringes from dust particel in the optical path.
convergence_treshold = 80%
mask_var_trsh = 0.01 # Mask variantion treshold

C) Make global BG mask 
bg_frame_num = 10 # Number of frames on which the background extracted. Final background is intersection of all (logical AND). E.g.: 1=middle frame, 2=first and last frame, 3=first, middle, and last frame,...

Open DHM sequence:
Chose file path of the DHM sequence in "LynceeTec format" ("bin", "bnr")
- If "bin", select any bin file in the folder. Bin files are single frames of the sequence.
- "bnr" - binary file containing the whole sequence.
"""
import sys
from threading import Thread
import Background_tools as BGt
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
from os import path
from numpy import uint8, zeros, ndarray, argwhere
from PIL import Image, ImageTk

class ConsoleRedirector(object):
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, message):
        self.text_widget.insert(tk.END, message)
        self.text_widget.see(tk.END)

    def flush(self):
        pass

def manage_widgets(frame,state):
    for child in frame.winfo_children():
        try:
            child.configure(state=state)
        except tk.TclError:
            print(child)
            pass   # some widgets (e.g. frames, labels) don't have a "state"

# Parameters:
prog_vers = "20260106"
values = {} # dictionnary for all parameter values
noisy_mask = None
global_mask = None

displaysize = 500

para_list = [
    "low_treshold_para",
    "noise_treshold",
    "extend_para",
    "polynom_order",
    "noise_delta",
    "convergence_treshold",
    "mask_var_trsh",
    "bg_frame_num"
]
NM_list = ["low_treshold_para","noise_treshold","extend_para"]
BG_list = ["polynom_order","noise_delta","convergence_treshold","mask_var_trsh","bg_frame_num"]

text_dict = {}
text_dict["low_treshold_para"] = "Low treshold"
text_dict["noise_treshold"] = "Noise treshold"
text_dict["extend_para"] = "Mask extension"
text_dict["polynom_order"] = "Polynomial order"
text_dict["noise_delta"] = "Background noise tolerance"
text_dict["convergence_treshold"] = "Convergence treshold"
text_dict["mask_var_trsh"] = "Mask variantion treshold"
text_dict["bg_frame_num"] = "Number of frames for global background"

hint_dict = {}
hint_dict["low_treshold_para"] = "degree, float > 0"
hint_dict["noise_treshold"] = "degree, float > 0"
hint_dict["extend_para"] = "pixel, int > 0"
hint_dict["polynom_order"] = "1, 2, 3, 4, 5, or 6"
hint_dict["noise_delta"] = "degree, float > 0"
hint_dict["convergence_treshold"] = "0-100 %"
hint_dict["mask_var_trsh"] = "float > 0"
hint_dict["bg_frame_num"] = "int > 0"

type_dict = {}
type_dict["low_treshold_para"] = "flt"
type_dict["noise_treshold"] = "flt"
type_dict["extend_para"] = "int"
type_dict["polynom_order"] = "int"
type_dict["noise_delta"] = "flt"
type_dict["convergence_treshold"] = "flt"
type_dict["mask_var_trsh"] = "flt"
type_dict["bg_frame_num"] = "int"

default_value_dict = {}
default_value_dict["low_treshold_para"] = 115
default_value_dict["noise_treshold"] = 15
default_value_dict["extend_para"] = 2
default_value_dict["polynom_order"] = 3
default_value_dict["noise_delta"] = 10
default_value_dict["convergence_treshold"] = 80
default_value_dict["mask_var_trsh"] = 0.01
default_value_dict["bg_frame_num"] = 10

para_dict = {}
for para in para_list:
    para_dict[para] = {}
    para_dict[para]["text"] = text_dict[para]
    para_dict[para]["hint"] = hint_dict[para]
    para_dict[para]["type"] = type_dict[para]
    para_dict[para]["defval"] = default_value_dict[para]


# Validation rules
def validate_inputs(values):
    errors = []
    
    if "file_path" in values.keys():
        # file_path must not be empty
        if not path.isfile(values["file_path"]):
            errors.append("Valid file path must be selected")
            
    if "low_treshold_para" in values.keys():
        # low_treshold_para: float > 0
        try:
            if float(values["low_treshold_para"]) <= 0:
                errors.append("low_treshold_para must be > 0")
        except ValueError:
            errors.append("low_treshold_para must be a float")
    if "noise_treshold" in values.keys():
        # noise_treshold: float > 0
        try:
            if float(values["noise_treshold"]) <= 0:
                errors.append("noise_treshold must be > 0")
        except ValueError:
            errors.append("noise_treshold must be a float")
    if "extend_para" in values.keys():
        # extend_para: int > 0
        try:
            if int(values["extend_para"]) <= 0:
                errors.append("extend_para must be a positive integer")
        except ValueError:
            errors.append("extend_para must be an integer")
    if "polynom_order" in values.keys():
        # polynom_order: int 1-6
        try:
            po = int(values["polynom_order"])
            if po < 1 or po > 6:
                errors.append("polynom_order must be between 1 and 6")
        except ValueError:
            errors.append("polynom_order must be an integer")
    if "noise_delta" in values.keys():
        # noise_delta: float > 0
        try:
            if float(values["noise_delta"]) <= 0:
                errors.append("noise_delta must be > 0")
        except ValueError:
            errors.append("noise_delta must be a float")
    if "convergence_treshold" in values.keys():
        # convergence_treshold: percentage 0-100
        try:
            ct = float(values["convergence_treshold"])
            if ct < 0 or ct > 100:
                errors.append("convergence_treshold must be between 0 and 100 (%)")
        except ValueError:
            errors.append("convergence_treshold must be a float (percentage)")
    if "mask_var_trsh" in values.keys():
        # mask_var_trsh: float > 0
        try:
            if float(values["mask_var_trsh"]) <= 0:
                errors.append("mask_var_trsh must be > 0")
        except ValueError:
            errors.append("mask_var_trsh must be a float")
    if "bg_frame_num" in values.keys():
        # bg_frame_num: int > 0
        try:
            if int(values["bg_frame_num"]) <= 0:
                errors.append("bg_frame_num must be an integer > 0")
        except ValueError:
            errors.append("bg_frame_num must be an integer")

    return errors

def browse_file():
    global values, noisy_mask, global_mask
    file_path = filedialog.askopenfilename()
    if file_path:
        file_entry.delete(0, tk.END)
        file_entry.insert(0, file_path)
        values = {}
        noisy_mask = None
        global_mask = None
        
def submit(btnstrg):
    thread = Thread(target = main_proceedure, args = (btnstrg,), daemon=True)
    thread.start()

def main_proceedure(btnstrg):
    global values, noisy_mask, global_mask
    
    if btnstrg == "sNM": # Save noisy mask to png file
        if not isinstance(noisy_mask, ndarray):
            messagebox.showerror(":Error", "Make noisy mask first!")
        else:
            # Make file name depending on file format
            file_name, file_type = path.splitext(values["file_path"]) # get file suffix
            if file_type == ".bin":
                mask_file = path.dirname(path.dirname(values["file_path"])) + "/" + path.basename(path.splitext(values["file_path"])[0]) + "_NoisyMask.png"
            elif file_type == ".bnr":
                mask_file = file_name + "_NoisyMask.png"
            
            Image.fromarray(noisy_mask.astype(uint8)*255, mode="L").save(mask_file)
        
    if btnstrg == "sGM":
        # Save global mask to png file
        if not isinstance(global_mask, ndarray):
            messagebox.showerror(":Error", "Make global mask first!")
        else:
            # Make file name depending on file format
            file_name, file_type = path.splitext(values["file_path"]) # get file suffix
            if file_type == ".bin":
                mask_file = path.dirname(path.dirname(values["file_path"])) + "/" + path.basename(path.splitext(values["file_path"])[0]) + "_GlobalMask.png"
            elif file_type == ".bnr":
                mask_file = file_name + "_GlobalMask.png"
            
            Image.fromarray(global_mask.astype(uint8)*255, mode="L").save(mask_file)

    if btnstrg == "NM":
        
        manage_widgets(main_frame,"disabled")
        
        # Check input:
        for para in NM_list:
            values[para] = para_dict[para]['entry'].get()
        values["file_path"] = file_entry.get()
        errors = validate_inputs(values)
        if errors:
            messagebox.showerror("Validation Error", "\n".join(errors))
        else:
            print("Input parameters validated!")
            
            # Convert input from string to float or int
            for para in NM_list:
                if para_dict[para]['type'] == "flt":
                    values[para] = float(values[para])
                if para_dict[para]['type'] == "int":
                    values[para] = int(values[para])
           
            # Proceed
            # A) make noisy mask
            noisy_mask = BGt.Make_noisy_mask(values["file_path"],values["low_treshold_para"]/360*2*3.14,values["noise_treshold"]/360*2*3.14,values["extend_para"]) # make mask to excludes pixels with very low phase values and noisy pixels by noise detection by phase treshold.
            
            update_canvas_from_array(canvasNM,noisy_mask)
        
        manage_widgets(main_frame,"normal")
    
    elif btnstrg == "GM":
        
        manage_widgets(main_frame,"disabled")
        
        if not isinstance(noisy_mask, ndarray):
            messagebox.showerror(":Error", "Make noisy mask first!")
        else:
            # Check input:
            for para in BG_list:
                values[para] = para_dict[para]['entry'].get()
            values["file_path"] = file_entry.get()
            errors = validate_inputs(values)
            if errors:
                messagebox.showerror("Validation Error", "\n".join(errors))
            else:
                print("Input parameters validated!")
                
                # Convert input from string to float or int
                for para in BG_list:
                    if para_dict[para]['type'] == "flt":
                        values[para] = float(values[para])
                    if para_dict[para]['type'] == "int":
                        values[para] = int(values[para])

                '''
                Method of background selection. 
                                - "clamp" : Clamp test image (y0) to not exceed background (bk) -> y0 = np.minimum(y0, bk)
                                            This ensures that the evolving test image never exceeds the current baseline.
                                - "mask" : Updat binary mask for every itineration, to exlude evolving foreground from surface fit.
                '''
                method = "mask" #TODO : currently only using method "mask"...
                
                global_mask = BGt.Sequence_BG_detection(values["file_path"],values["polynom_order"],method,noisy_mask,values["noise_delta"]/360*2*3.14,values["convergence_treshold"],values["mask_var_trsh"],values["bg_frame_num"])
                update_canvas_from_array(canvasGM,global_mask)
        
        manage_widgets(main_frame,"normal")

    elif btnstrg == "BGtZ":
   
        manage_widgets(main_frame,"disabled")
        
        goo = True
        if not isinstance(noisy_mask, ndarray):
            messagebox.showerror(":Error", "Make noisy mask and a global background mask first!")
            goo = False
        else:
            if not isinstance(global_mask, ndarray):
                messagebox.showerror(":Error", "Make a global background mask first!")
                goo = False
        
        if goo:
            
            file_name, file_extension = path.splitext(values["file_path"])
           
            '''
            Method of background selection. 
                            - "clamp" : Clamp test image (y0) to not exceed background (bk) -> y0 = np.minimum(y0, bk)
                                        This ensures that the evolving test image never exceeds the current baseline.
                            - "mask" : Updat binary mask for every itineration, to exlude evolving foreground from surface fit.
            '''
            method = "mask" #TODO : currently only using method "mask"...
                
            # Save configuration and global mask to file
            # Make file name depending on file format
            file_name, file_type = path.splitext(values["file_path"]) # get file suffix
            if file_type == ".bin":
                config_file = path.dirname(path.dirname(values["file_path"])) + "/" + "BGtZparas.txt"
            elif file_type == ".bnr":
                config_file = file_name + "_BGtZparas.txt"

            with open(config_file, 'w') as fileID:
            
                fileID.write("Background to Zero configuration:\n\n")
                fileID.write('BGtZ done with programme "BGt DHM Background to Zero" (BGtZero_user_GUI.py)\n')
                fileID.write("Autor: tonreg, team UMI, CNP-CHUV Lausanne\n")
                fileID.write(f"Version: {prog_vers}\n")
                
                fileID.write("Parameter:\n")
                for para in para_list:
                    fileID.write(f"{para}: {values[para]}\n")

                fileID.write("\nGlobal background mask:\n")
                fileID.write(f"Image height and width: {global_mask.shape}\n")
                fileID.write("Background mask coordinates (row, column):\n")
                
                BG_coords = argwhere(global_mask == False)
                
                for i in range(len(BG_coords)):
                    fileID.write(f"({BG_coords[i][0]},{BG_coords[i][1]}),")
                fileID.write(f"({BG_coords[-1][0]},{BG_coords[-1][1]})")
        
            # Flatten every frame by fitting background using global background mask.
            BGt.Sequence_BGtoZ(values["file_path"],values["polynom_order"],global_mask)

        manage_widgets(main_frame,"normal")

def update_canvas_from_array(canvas,array):
    global img_tk_list
    
    true_color = (253, 231, 36)
    false_color = (68, 1, 84)
    
    # Create RGB image from boolean array
    h, w = array.shape
    rgb_array = zeros((h, w, 3), dtype=uint8)
    rgb_array[array] = true_color
    rgb_array[~array] = false_color

    img = Image.fromarray(rgb_array, mode="RGB")
    img = img.resize((600, 600), Image.NEAREST)
    img_tk = ImageTk.PhotoImage(img)
    
    # Keep reference to avoid garbage collection
    img_tk_list.append(img_tk)
    
    canvas.config(width=600, height=600)
    canvas.create_image(0, 0, anchor="nw", image=img_tk)


root = tk.Tk()
root.title("Background to Zero")
root.iconbitmap("road_roller.ico")
root.geometry("1750x700")


# Create the Notebook (tab container)
notebook = ttk.Notebook(root)
notebook.pack(fill="both", expand=True)

# Main page
page1 = ttk.Frame(notebook)
notebook.add(page1, text="Main page")

# Console page
page2 = ttk.Frame(notebook)
notebook.add(page2, text="Console")


# Main page layout

# Do stuff to add scroll bars to window:
background_frame = tk.Frame(page1)
background_frame.pack(fill=tk.BOTH, expand=1)
# Canvas
my_canvas = tk.Canvas(background_frame)
my_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
# Vertical scrollbar
v_scrollbar = ttk.Scrollbar(background_frame, orient=tk.VERTICAL, command=my_canvas.yview)
v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
# Horizontal scrollbar
h_scrollbar = ttk.Scrollbar(page1, orient=tk.HORIZONTAL, command=my_canvas.xview)
h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
# Configure canvas for scrollbars
my_canvas.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
my_canvas.bind('<Configure>', lambda e: my_canvas.configure(scrollregion=my_canvas.bbox("all")))
# Frame inside canvas
main_frame = tk.Frame(my_canvas)
# Add frame to canvas
my_canvas.create_window((0, 0), window=main_frame, anchor="nw")


# File selection row
tk.Label(main_frame, text="Select File").grid(row=0, column=0, sticky="w", padx=5, pady=5)
file_entry = tk.Entry(main_frame, width=30)
file_entry.grid(row=0, column=1, padx=0, pady=5)
browse_btn = tk.Button(main_frame, text="Browse", command=browse_file)
browse_btn.grid(row=0, column=2, padx=5, pady=5)

# Noisy mask parameters
tk.Label(main_frame, text="Make-Noisy-Mask Parameters:").grid(row=1, column=0, sticky="w", padx=5, pady=5)

for i, para in enumerate(NM_list):
    tk.Label(main_frame, text=f"{para_dict[para]['text']} ({para_dict[para]['hint']})").grid(row=i+2, column=0, sticky="w", padx=5, pady=5)
    entry = tk.Entry(main_frame, width=10)
    entry.grid(row=i+2, column=1, padx=5, pady=5)
    para_dict[para]['entry'] = entry
    entry.insert(0, para_dict[para]['defval'])

NM_button = tk.Button(main_frame, text="Make Noisy Mask", command=lambda:submit("NM"))
NM_button.grid(row=len(NM_list)+1, column=2, pady=10)

saveNM_button = tk.Button(main_frame, text="Save Noisy Mask as png", command=lambda:submit("sNM"))
saveNM_button.grid(row=len(NM_list)+2, column=2, padx=5, pady=5)

NM_row_offset = len(NM_list)+3

#Background-Estimation Parameters
tk.Label(main_frame, text="Background-Estimation Parameters:").grid(row=NM_row_offset, column=0, sticky="w", padx=5, pady=5)

for i, para in enumerate(BG_list):
    tk.Label(main_frame, text=f"{para_dict[para]['text']} ({para_dict[para]['hint']})").grid(row=NM_row_offset+1+i, column=0, sticky="w", padx=5, pady=5)
    entry = tk.Entry(main_frame, width=10)
    entry.grid(row=NM_row_offset+1+i, column=1, padx=5, pady=5)
    para_dict[para]['entry'] = entry
    entry.insert(0, para_dict[para]['defval'])
    
BG_button = tk.Button(main_frame, text="Make Global Mask", command=lambda:submit("GM"))
BG_button.grid(row=NM_row_offset + len(BG_list) + 1, column=2, pady=10)

saveGM_button = tk.Button(main_frame, text="Save Global Mask as png", command=lambda:submit("sGM"))
saveGM_button.grid(row=NM_row_offset + len(BG_list) + 2, column=2, padx=5, pady=5)

BG_row_offset = row=NM_row_offset + len(BG_list) + 3

# tk.Label(main_frame, text="   ").grid(row=BG_row_offset, column=0, sticky="w", padx=5, pady=5)

# tk.Label(main_frame, text="   ").grid(row=BG_row_offset + 3, column=0, sticky="w", padx=5, pady=5)

submit_btn = tk.Button(main_frame, text="Background to ZERO", command=lambda:submit("BGtZ"), height = 3, width = 30)
submit_btn.grid(row=NM_row_offset + len(BG_list) + 2, column=0, padx = 30, pady=0)

# Canvas for numpy boolean array visualization
canvasNM = tk.Canvas(main_frame, bg="white", width=displaysize, height=displaysize)
canvasNM.grid(row=0, column=3, rowspan=15, padx=10, pady=10)

# Canvas for numpy boolean array visualization
canvasGM = tk.Canvas(main_frame, bg="white", width=displaysize, height=displaysize)
canvasGM.grid(row=0, column=4, rowspan=15, padx=10, pady=10)

# Keep references to images
img_tk_list = []

# -------------------------------------------------------------------
# Console text field
console_box = scrolledtext.ScrolledText(page2, width=200, height=30)
console_box.grid(row=1, column=0, sticky="we")

# redirect output
sys.stdout = ConsoleRedirector(console_box)
sys.stderr = ConsoleRedirector(console_box)


root.mainloop()