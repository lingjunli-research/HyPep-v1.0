# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:29:25 2022

@author: lawashburn
"""

from tkinter import *
from PIL import ImageTk, Image
from tkinter.filedialog import askopenfilename
import sys
from itertools import islice
from subprocess import Popen, PIPE
from textwrap import dedent
from tkinter import messagebox
from tkinter import filedialog
import os
import pickle
from tkinter import ttk
import threading
import time

root = Tk()
root.title('HyPep 1.0')
root.iconbitmap(r"hypep_icon.ico")
root.geometry('700x525')

frame= Frame(root, width= 20, height= 700, bg = '#2F4FAA')
frame.pack(side=LEFT,  anchor=N)

logo_frame= Frame(root, width= 700, height= 700, bg = 'purple')
logo_frame.pack(side=RIGHT, anchor=N)

###
fdr_entry_val =  StringVar()
disc_shuffle_var =  IntVar()
disc_reverse_var = IntVar()
disc_random_var =  IntVar()
disc_hybrid_var = IntVar()
fdr_alg_choice = IntVar()
shuffle_var =  IntVar()
reverse_var = IntVar()
random_var =  IntVar()
hybrid_var = IntVar()
amm_max_prec_charge = StringVar()
amm_max_frag_charge = StringVar()
amm_promex_cutoff_report = StringVar()
amm_prec_err_cutoff_entry_report = StringVar()
amm_frag_cutoff_report = StringVar()
disc_mode_choice = IntVar()
global malc_entry
global disc_alc_entry_val_choice
global motif_path
global hypep_score_cutoff_fdr
global db_path
global out_path
###




discov_mode_decision = []

#SHS mode choice storage
id_fdr_thresh_decision = []
id_reverse_alg_decision =  []
id_hybrid_alg_decision  = []
id_random_alg_decision = []
id_shuffled_ag_decision = []
id_sws_decision = []

#Discovery mode choice storage

discovery_alc_thresh_decision = []
discovery_motif_path_input_decision = []
discovery_malc_thresh_decision = []
discovery_fdr_thresh_decision = []


#AMM mode choice storage
prec_cutoff_thresh_decision = []
frag_cutoff_thresh_decision =  []
promex_score_thresh_decision  = []
prec_max_charge_decison = []
frag_max_charge_decison = []

SWS_OPTIONS = [
"1",
"2",
"3",
"4",
"5",
"6",
"7",
"8",
"9",
"10"
] #etc

LOOPS_OPTIONS = [
"1",
"2",
"3",
"4",
"5",
] #etc

def hide_all_frames():
    logo_frame.pack_forget()
    about_click_frame.pack_forget()
    modules_click_frame.pack_forget()
    files_click_frame.pack_forget()
    run_click_frame.pack_forget()
    help_click_frame.pack_forget()
def about_frame():
    hide_all_frames()
    about_click_frame.pack(side = RIGHT)
def add_modules():
    hide_all_frames()
    modules_click_frame.pack(side = RIGHT)  
def select_files():
    hide_all_frames()
    files_click_frame.pack(side = RIGHT)
def run_hypep():
    hide_all_frames()
    run_click_frame.pack(side = RIGHT)
def help_form():
    hide_all_frames()
    help_click_frame.pack(side = RIGHT)
def get_motif_file_path(): 
    """ Function provides the database full file path."""
    return path_motif_file


def param_form_save():
    #these "get" values are not for variable storage, just for checking if the input values are valid
    #This is the command for when the user selects "save preferences" on the parameter input page
    num_fdr = fdr_entry_val.get()
    num_prec = amm_max_prec_charge.get()
    num_frag = amm_max_frag_charge.get()
    num_prec_error = amm_prec_err_cutoff_entry_report.get()
    num_frag_error = amm_frag_cutoff_report.get()
    num_promex = amm_promex_cutoff_report.get()
    
    if num_fdr.isalpha():
        messagebox.showerror('Invalid FDR Entry','Only numeric values are allowed for FDR!')
    if num_prec.isalpha():
        messagebox.showerror('Invalid Precursor Charge Entry','Only numeric values are allowed for charges!')
    if num_frag.isalpha():
        messagebox.showerror('Invalid Fragment Charge Entry','Only numeric values are allowed for charges!')
    if num_prec_error.isalpha():
        messagebox.showerror('Invalid Precursor Error Entry','Only numeric values are allowed for error!')
    if num_frag_error.isalpha():
        messagebox.showerror('Invalid Fragment Error Entry','Only numeric values are allowed for error!')
    if num_promex.isalpha():
        messagebox.showerror('Invalid Promex Score Entry','Only numeric values are allowed for promex score!')

def param_file_input_save():
    #these "get" values are not for variable storage, just for checking if the input values are valid
    #This is the command for when the user selects "save preferences" on the file selection input page
    space_sample_name = input_sample_file_text.get()
    if space_sample_name.isspace():
        messagebox.showerror('Invalid Sample Name','Sample name cannot have spaces!')

def param_log_export():
    to_discover = disc_mode_choice.get()
    if to_discover == 1:
        discovery_report = 'Perform discovery mode? Yes'
        malc_current = discovery_malc_thresh_decision[-1]
        malc_report = 'Minimum average local confidence value: ' + malc_current
        disc_alc_current = discovery_alc_thresh_decision[-1]
        alc_report = 'Average local confidence value: ' + disc_alc_current
        motif_current = discovery_motif_path_input_decision[-1]
        motif_report = 'Path to motif file: ' + motif_current
        disc_fdr_current = discovery_fdr_thresh_decision[-1]
        disc_fdr_report = 'Discovery mode FDR cutoff: ' + disc_fdr_current
    else:
        discovery_report = 'Perform discovery mode? No'
        malc_report = 'Minimum average local confidence value: N/A'
        alc_report = 'Average local confidence value: N/A'
        motif_report = 'Path to motif file: N/A'
    
    fdr_selection = fdr_alg_choice.get()
    fdr_selection = str(fdr_selection)
    if fdr_selection == '1':
        fdr_algorithm_report = 'FDR algorithm choice: Reverse'
    elif fdr_selection == '2':
        fdr_algorithm_report = 'FDR algorithm choice: Shuffle'
    elif fdr_selection == '3':
        fdr_algorithm_report = 'FDR algorithm choice: Random'
    elif fdr_selection == '4':
        fdr_algorithm_report = 'FDR algorithm choice: Hybrid'
    
    fdr_percent = fdr_entry_val.get()
    fdr_percent = str(fdr_percent)
    fdr_percent_report = 'FDR Threshold (%): ' + fdr_percent
    
    window_size = sws_variable.get()
    window_size = str(window_size)
    window_size_report = 'Sliding window size: ' + window_size
    
    max_p_charge = amm_max_prec_charge.get()
    max_p_charge = str(max_p_charge)
    max_p_charge_report = 'Maximum precursor charge: +' + max_p_charge
    
    max_f_charge = amm_max_frag_charge.get()
    max_f_charge = str(max_f_charge)
    max_f_charge_report = 'Maximum fragment charge: +' + max_f_charge
    
    p_err = amm_prec_err_cutoff_entry_report.get()
    p_err = str(p_err)
    p_err_report = 'Precursor error threshold (ppm): ' + p_err
    
    f_err = amm_frag_cutoff_report.get()
    f_err = str(f_err)
    f_err_report = 'Fragment error threshold (Da): ' + f_err
    
    pro_mex = amm_promex_cutoff_report.get()
    pro_mex = str(pro_mex)
    pro_mex_report = 'Promex threshold: ' + pro_mex
    
    loop = amm_loops_variable.get()
    loop = str(loop)
    loop_report = 'Number of AMM matching loops: ' + loop
    
    sample_ID = input_sample_file_text.get()
    sample_ID_report = 'Sample name: ' + sample_ID
    
    db_path = input_db_text.get()
    db_path_report = 'Database path: ' + db_path
    
    peaks_path = input_peaks_out_text.get()
    peaks_path_report = 'PEAKS output file path: ' + peaks_path
    
    topfd_path = input_topfd_file_text.get()
    topfd_path_report = 'TopFD output file path: ' + topfd_path
    
    rawconverter_path = input_rawconverter_file_text.get()
    rawconverter_path_report = 'RawConverter output file path: ' + rawconverter_path
    
    pp_path = input_pp_file_text.get()
    pp_path_report = 'ProteinProspector theoretical fragment directory path: ' + pp_path
    
    out_path = input_out_dir_text.get()
    out_path_report = 'Ouput directory path: ' + out_path
    
    param_file_entries = [sample_ID_report,discovery_report,malc_report,alc_report,fdr_algorithm_report,
                          fdr_percent_report,window_size_report,max_p_charge_report,max_f_charge_report,
                          p_err_report,f_err_report,pro_mex_report,loop_report,motif_report,db_path_report,
                          peaks_path_report,topfd_path_report,rawconverter_path_report,pp_path_report,
                          out_path_report]
    
    param_file_path = out_path + '\\' + sample_ID + '_parameter_file.txt'
    with open(param_file_path,'a') as f:
        f.writelines('\n'.join(param_file_entries))

    min_alc_pick_path = 'min_alc.pkl'
    discovery_alc_pick_path = 'discovery_alc.pkl'
    motif_path_pick_path = 'motif_path.pkl'
    discovery_fdr_score_pick_path = 'discovery_fdr_score.pkl'
    database_path_pick_path = 'database_path.pkl'
    fdr_pick_path = 'fdr.pkl'
    max_precursor_z_pick_path = 'max_precursor_z.pkl'
    max_fragment_z_pick_path = 'max_fragment_z.pkl'
    precursor_error_pick_path = 'precursor_error.pkl'
    fragment_error_pick_path = 'fragment_error.pkl'
    promex_score_pick_path = 'promex_score.pkl'
    loops_pick_path = 'loops.pkl'
    sample_name_pick_path = 'sample_name.pkl'    
    peaks_path_pick_path = 'peaks_path.pkl'
    topfd_path_pick_path = 'top_fd_path.pkl'
    raw_converter_path_pick_path = 'RawConverter_path.pkl'
    theoretical_fragment_path_pick_path = 'theoretical_fragment_path.pkl'
    output_directory_path_pick_path = 'output_directory.pkl'
    window_size_pick_path = 'window_size.pkl'
    start_discovery_mode_pick_path = 'discovery_mode.pkl'
    fdr_algorithm_pick_path = 'fdr_algorithm.pkl'
     #A new file will be created
    #with open(min_alc_pick_path, 'wb') as file_a:
    #    pickle.dump(malc_current, file_a)
    #with open(discovery_alc_pick_path, 'wb') as file_b:
    #    pickle.dump(disc_alc_current, file_b)
    #with open(motif_path_pick_path, 'wb') as file_c:
    #    pickle.dump(motif_current, file_c)
    #with open(discovery_fdr_score_pick_path, 'wb') as file_d:
    #    pickle.dump(disc_fdr_current, file_d)
    with open(database_path_pick_path, 'wb') as file_e:
        pickle.dump(db_path, file_e)
    with open(fdr_pick_path, 'wb') as file_f:
        pickle.dump(fdr_percent, file_f)
    with open(max_precursor_z_pick_path, 'wb') as file_g:
        pickle.dump(max_p_charge, file_g)
    with open(max_fragment_z_pick_path, 'wb') as file_h:
        pickle.dump(max_f_charge, file_h)
    with open(precursor_error_pick_path, 'wb') as file_i:
        pickle.dump(p_err, file_i)
    with open(fragment_error_pick_path, 'wb') as file_j:
        pickle.dump(f_err, file_j)
    with open(promex_score_pick_path, 'wb') as file_k:
        pickle.dump(pro_mex, file_k)
    with open(loops_pick_path, 'wb') as file_l:
        pickle.dump(loop, file_l)
    with open(sample_name_pick_path, 'wb') as file_m:
        pickle.dump(sample_ID, file_m)
    with open(peaks_path_pick_path, 'wb') as file_n:
        pickle.dump(peaks_path, file_n)
    with open(topfd_path_pick_path, 'wb') as file_o:
        pickle.dump(topfd_path, file_o)
    with open(raw_converter_path_pick_path, 'wb') as file_p:
        pickle.dump(rawconverter_path, file_p)
    with open(theoretical_fragment_path_pick_path, 'wb') as file_q:
        pickle.dump(pp_path, file_q)
    with open(output_directory_path_pick_path, 'wb') as file_q:
        pickle.dump(out_path, file_q)    
    with open(window_size_pick_path, 'wb') as file_s:
        pickle.dump(window_size, file_s)  
    with open(start_discovery_mode_pick_path, 'wb') as file_t:
        pickle.dump(to_discover, file_t)    
    with open(fdr_algorithm_pick_path, 'wb') as file_u:
        pickle.dump(fdr_selection, file_u)  
def disc_mode_param_form():
    disc_alc_entry_report = StringVar()
    input_motif_file_text = StringVar()
    disc_mlc_entry_report = StringVar()  
    hypep_score_cutoff_fdr_report = StringVar()
    
    def set_path_motif_file_field():
        path_motif_file = askopenfilename(filetypes=[("CSV Files","*.csv")]) 
        input_motif_file_text.set(path_motif_file)
    def get_motif_file_path(): 
        """ Function provides the database full file path."""
        return path_motif_file
    
    def stay_on_top():
        win_disc.lift()
        win_disc.after(2000,stay_on_top)
    
    win_disc = Toplevel(root)
    stay_on_top()
    win_disc.title('Discovery Mode Parameter Selection')
    win_disc.iconbitmap(r"hypep_icon.ico")
    win_disc.geometry('500x310')
    disc_parameters = Canvas(win_disc, width= 500, height= 310)
    disc_parameters.pack()
    disc_mode_title = Canvas(disc_parameters, width= 500, height= 50)
    disc_mode_title.pack(side=TOP)
    disc_mode_title.create_text(250, 25, text="Discovery Mode Parameters", fill="#2F4FAA", font=(Font_tuple,20),justify=CENTER)
    disc_mode_title.pack(side=TOP)
    disc_mlc_percent_set = Canvas(disc_parameters, width= 500, height= 50)
    disc_mlc_percent_set.pack(side=TOP)
    disc_mlc_percent_title = Canvas(disc_mlc_percent_set,width=250,height=50)
    disc_mlc_percent_title.pack(side=LEFT)
    disc_mlc_percent_title.create_text(100, 25, text="Minimum ALC", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
    disc_mlc_percent_title.pack()
    disc_mlc_entry = Canvas(disc_mlc_percent_set,width=250,height=50)
    disc_mlc_entry.pack(side=RIGHT)
    disc_mlc_entry_space = Entry(disc_mlc_entry,textvariable=disc_mlc_entry_report)
    disc_mlc_entry_space.pack()
    disc_alc_percent_set = Canvas(disc_parameters, width= 500, height= 50)
    disc_alc_percent_set.pack(side=TOP)
    disc_alc_percent_title = Canvas(disc_alc_percent_set,width=250,height=50)
    disc_alc_percent_title.pack(side=LEFT)
    disc_alc_percent_title.create_text(100, 25, text="Minimum local confidence", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
    disc_alc_percent_title.pack()
    
    hypep_score_cutoff_fdr_set = Canvas(disc_parameters, width= 500, height= 50)
    hypep_score_cutoff_fdr_set.pack(side=TOP)
    
    hypep_score_cutoff_fdr_set_title = Canvas(hypep_score_cutoff_fdr_set,width=310,height=50)
    hypep_score_cutoff_fdr_set_title.pack(side=LEFT)
    hypep_score_cutoff_fdr_set_title.create_text(160, 25, text="HyPep Score Cutoff at", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
    hypep_score_cutoff_fdr_set_title.pack()
    
    hypep_score_cutoff_fdr_set_title2 = Canvas(hypep_score_cutoff_fdr_set,width=250,height=50)
    hypep_score_cutoff_fdr_set_title2.pack(side=RIGHT)
    
    hypep_score_cutoff_fdr_entry = Canvas(hypep_score_cutoff_fdr_set_title2,width=125,height=50)
    hypep_score_cutoff_fdr_entry.pack(side=LEFT)
    hypep_score_cutoff_fdr_entry_space = Entry(hypep_score_cutoff_fdr_entry,textvariable=hypep_score_cutoff_fdr_report)
    hypep_score_cutoff_fdr_entry_space.pack()
    
    hypep_score_cutoff_fdr_text2 = Canvas(hypep_score_cutoff_fdr_set_title2,width=125,height=50)
    hypep_score_cutoff_fdr_text2.pack(side=RIGHT)
    hypep_score_cutoff_fdr_text2.create_text(25, 25, text="% FDR", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
    hypep_score_cutoff_fdr_text2.pack(side=RIGHT)
    
    
    disc_alc_entry = Canvas(disc_alc_percent_set,width=250,height=50)
    disc_alc_entry.pack(side=RIGHT)
    disc_alc_entry_space = Entry(disc_alc_entry,textvariable=disc_alc_entry_report)
    disc_alc_entry_space.pack(side=LEFT)
    motif_file_entry_frame = Canvas(disc_parameters, width= 450, height= 50)
    motif_file_entry_frame.pack(side=TOP)
    motif_file_entry_text_frame = Canvas(motif_file_entry_frame, width= 103, height= 50)
    motif_file_entry_text_frame.pack(side=LEFT)
    motif_file_entry_text_frame.create_text(53, 25, text="Motif File", fill="#2F4FAA", font=(Font_tuple,10),justify=CENTER)
    motif_file_entry_text_frame.pack()
    motif_file_entry_choice_frame = Canvas(motif_file_entry_frame, width= 300, height= 50)
    motif_file_entry_choice_frame.pack(side=RIGHT)
    motif_file_choice_browse_entry = Entry(motif_file_entry_choice_frame, textvariable = input_motif_file_text, width = 45)
    motif_file_choice_browse_entry.pack(side=LEFT)
    motif_file_choice_browse = Button(motif_file_entry_choice_frame, text = "Browse", command = set_path_motif_file_field)
    motif_file_choice_browse.pack(side=RIGHT)
    
    
    
    def disc_mode_form_kill():
        win_disc.destroy()
        disc_alc_entry_val_choice = disc_alc_entry_report.get()
        malc_entry = disc_mlc_entry_report.get()
        hypep_score_cutoff_fdr = hypep_score_cutoff_fdr_report.get()
        motif_path = input_motif_file_text.get()
        discovery_alc_thresh_decision.append(disc_alc_entry_val_choice)
        discovery_motif_path_input_decision.append(motif_path)
        discovery_malc_thresh_decision.append(malc_entry)
        discovery_fdr_thresh_decision.append(hypep_score_cutoff_fdr)
        
    disc_save_button = Button(win_disc, width=218,height=50,text='Save preferences',fg='white',relief='flat',borderwidth=5, 
                         bg='#2F4FAA',font=(Font_tuple,15),command=disc_mode_form_kill)
    disc_save_button.pack(side=BOTTOM)    



im1 = Image.open(r"About.png")
resized_im1 = im1.resize((218,50))
button1 = Button(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA',command=about_frame)
ph_im1 = ImageTk.PhotoImage(resized_im1) # <----------
button1.config(image=ph_im1)
button1.pack(fill=X, expand=1, anchor=N)

im2 = Image.open(r"Add Modules2.png")
resized_im2 = im2.resize((218,50))
button2 = Button(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA',command=add_modules)
ph_im2 = ImageTk.PhotoImage(resized_im2) # <----------
button2.config(image=ph_im2)
button2.pack(fill=X, expand=1, anchor=N)

im3 = Image.open(r"Select Files.png")
resized_im3 = im3.resize((218,50))
button3 = Button(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA' ,command=select_files)
ph_im3 = ImageTk.PhotoImage(resized_im3) # <----------
button3.config(image=ph_im3)
button3.pack(fill=X, expand=1, anchor=N)


button4 = Canvas(frame, width=218,height=50, bg='#2F4FAA',highlightbackground = '#2F4FAA' )
button4.pack(fill=X, expand=1, anchor=N)

im5 = Image.open(r"Run Analysis.png")
resized_im5 = im5.resize((218,50))
button5 = Button(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA',command=run_hypep)
ph_im5 = ImageTk.PhotoImage(resized_im5) # <----------
button5.config(image=ph_im5)
button5.pack(fill=X, expand=1, anchor=N)


button6 = Canvas(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA',highlightbackground = '#2F4FAA' )
button6.pack(fill=X, expand=1, anchor=N, pady=50)

im7 = Image.open(r"Help.png")
resized_im7 = im7.resize((218,50))
button7 = Button(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA',command=help_form)
ph_im7 = ImageTk.PhotoImage(resized_im7) # <----------
button7.config(image=ph_im7)
button7.pack(fill=X, expand=1, anchor=N)


logo = Image.open(r"hypep_logo.png")
resized_logo = logo.resize((500,500))
logo_pos = Label(logo_frame, width=500,height=500,relief='flat',borderwidth=5)
logo2 = ImageTk.PhotoImage(resized_logo) # <----------
logo_pos.config(image=logo2)
logo_pos.pack(expand= YES,fill= BOTH, anchor=CENTER)

Font_tuple = ("Corbel Light", 20)

#about window
about_click_frame= Canvas(root, width= 700, height= 500)
about_click_frame.pack()
title_text_box= Canvas(about_click_frame, width= 500, height= 100)
title_text_box.pack(anchor=N)
title_text_box.create_text(250, 50, text="About HyPep", fill="#2F4FAA", font=(Font_tuple,25),justify=CENTER)
title_text_box.pack(anchor=N)
body_text_box= Canvas(about_click_frame, width= 500, height= 500)
body_text_box.pack(anchor=CENTER)
body_text_box.create_text(225, 100, text="HyPep is a hybrid neuropeptide \nidentification software composed of a sequence homology search \nalgorithm and an accurate \nmass matching algorithm", 
                              fill="#2F4FAA", font=Font_tuple, justify=LEFT, width = 400)
body_text_box.pack(expand= YES,fill= BOTH, anchor=CENTER)

#Add modules window
modules_click_frame= Canvas(root, width= 460, height= 525)
modules_click_frame.pack(side=TOP)
modules_title_frame = Canvas(modules_click_frame,width=460,height=50)
modules_title_frame.pack(side=TOP)
modules_title_frame.create_text(250, 25, text="Add Modules", fill="#2F4FAA", font=(Font_tuple,25),justify=CENTER)
modules_title_frame.pack(side=TOP)
#discovery_frame = Canvas(modules_click_frame,width=460,height=25)
#discovery_frame.pack(side=TOP)
#discovery_frame_title = Canvas(discovery_frame,width=200,height=25)
#discovery_frame_title.pack(side=LEFT)
#discovery_frame_title.create_text(150, 15, text="Discovery Mode?", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
#discovery_frame_title.pack()
#discovery_check_frame = Canvas(discovery_frame,width=200,height=25)
#discovery_check_frame.pack(side=RIGHT)
#discovery_check_frame_checkbox = Checkbutton(discovery_check_frame, variable=disc_mode_choice, command=disc_mode_param_form,font=(Font_tuple,10),width=2,height=3)
#discovery_check_frame_checkbox.pack(side=TOP)
shs_frame = Canvas(modules_click_frame,width=460,height=50)
shs_frame.pack(side=TOP)
shs_frame.create_text(250, 25, text="SHS Parameters", fill="#2F4FAA", font=(Font_tuple,14),justify=CENTER)
shs_frame.pack(side=TOP)
shs_fdr_alg_frame = Canvas(modules_click_frame,width=460,height=25)
shs_fdr_alg_frame.pack(side=TOP)
shs_fdr_alg_title_frame = Canvas(shs_fdr_alg_frame,width=200,height=25)
shs_fdr_alg_title_frame.pack(side=LEFT)
shs_fdr_alg_title_frame.create_text(100, 15, text="FDR Algorithm", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
shs_fdr_alg_title_frame.pack()
shs_fdr_alg_check_frame = Canvas(shs_fdr_alg_frame,width=200,height=25)
shs_fdr_alg_check_frame.pack(side=RIGHT)
reverse_mode_check = Radiobutton(shs_fdr_alg_check_frame, text="Reverse", value='1',variable=fdr_alg_choice,font=(Font_tuple,8))
reverse_mode_check.pack(side=RIGHT)
shuffle_mode_check = Radiobutton(shs_fdr_alg_check_frame, text="Shuffle", value='2',variable=fdr_alg_choice,font=(Font_tuple,8))
shuffle_mode_check.pack(side=RIGHT)
random_mode_check = Radiobutton(shs_fdr_alg_check_frame, text="Random", value='3',variable=fdr_alg_choice,font=(Font_tuple,8))
random_mode_check.pack(side=RIGHT)
hybrid_mode_check = Radiobutton(shs_fdr_alg_check_frame, text="Hybrid", value='4',variable=fdr_alg_choice,font=(Font_tuple,8))
hybrid_mode_check.pack(side=RIGHT)
shs_fdr_frame = Canvas(modules_click_frame,width=460,height=25)
shs_fdr_frame.pack(side=TOP)
shs_fdr_title_frame = Canvas(shs_fdr_frame,width=200,height=25)
shs_fdr_title_frame.pack(side=LEFT)
shs_fdr_title_frame.create_text(100, 15, text="FDR % Threshold", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
shs_fdr_title_frame.pack()
shs_fdr_check_frame = Canvas(shs_fdr_frame,width=200,height=25)
shs_fdr_check_frame.pack(side=RIGHT)
fdr_entry_space = Entry(shs_fdr_check_frame,textvariable=fdr_entry_val)
fdr_entry_space.pack()
shs_sws_frame = Canvas(modules_click_frame,width=460,height=25)
shs_sws_frame.pack(side=TOP)
shs_sws_title_frame = Canvas(shs_sws_frame,width=200,height=25)
shs_sws_title_frame.pack(side=LEFT)
shs_sws_title_frame.create_text(100, 15, text="Sliding window size", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
shs_sws_title_frame.pack()
shs_sws_check_frame = Canvas(shs_sws_frame,width=200,height=25)
shs_sws_check_frame.pack(side=RIGHT)
sws_variable = StringVar(shs_sws_check_frame)
sws_variable.set(SWS_OPTIONS[2]) # default value
sws_dropdown = OptionMenu(shs_sws_check_frame, sws_variable, *SWS_OPTIONS)
sws_dropdown.pack()
amm_frame = Canvas(modules_click_frame,width=460,height=50)
amm_frame.pack(side=TOP)
amm_frame.create_text(250, 25, text="AMM Parameters", fill="#2F4FAA", font=(Font_tuple,14),justify=CENTER)
amm_frame.pack(side=TOP)
amm_prec_maxz_frame = Canvas(modules_click_frame,width=460,height=25)
amm_prec_maxz_frame.pack(side=TOP)
amm_prec_maxz_title_frame = Canvas(amm_prec_maxz_frame,width=200,height=25)
amm_prec_maxz_title_frame.pack(side=LEFT)
amm_prec_maxz_title_frame.create_text(100, 15, text="Max Precursor Charge", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
amm_prec_maxz_title_frame.pack()
amm_prec_maxz_check_frame = Canvas(amm_prec_maxz_frame,width=200,height=25)
amm_prec_maxz_check_frame.pack(side=RIGHT)
prec_maxz = Entry(amm_prec_maxz_check_frame, textvariable=amm_max_prec_charge)
prec_maxz.pack()
amm_frag_maxz_frame = Canvas(modules_click_frame,width=460,height=25)
amm_frag_maxz_frame.pack(side=TOP)
amm_frag_maxz_title_frame = Canvas(amm_frag_maxz_frame,width=200,height=25)
amm_frag_maxz_title_frame.pack(side=LEFT)
amm_frag_maxz_title_frame.create_text(100, 15, text="Max Fragment Charge", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
amm_frag_maxz_title_frame.pack()
amm_frag_maxz_check_frame = Canvas(amm_frag_maxz_frame,width=200,height=25)
amm_frag_maxz_check_frame.pack(side=RIGHT)
amm_prec_maxz_check_frame.pack(side=RIGHT)
frag_maxz = Entry(amm_frag_maxz_check_frame, textvariable=amm_max_frag_charge)
frag_maxz.pack()
amm_prec_err_frame = Canvas(modules_click_frame,width=460,height=25)
amm_prec_err_frame.pack(side=TOP)
amm_prec_err_title_frame = Canvas(amm_prec_err_frame,width=200,height=25)
amm_prec_err_title_frame.pack(side=LEFT)
amm_prec_err_title_frame.create_text(100, 15, text="Precursor Error Cutoff (ppm)", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
amm_prec_err_title_frame.pack()
amm_prec_err_check_frame = Canvas(amm_prec_err_frame,width=200,height=25)
amm_prec_err_check_frame.pack(side=RIGHT)
amm_prec_err_cutoff_entry_space = Entry(amm_prec_err_check_frame, textvariable=amm_prec_err_cutoff_entry_report)
amm_prec_err_cutoff_entry_space.pack()
amm_frag_err_frame = Canvas(modules_click_frame,width=460,height=25)
amm_frag_err_frame.pack(side=TOP)
amm_frag_err_title_frame = Canvas(amm_frag_err_frame,width=200,height=25)
amm_frag_err_title_frame.pack(side=LEFT)
amm_frag_err_title_frame.create_text(100, 15, text="Fragment Error Cutoff (Da)", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
amm_frag_err_title_frame.pack()
amm_frag_err_check_frame = Canvas(amm_frag_err_frame,width=200,height=25)
amm_frag_err_check_frame.pack(side=RIGHT)
amm_frag_err_check_frame_space = Entry(amm_frag_err_check_frame, textvariable=amm_frag_cutoff_report)
amm_frag_err_check_frame_space.pack()
amm_promex_frame = Canvas(modules_click_frame,width=460,height=25)
amm_promex_frame.pack(side=TOP)
amm_promex_title_frame = Canvas(amm_promex_frame,width=200,height=25)
amm_promex_title_frame.pack(side=LEFT)
amm_promex_title_frame.create_text(100, 15, text="Promex Score Cutoff", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
amm_promex_title_frame.pack()
amm_promex_check_frame = Canvas(amm_promex_frame,width=200,height=25)
amm_promex_check_frame.pack(side=RIGHT)
amm_promex_cutoff_entry_space = Entry(amm_promex_check_frame, textvariable=amm_promex_cutoff_report)
amm_promex_cutoff_entry_space.pack()
amm_loops_frame = Canvas(modules_click_frame,width=460,height=25)
amm_loops_frame.pack(side=TOP)
amm_loops_title_frame = Canvas(amm_loops_frame,width=200,height=25)
amm_loops_title_frame.pack(side=LEFT)
amm_loops_title_frame.create_text(100, 15, text="# matching loops", fill="#2F4FAA", font=(Font_tuple,10),justify=RIGHT)
amm_loops_title_frame.pack()
amm_loops_check_frame = Canvas(amm_loops_frame,width=200,height=25)
amm_loops_check_frame.pack(side=RIGHT)
amm_loops_variable = StringVar(amm_loops_check_frame)
amm_loops_variable.set(LOOPS_OPTIONS[4]) # default value
amm_loops_dropdown = OptionMenu(amm_loops_check_frame, amm_loops_variable, *LOOPS_OPTIONS)
amm_loops_dropdown.pack()
##
#save_button_frame= Canvas(modules_click_frame, width= 460, height= 50)
save_button_frame= Canvas(modules_click_frame, width= 460, height= 100)
save_button_frame.pack(side=TOP)
#save_button = Button(save_button_frame, width=450,height=1,text='Save preferences',fg='white',relief='flat',borderwidth=5, 
#                    bg='#2F4FAA',font=(Font_tuple,15),command=param_form_save)
save_button = Button(save_button_frame, width=450,height=4,text='Save preferences',fg='white',relief='flat',borderwidth=5, 
                    bg='#2F4FAA',font=(Font_tuple,15),command=param_form_save)
save_button.pack(side=BOTTOM)
################


#Select files page
sample_name =  StringVar()
path_db = StringVar()
input_db_text = StringVar()
input_out_dir_text = StringVar()
input_peaks_out_text =  StringVar()
input_motif_file_text =  StringVar()
input_topfd_file_text =  StringVar()
input_rawconverter_file_text =  StringVar()
input_pp_file_text =  StringVar()
input_sample_file_text =  StringVar()
shs_select = IntVar()
amm_select = IntVar()
files_click_frame= Canvas(root, width= 450, height= 525)
files_click_frame.pack(side=TOP)

def set_path_database_field():
    path_db = askopenfilename(filetypes=[("CSV Files","*.csv")]) 
    input_db_text.set(path_db)
def get_database_path(): 
    """ Function provides the database full file path."""
    return path_db
def set_path_out_dir_field():
    path_out_dir = filedialog.askdirectory() 
    input_out_dir_text.set(path_out_dir)
def get_out_dir_path(): 
    """ Function provides the database full file path."""
    return path_out_dir
def set_path_peaks_out_field():
    path_peaks_out = askopenfilename(filetypes=[("CSV Files","*.csv")]) 
    input_peaks_out_text.set(path_peaks_out)
def get_peaks_out_path(): 
    """ Function provides the database full file path."""
    return path_peaks_out
def set_path_topfd_file_field():
    path_topfd_file = askopenfilename(filetypes=[("CSV Files","*.csv")]) 
    input_topfd_file_text.set(path_topfd_file)
def get_topfd_file_path(): 
    """ Function provides the database full file path."""
    return path_topfd_file
def set_path_rawconverter_file_field():
    path_rawconverter_file = askopenfilename(filetypes=[("MS2 Files","*.MS2")]) 
    input_rawconverter_file_text.set(path_rawconverter_file)
def get_rawconverter_file_path(): 
    """ Function provides the database full file path."""
    return path_rawconverter_file
def set_path_pp_file_field():
    path_pp_file = filedialog.askdirectory() 
    input_pp_file_text.set(path_pp_file)
def get_pp_file_path(): 
    """ Function provides the database full file path."""
    path_pp_file = path_pp_file

files_title_text_box= Canvas(files_click_frame, width= 450, height= 75)
files_title_text_box.pack(side=TOP)
files_title_text_box.create_text(250, 25, text="Select Files", fill="#2F4FAA", font=(Font_tuple,20),justify=CENTER)
files_title_text_box.pack(side=TOP)


check_postion = Canvas(files_click_frame, width= 450, height= 375)
check_postion.pack(side=TOP)

sample_entry_frame = Canvas(check_postion, width= 450, height= 55)
sample_entry_frame.pack(side=TOP)
sample_entry_text_frame = Canvas(sample_entry_frame, width= 103, height= 50)
sample_entry_text_frame.pack(side=LEFT)
sample_entry_text_frame.create_text(53, 25, text="Sample Name", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
sample_entry_text_frame.pack()
sample_entry_choice_frame = Canvas(sample_entry_frame, width= 300, height= 50)
sample_entry_choice_frame.pack(side=RIGHT)
sample_choice_browse_entry = Entry(sample_entry_choice_frame, textvariable = input_sample_file_text, width = 53)
sample_choice_browse_entry.pack(side=LEFT)

db_entry_frame = Canvas(check_postion, width= 450, height= 55)
db_entry_frame.pack(side=TOP)
db_entry_text_frame = Canvas(db_entry_frame, width= 103, height= 50)
db_entry_text_frame.pack(side=LEFT)
db_entry_text_frame.create_text(53, 25, text="Database", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
db_entry_text_frame.pack()
db_entry_choice_frame = Canvas(db_entry_frame, width= 300, height= 50)
db_entry_choice_frame.pack(side=RIGHT)
db_choice_browse_entry = Entry(db_entry_choice_frame, textvariable = input_db_text, width = 45)
db_choice_browse_entry.pack(side=LEFT)
db_choice_browse = Button(db_entry_choice_frame, text = "Browse", command = set_path_database_field)
db_choice_browse.pack(side=RIGHT)


peaks_out_entry_frame = Canvas(check_postion, width= 450, height= 55)
peaks_out_entry_frame.pack(side=TOP)
peaks_out_entry_text_frame = Canvas(peaks_out_entry_frame, width= 103, height= 50)
peaks_out_entry_text_frame.pack(side=LEFT)
peaks_out_entry_text_frame.create_text(53, 25, text="PEAKS Output File", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
peaks_out_entry_text_frame.pack()
peaks_out_entry_choice_frame = Canvas(peaks_out_entry_frame, width= 300, height= 50)
peaks_out_entry_choice_frame.pack(side=RIGHT)
peaks_out_choice_browse_entry = Entry(peaks_out_entry_choice_frame, textvariable = input_peaks_out_text, width = 45)
peaks_out_choice_browse_entry.pack(side=LEFT)
peaks_out_choice_browse = Button(peaks_out_entry_choice_frame, text = "Browse", command = set_path_peaks_out_field)
peaks_out_choice_browse.pack(side=RIGHT) 

topfd_file_entry_frame = Canvas(check_postion, width= 450, height= 55)
topfd_file_entry_frame.pack(side=TOP)
topfd_file_entry_text_frame = Canvas(topfd_file_entry_frame, width= 103, height= 50)
topfd_file_entry_text_frame.pack(side=LEFT)
topfd_file_entry_text_frame.create_text(53, 25, text="TopFD File", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
topfd_file_entry_text_frame.pack()
topfd_file_entry_choice_frame = Canvas(topfd_file_entry_frame, width= 300, height= 50)
topfd_file_entry_choice_frame.pack(side=RIGHT)
topfd_file_choice_browse_entry = Entry(topfd_file_entry_choice_frame, textvariable = input_topfd_file_text, width = 45)
topfd_file_choice_browse_entry.pack(side=LEFT)
topfd_file_choice_browse = Button(topfd_file_entry_choice_frame, text = "Browse", command = set_path_topfd_file_field)
topfd_file_choice_browse.pack(side=RIGHT)

rawconverter_file_entry_frame = Canvas(check_postion, width= 450, height= 55)
rawconverter_file_entry_frame.pack(side=TOP)
rawconverter_file_entry_text_frame = Canvas(rawconverter_file_entry_frame, width= 103, height= 50)
rawconverter_file_entry_text_frame.pack(side=LEFT)
rawconverter_file_entry_text_frame.create_text(53, 25, text="RawConverter File", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
rawconverter_file_entry_text_frame.pack()
rawconverter_file_entry_choice_frame = Canvas(rawconverter_file_entry_frame, width= 300, height= 50)
rawconverter_file_entry_choice_frame.pack(side=RIGHT)
rawconverter_file_choice_browse_entry = Entry(rawconverter_file_entry_choice_frame, textvariable = input_rawconverter_file_text, width = 45)
rawconverter_file_choice_browse_entry.pack(side=LEFT)
rawconverter_file_choice_browse = Button(rawconverter_file_entry_choice_frame, text = "Browse", command = set_path_rawconverter_file_field)
rawconverter_file_choice_browse.pack(side=RIGHT)

pp_file_entry_frame = Canvas(check_postion, width= 450, height= 55)
pp_file_entry_frame.pack(side=TOP)
pp_file_entry_text_frame = Canvas(pp_file_entry_frame, width= 103, height= 50)
pp_file_entry_text_frame.pack(side=LEFT)
pp_file_entry_text_frame.create_text(53, 25, text="ProteinProspector\nDirectory", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
pp_file_entry_text_frame.pack()
pp_file_entry_choice_frame = Canvas(pp_file_entry_frame, width= 300, height= 50)
pp_file_entry_choice_frame.pack(side=RIGHT)
pp_file_choice_browse_entry = Entry(pp_file_entry_choice_frame, textvariable = input_pp_file_text, width = 45)
pp_file_choice_browse_entry.pack(side=LEFT)
pp_file_choice_browse = Button(pp_file_entry_choice_frame, text = "Browse", command = set_path_pp_file_field)
pp_file_choice_browse.pack(side=RIGHT)

out_dir_entry_frame = Canvas(check_postion, width= 450, height= 55)
out_dir_entry_frame.pack(side=TOP)
out_dir_entry_text_frame = Canvas(out_dir_entry_frame, width= 103, height= 50)
out_dir_entry_text_frame.pack(side=LEFT)
out_dir_entry_text_frame.create_text(53, 25, text="Output directory", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
out_dir_entry_text_frame.pack()
out_dir_entry_choice_frame = Canvas(out_dir_entry_frame, width= 300, height= 50)
out_dir_entry_choice_frame.pack(side=RIGHT)
out_dir_choice_browse_entry = Entry(out_dir_entry_choice_frame, textvariable = input_out_dir_text, width = 45)
out_dir_choice_browse_entry.pack(side=LEFT)
out_dir_choice_browse = Button(out_dir_entry_choice_frame, text = "Browse", command = set_path_out_dir_field)
out_dir_choice_browse.pack(side=RIGHT)

save_button_frame= Canvas(files_click_frame, width= 450, height= 75)
save_button_frame.pack(side=BOTTOM)
save_button = Button(save_button_frame, width=450,height=2,text='Save preferences',fg='white',relief='flat',borderwidth=5, 
                    bg='#2F4FAA',font=(Font_tuple,15),command=param_file_input_save)
save_button.pack(side=BOTTOM)

def RunHyPep():
    time.sleep(5)
    os.system('python command_center.py')

def RunSelectScripts():
    time.sleep(5)
    os.system('python command_center_after_prot_pros.py')

def combine_funcs(*funcs):
    def combined_func(*args, **kwargs):
        for f in funcs:
            f(*args, **kwargs)
    return combined_func



#Run page
run_click_frame= Canvas(root, width= 450, height= 525)

pb = ttk.Progressbar(
    run_click_frame,
    orient='horizontal',
    mode='indeterminate',
    length=280
)
pb.pack(side=BOTTOM)

#save_button = Button(run_click_frame, width=450,height=1,text='Run Analysis!',fg='white',relief='flat',borderwidth=5, 
#                    bg='#2F4FAA',font=(Font_tuple,15),command = combine_funcs(param_log_export,RunHyPep))

save_button = Button(run_click_frame, width=450,height=1,text='Run Analysis!',fg='white',relief='flat',borderwidth=5, 
                    bg='#2F4FAA',font=(Font_tuple,15),command = threading.Thread(target=combine_funcs(param_log_export,pb.start,RunHyPep)).start)
save_button.pack(side=BOTTOM)

#save_button = Button(run_click_frame, width=450,height=1,text='Run Analysis!',fg='white',relief='flat',borderwidth=5, 
#                    bg='#2F4FAA',font=(Font_tuple,15),command = combine_funcs(param_log_export,RunHyPep))

advanced_button = Button(run_click_frame, width=450,height=1,text='Advanced settings',fg='white',relief='flat',borderwidth=5, 
                    bg='#2F4FAA',font=(Font_tuple,15),command = threading.Thread(target=combine_funcs(param_log_export,pb.start,RunSelectScripts)).start)
advanced_button.pack(side=BOTTOM)


#Help page
help_click_frame= Canvas(root, width= 450, height= 525)
help_click_frame.pack(side=TOP)
help_text_box= Canvas(help_click_frame, width= 500, height= 100)
help_text_box.pack(anchor=N)
help_text_box.create_text(250, 50, text="Help", fill="#2F4FAA", font=(Font_tuple,25),justify=CENTER)
help_text_box.pack(anchor=N)
help_body_text_box= Canvas(help_click_frame, width= 500, height= 500)
help_body_text_box.pack(anchor=CENTER)
help_body_text_box.create_text(225, 100, text="Please direct all HyPep questions to Lauren Fields, lawashburn@wisc.edu", 
                              fill="#2F4FAA", font=Font_tuple, justify=LEFT, width = 400)
help_body_text_box.pack(expand= YES,fill= BOTH, anchor=CENTER)
root.mainloop()