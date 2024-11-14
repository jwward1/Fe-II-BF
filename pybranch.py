"""
A script for calculating branching fractions from sets of intensity ratios
measured in various emission spectra.
"""
import logging
import math
from branchgrep import grep_open
from branchrc import set_parameters
from glob import glob
from itertools import product
from os import popen
from tkinter import *
from tkinter.scrolledtext import ScrolledText
from tkinter import simpledialog

calc_file = "FeII_waveno.E1"

def get_calculations(E1):
	  #
	  # This function gets the calculated A values for calculating a residual
	  #
	  
    for line in open(calc_file, 'r').readlines():
        line = line.split()
        waveno = line[0]
        upper_level = line[10]
        E1[upper_level,waveno] = float( line[1])
    return(E1)
     

class BranchingFractionCalc(Frame):
    """
    The frame the user interacts with.
    """
    # Set other initial variables
    all_spectrum_files = glob(set_parameters()[8])
    life_unc = set_parameters()[7]

    def __init__(self, parent=None):
        """
        Description:
            Initializes frame, sets title and resolution.

        Inputs:
            parent - Frame containing Results window. Default is None.
        """
        Frame.__init__(self, parent=None)
        self.master.title("Branching Fraction Calculator")
        self.master.geometry("1200x700")
        self.master.option_add('*Font','mono 8')
        self.pack(expand=YES, fill=BOTH)
        self.make_widgets()
        self.get_lifetimes()
       
    def get_lifetimes(self):
        """
        Description:
            Imports level lifetimes from .lev files.
        """
        self.lifetimes = {}
        self.upper_value = {}
        lifetime_files = glob('*.lev')


        for lifetime_file in lifetime_files:

            for line in open(lifetime_file, 'r').readlines():
                line = line.split()
                upper_level = line[6]
                self.upper_value[upper_level] = line[2]
                try: lifetime = eval(line[5])
                except: lifetime = None
                self.lifetimes[upper_level] = lifetime
                self.upper_value[upper_level] = line[2]

    def log(self, text):
        """
        Description:
            Writes text to log file and displays line in ScrolledText on main display.

        Inputs:
            text - text to be written to log and put on main display.
        """
        logging.info(text)
        self.display_text.insert('end', text+'\n')

    def set_output(self, event):
        """
        Description:
            Sets log_name to entry contents after pressing Enter key (event).
        """
        log_name = self.file_name_entry.get()
        logging.basicConfig(filemode='a', filename=log_name+'.log', level=logging.INFO, format='%(message)s')
        self.log(f'Log file name set to {log_name}.log\n')
        event.widget["foreground"] = "Black"

    def get_upper_level(self, event):
        """
        Description:
            Sets self.upper_level to entry contents after pressing Enter key (event).
            First displays table of previously identified transitions and lifetime
            for given upper level, then displays the wavenumbers of transitions,
            lower levels, SNRs, Intensities, LS configurations, and notes for observed lines.
        """
        self.upper_level = self.upper_level_entry.get()
        event.widget["foreground"] = "Black"
        self.log(f'Upper level: {self.upper_level}')
        self.show_identified_lines()
        self.get_snrs_and_ints()
        self.show_values()

    def show_identified_lines(self):
        """
        Description
            Finds and displays all previously identified lines associated with
            upper level in the .sort file.
        """

        self.known_lines = []
        energy_levels_file = set_parameters()[6]

        # Check level file for given level key.
        for line in open(energy_levels_file, 'r').readlines():
            line = line.split()

            if len(line) > 0:

                level_file_key = line[6]

                # If level Key is found -> return upper level conf key
                if(level_file_key == self.upper_level):

                    upper_energy_key = line[6]

                    # JJL
                    self.known_lines = grep_open(upper_energy_key, self.id_lines_file)
                    # self.known_lines = popen(f"grep {upper_energy_key} {self.id_lines_file}", 'r').readlines()

                    break

        self.log(f'Lifetime of upper level: {self.lifetimes[self.upper_level]} ns\n')
        self.log(f'Identified transitions from upper level {self.upper_level}:')
        self.log(''.join(['-']*67))
        self.log(f'| Intensity | Wavenumber | L Conf       | U Conf        | Notes   |')
        self.log(''.join(['-']*67))

        for transition in self.known_lines:

            transition = transition.split()
            transition_key = transition[7]

            if (transition_key == upper_energy_key):
                wavenumber = eval(transition[4])
                l_conf     = transition[6]
                u_conf     = transition[7]
                try: 
                    note = (transition[0],transition[1],transition[8],transition[9])
                except: 
                    note=(transition[0],transition[1])
                try:
                    self.log(f'| {note[0]:>1s} {note[1]:>5s}   | {wavenumber:>10.3f} | {l_conf:>12s} | {u_conf:>12s}  | {note[2]:>1s} {note[3]:>2s}    |') 
                except:
                    self.log(f'| {note[0]:>1s} {note[1]:>5s}   | {wavenumber:>10.3f} | {l_conf:>12s} | {u_conf:>12s}  |         |')

        self.log(''.join(['-']*67)+'\n')

    def get_snrs_and_ints(self):
        """
        Description:
            Creates wavenumbers list, and SNRs, Intensities, and transition IDs
            dictionaries for upper level.
        """
        self.snrs = {} # {(spectrum, lower wavenumber) : snr}
        self.intensities = {}  # {(spectrum, lower wavenumber) : intensity}
        self.unc={}
        self.transition_ids = {}  # {transition wavenumber : lower wavenumber}
        root_npts={}
        self.w_maxI={}            # {(spectrum): Wavenumber of transition with max intensity in a spectrum}
        calunc_per_1000={}

        # Get all lines in files containing this level and cross reference against known linelist for wavenumbers.

        energy_levels_file = set_parameters()[6]

        # Check level file for given lower level key.
        levels = open(energy_levels_file, 'r').readlines()

        for spectrum in self.all_spectrum_files:

            # JJL
            # transitions = popen(f'grep -e "- {self.upper_level}" {spectrum}', 'r')

            # for transition in transitions.readlines():
            with open(spectrum) as f:
                params = f.readline().split()
            resoln=float(params[0])  # Used for no. points/fwhm
            lowE=float(params[1])    # Used to get range for 7% unc.
            hiE=float(params[2])     
                                     # 7% from +- 5% discrepancy in D2 measurements
                                     # converts to unc of 2x5/sqrt(3) = 5.8%
                                     # added onto lamp unc. of 3.5% = 6.8% - round to 7 %
            calunc_per_1000[spectrum] = 70/(hiE-lowE) # Calibration uncertainty per 1000 cm-1.
                                              # Estimated from 7%/(sig_hi - sig_low)

            transitions = grep_open("- "+self.upper_level, spectrum)

            self.w_maxI[spectrum] = 0
            maxI = 0
            
            for transition in transitions:

                transition = transition.split()
                lower_level = transition[6]
                lower_level = lower_level.strip('*')
                self.snrs[spectrum, lower_level] = eval(transition[0])
                self.intensities[spectrum, lower_level] = eval(transition[1])
                if float(transition[1])> maxI :        # Find the wavenumber of the strongest line
                    maxI = float(transition[1])        # Used for calibration uncertainty
                    self.w_maxI[spectrum] = float(transition[4])
                    
                root_npts[spectrum,lower_level]=(0.001*float(transition[2])/resoln)    # no. points/FWHM.
                self.unc[spectrum,lower_level]=2.25/(self.snrs[spectrum,lower_level]**2 *root_npts[spectrum,lower_level])
                                                                             # Variance
                for level in levels:
                    level = level.split()

                    if (level[6] == lower_level):
                        lower_level_key = level[6].strip('*')
#                print(lower_level, lower_level_key,spectrum,self.snrs[spectrum,lower_level])

                for line in self.known_lines:
                    line = line.split()
#                    print(line)
                    if (line[6].strip('*') == lower_level_key):
                        self.transition_ids[eval(line[4])] = lower_level  
        self.wavenumbers = sorted(self.transition_ids.keys(), reverse=True)
   

        for spectrum,lower_level in self.unc :
#            print(spectrum,lower_level,self.transition_ids)
            wnum =  (list(self.transition_ids.keys())[list(self.transition_ids.values()).index(lower_level)])
            calunc = calunc_per_1000[spectrum]*(wnum-self.w_maxI[spectrum])/1000
            self.unc[spectrum,lower_level] = math.sqrt(calunc*calunc+self.unc[spectrum,lower_level])
#            print(spectrum,wnum,self.w_maxI[spectrum],calunc,self.unc[spectrum,lower_level])


    def show_values(self):
        """
        Description:
            Outputs SNRs and Intensities for upper level in each .I file
        """
        file_num = len(self.all_spectrum_files)
        width = 27+18*file_num

        self.log(f'SNRs and intensities of observed transitions from '
                 f'upper level {self.upper_level}:')
        self.log(''.join(['-']*width))

        file_name_line = f'| File Name:              | '
        for spectrum in self.all_spectrum_files: 
        	file_name_line += f'{spectrum[0:13]:>14s}  | '

        self.log(file_name_line)
        self.log(''.join(['-']*width))
        self.log(f'| Wavenumber | L Level    | '+f'SNR | Intensity | '*file_num)
        self.log(''.join(['-']*width))

        for wavenumber in self.wavenumbers:

            data_line = f"| {wavenumber:>10.3f} | {self.transition_ids[wavenumber]:>10s} | "

            for spectrum in self.all_spectrum_files:

                try:
                    #print(self.snrs[spectrum, self.transition_ids[wavenumber]], self.intensities[spectrum, self.transition_ids[wavenumber]])
                    data_line += (f"{round(self.snrs[spectrum, self.transition_ids[wavenumber]]):>3d} | "
                                  f"{round(self.intensities[spectrum, self.transition_ids[wavenumber]]):9d} | ")

                except: data_line += f"--- | --------- | "

            self.log(data_line)

        self.log(''.join(['-']*width)+'\n')

    def get_reference_level(self, event):
        """
        Description:
            Sets self.reference_level to entry contents after pressing Enter key (event).
            Normalizes the Intensities of each line with respect to the reference level,
            with its intensity set to 1000, then displays the wavenumbers of transitions,
            lower levels, SNRs, Intensities, LS configurations, and notes for observed lines.
        """
        self.reference_level = self.reference_level_entry.get()
        event.widget["foreground"] = "Black"
        self.normalize_spectrum()
        self.log(f"Normalizing each spectrum with respect to level {self.reference_level} = 1000.\n")
        self.show_values()

    def normalize_spectrum(self):
        """
        Normalizes each emission spectrum file's intensities, with self.reference_level
        having an intensity of 1000.
        """
        for spectrum in self.all_spectrum_files:

            try:

                normal = self.intensities[spectrum, self.reference_level]

            except:

                self.log(f"Spectrum {spectrum} does not contain level {self.reference_level}\n")
                continue

            for lower_level in self.transition_ids.values():

                try: self.intensities[spectrum, lower_level] = 1000*self.intensities[spectrum, lower_level]/normal
                except KeyError: pass

    def get_reference_file(self, event):
        """
        Description:
            Sets self.reference_file to entry contents after pressing Enter key (event).
        """
        self.reference_file = self.reference_file_entry.get()
        event.widget["foreground"] = "Black"
        self.log(f"Reference spectrum file selected: {self.reference_file}\n")

    def get_normal_level(self, event):

        self.normal_level = self.normal_level_entry.get()
        event.widget["foreground"] = "Black"
        self.log(f"Normalized levels using level {self.normal_level} from spectrum {self.reference_file}")
        comment = simpledialog.askstring(title="Comment",prompt="Why did you renormalize this spectrum?")
        self.log(comment)
        self.normalize_all_spectra()
        self.show_values()

    def normalize_all_spectra(self):
        """
        Description:
            Finds the weighted average intensity for the normal level for all files
            containing it, and uses this to renormalize all other intensities.
        """
        weighted_int = 0
        sum_weight = 0
        sqsum_int = 0

        for spectrum in self.all_spectrum_files:

            if (((spectrum,self.normal_level) in self.intensities.keys()) and ((spectrum,self.reference_level)) in self.intensities.keys()):
                
                # Intensity weighted
                weight = 1/(self.unc[spectrum,self.normal_level])**2
                if weight > 555:           # 555 = 1/(0.06**2/2)
                                           # Cap on unc. from calibration, as line ratio.
                                           # So all lines with SNR>24 get same weight.
                                           # Omit lamp unc. here as it's the same for all spectra.
                    weight = 555
                weighted_int += self.intensities[spectrum, self.normal_level] * weight
                sum_weight += weight
                sqsum_int += self.intensities[spectrum, self.normal_level]**2

        # Take weighted average/weight as effective SNR for transfer level
        # Renormalize SNRs by 1/sqrt(1/SNR(meas)^2 + 1/SNR(calc)^2)
        # -- this is root sum of squares of uncertainties

        avg_int = weighted_int / sum_weight
        sqsum_int = sqsum_int**(0.5) / sum_weight

        unc_renorm = (1/(self.unc[self.reference_file,self.normal_level]**2) + sum_weight)**(-0.5)
        int_renorm = avg_int / self.intensities[self.reference_file, self.normal_level]
        
        for lower_level in self.transition_ids.values():

            try:

                self.intensities[self.reference_file, lower_level] *= int_renorm

                # If level is transfer level, then SNR = unc_renorm
                if lower_level == self.normal_level:

                    self.unc[self.reference_file, lower_level] = unc_renorm
                

                # Otherwise add uncertainties to unc. of ref. level
                else:

                    self.unc[self.reference_file, lower_level] = (self.unc[spectrum, lower_level]**(2) + unc_renorm**(2))**(0.5)

            except Exception as e: pass
        

    def Delfil(self,event):        # Put back ability to delete lines (GN, July 21)

        self.delfile = self.Dfilenentry.get()
        event.widget["foreground"] = "Black"
        self.log(f"Deleting value from file {self.delfile}")
        
    def DelLine(self,event):       # Put back ability to delete lines (GN, July 21)
        self.delval =self.Delentry.get()
        del self.snrs[self.delfile,self.delval]
        del self.intensities[self.delfile,self.delval]
        self.log(f'{self.delfile}, {self.delval} deleted from list' )        
        comment = simpledialog.askstring(title="Comment",prompt="Why did you delete it? ")
        self.log(comment)
        self.show_values()

    def get_id_lines(self, event):

        self.id_lines_file = self.id_lines_entry.get()
        event.widget["foreground"] = "Black"
        self.log(f"ID'd lines file selected: {self.id_lines_file}\n")

    def display(self):
        """
        Description:
            Displays the wavenumbers of transitions, lower levels, SNRs,
            Intensities, LS configurations, and notes for observed lines.
        """
        total_int = 0
        level_keys = self.transition_ids.values()
        avg_int = {}
        sum_weight = {}
        sqsum = {}
        variance = {}
        stdev = {}
        waveno = {}
        upper_level = {}
        E1 = {}
        jupp = {}

        BFsq = 0

        self.log('Residuals:')

        for level_key in level_keys:

            avg_int[level_key] = 0
            sum_weight[level_key] = 0
            sqsum[level_key] = 0

            for spectrum in self.all_spectrum_files:

                if (spectrum, level_key) in self.intensities.keys():
                    weight = 1/(self.unc[spectrum,level_key])**2
                    if weight > 555:           # 555 = 1/(0.06**2/2)
                                               # Cap on unc. from calibration, as line ratio.
                                               # So all lines with SNR>24 get same weight.
                                               # Omit lamp unc. here as it's the same for all spectra.
                        weight = 555
                        
                    avg_int[level_key] += self.intensities[spectrum, level_key] * weight
                    sum_weight[level_key]  += weight
                    sqsum[level_key]   += self.intensities[spectrum, level_key]**(2)

            avg_int[level_key] /= sum_weight[level_key]

            # This line is fractional weighted standard deviation of intensities (u(I)/I in Sikstrom eqn 6)
            # Calibration uncertainty has already been taken into account in calculation of the individual uncs.
            stdev[level_key] = 1/math.sqrt(sum_weight[level_key])

            total_int += avg_int[level_key]
        
      
#
# Calculate the residual
#  
        get_calculations(E1)
        sumA=0
        delval=[]
        theoval={}
        for (x,y) in E1.keys():
            if x==self.upper_level:
                for wno in self.wavenumbers:
                    if -0.2  < wno - float(y) < 0.2:
                        delval.append((x,y))  
                        theoval[x,wno]=E1[x,y]
        for (x,y) in E1.keys():            
            if x==self.upper_level and (x,y) not in delval:
                sumA += E1[x,y]
                self.log(f'{y:>10s} | {E1[x,y]:>7.4f}')
#                print(x, y, E1[x,y])      # Do we want all these in the output?
        frac_resid = sumA*self.lifetimes[self.upper_level]/1000.
        total_int *= (1+frac_resid)
#        print (theoval)
       
#
# Print the header for the results
#
        self.log(''.join(['-']*108)+'\n')
        self.log('Summary:')
        self.log(f'Level = {self.upper_level:>10s}, Lifetime = {self.lifetimes[self.upper_level]} ns, residual =, {frac_resid*100.:>6.3f} %')
        self.log(''.join(['-']*108))
        self.log(f'| Wavenumber | L Level    | Weight | B. Fraction | Uncertainty | Trans. Prob. | Uncertainty |  Theoretical |')
        self.log(''.join(['-']*108))

        for level_key in level_keys:            # This is 2nd term in eqn 7 of Sikstrom. 
            BFsq += (avg_int[level_key]/total_int)**2 * stdev[level_key]**2

# Correction for the residual (50 % uncertainty estimate on residual)
        BFsq += frac_resid**2 * 0.25

        for wavenumber in self.wavenumbers:

            level_key = self.transition_ids[wavenumber]
                                                       # First two terms of eqn 7 of Sikstrom (BF uncertainty)
            variance[level_key] = stdev[level_key]**2 * (1 - 2*avg_int[level_key]/total_int)+ BFsq
#            print('unc: ',math.sqrt(variance[level_key]))
            avg_int[level_key] /= total_int
            aval = (1000 * avg_int[level_key]) / (self.lifetimes[self.upper_level])
            unc_aval = 100 * ((variance[level_key]) + self.life_unc**(2))**(0.5)
#            print(self.upper_level, wavenumber)
            try: 
                theo_val = theoval[self.upper_level,wavenumber]
            except:
                theo_val = 0

            self.log(f'| {wavenumber:>10.3f} | {level_key:>10s} | {round(sum_weight[level_key]):>6d} | '
                     f'{100*avg_int[level_key]:>10.3f}% | {100*math.sqrt(variance[level_key]):>10.1f}% | '
                     f'{aval:>12.3f} | {unc_aval:>10.0f}% | {theo_val:>12.3f} |')


        self.log(''.join(['-']*108)+'\n')  
        
        comment = simpledialog.askstring(title="Comment",prompt="Comments on results")
        self.log(comment)




 #   def save(self):
 #       """
 #       Description:
 #           Logs the upper level, reference level, reference file, normalizing level,
 #           and normalizing intensity value (which is set to 1000.).
 #       """
 #       self.log(f"Upper level: {self.upper_level}")
 #       self.log(f"Reference level: {self.reference_level}")
 #       self.log(f"Reference file: {self.reference_file}")
 #       self.log(f"Normalized level: {self.normal_level}")
 #       self.log(f"Normalizing value: {1000.}\n")

    def make_widgets(self):
        """
        Description:
            Creates the widgets the user interacts with.
        """

        self.display_frame = Frame(self)
        self.display_frame.pack(side=RIGHT, fill=BOTH, expand=YES)
        self.display_scroll_y = Scrollbar(self.display_frame, orient=VERTICAL)
        self.display_scroll_y.pack(side=RIGHT, fill=Y)
        self.display_text = Text(self.display_frame, wrap=NONE, width=0)
        self.display_text.pack(side=TOP, fill=BOTH, expand=YES)
        self.display_scroll_y.config(command=self.display_text.yview)
        self.display_scroll_x = Scrollbar(self.display_frame, orient=HORIZONTAL, command=self.display_text.xview)
        self.display_scroll_x.pack(side=TOP, fill=X)

        self.display_text.config(xscrollcommand=self.display_scroll_x.set,
                                 yscrollcommand=self.display_scroll_y.set)


        self.widgets_frame = Frame(self)
        self.widgets_frame.pack(side=LEFT, fill=BOTH)
        self.toolbar_frame = Frame(self.widgets_frame, height=30)
        self.toolbar_frame.pack(side=TOP, fill=X)
        Button(self.toolbar_frame, text='Quit', command=self.quit).pack(side=LEFT, fill=X, expand=YES)
#        Button(self.toolbar_frame, text='Save', command=self.save).pack(side=LEFT, fill=X, expand=YES)
        Button(self.toolbar_frame, text='Results', command=self.display).pack(side=LEFT, fill=X, expand=YES)

        self.label_frame = Frame(self.widgets_frame)
        self.label_frame.pack(side=LEFT, fill=BOTH)
        self.entry_frame = Frame(self.widgets_frame)
        self.entry_frame.pack(side=RIGHT, fill=BOTH)

        Label(self.label_frame, text="Save As:").pack(side=TOP, fill=X, pady=2.5)
        self.file_name_entry = Entry(self.entry_frame, foreground="Blue")
        self.file_name_entry.pack(side=TOP, fill=X, pady=2)
        self.file_name_entry.insert(0, set_parameters()[0])
        self.file_name_entry.bind("<Key-Return>", self.set_output)

        Label(self.label_frame, text="ID'd lines file:").pack(side=TOP, fill=X, pady=2.5)
        self.id_lines_entry = Entry(self.entry_frame, foreground="Blue")
        self.id_lines_entry.pack(side=TOP, fill=X, pady=2)
        self.id_lines_entry.insert(0, set_parameters()[5])
        self.id_lines_entry.bind("<Key-Return>", self.get_id_lines)

        Label(self.label_frame, text="Upper Level:").pack(side=TOP, fill=X, pady=2.5)
        self.upper_level_entry = Entry(self.entry_frame, foreground="Blue")
        self.upper_level_entry.pack(side=TOP, fill=X, pady=2)
        self.upper_level_entry.insert(0, set_parameters()[1])
        self.upper_level_entry.bind("<Key-Return>", self.get_upper_level)

        Label(self.label_frame, text="Reference Level:").pack(side=TOP, fill=X, pady=2.5)
        self.reference_level_entry = Entry(self.entry_frame, foreground="Blue")
        self.reference_level_entry.pack(side=TOP, fill=X, pady=2)
        self.reference_level_entry.insert(0, set_parameters()[2])
        self.reference_level_entry.bind("<Key-Return>", self.get_reference_level)

        Label(self.label_frame, text="Reference File:").pack(side=TOP, fill=X, pady=2.5)
        self.reference_file_entry = Entry(self.entry_frame, foreground="Blue")
        self.reference_file_entry.pack(side=TOP, fill=X, pady=2)
        self.reference_file_entry.insert(0, set_parameters()[3])
        self.reference_file_entry.bind("<Key-Return>", self.get_reference_file)

        Label(self.label_frame, text="Normalized Level:").pack(side=TOP, fill=X, pady=2.5)
        self.normal_level_entry = Entry(self.entry_frame, foreground="Blue")
        self.normal_level_entry.pack(side=TOP, fill=X, pady=2)
        self.normal_level_entry.insert(0, set_parameters()[4])
        self.normal_level_entry.bind("<Key-Return>", self.get_normal_level)

# Put back ability to delete lines (GN, July21)        
        Label(self.label_frame, text="Delete file:").pack(side=TOP, fill=X, pady=2.5)
        self.Dfilenentry = Entry(self.entry_frame, foreground="Blue")
        self.Dfilenentry.pack(side=TOP, fill=X, pady=2)
        self.Dfilenentry.insert(0, set_parameters()[5])
        self.Dfilenentry.bind("<Key-Return>", self.Delfil)
        
        Label(self.label_frame, text="Delete level:").pack(side=TOP, fill=X, pady=2.5)
        self.Delentry = Entry(self.entry_frame, foreground="Blue")
        self.Delentry.pack(side=TOP, fill=X, pady=2)
        self.Delentry.insert(0, set_parameters()[6])
        self.Delentry.bind("<Key-Return>", self.DelLine)

     
if __name__ == '__main__':

    BranchingFractionCalc().mainloop()
