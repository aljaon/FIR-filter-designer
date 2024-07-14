import numpy as np
import matplotlib.pyplot as plt
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
import tkinter as tk
from tkinter import HORIZONTAL

max_taps = 1000
max_sample_rate = 1000000
max_cutoff = 1000000
max_interger_precision = 1000000

class tk_scale_variable:
    def __init__(self, frame, name:str, min, max):
        
        #Name and range
        self.name = name
        self.max = max
        self.min = min
        
        #Frame
        self.frame = tk.Frame(frame)
        self.frame.pack()
        
        #Variable
        self.tk_variable = tk.IntVar()
        self.tk_variable.set(min)
        self.tk_variable_str = tk.StringVar()
        self.tk_variable_str.set(f"{name} {self.tk_variable.get()}")
        
        self.tk_variable_label = ttk.Label(self.frame, bootstyle="default", textvariable = self.tk_variable_str)
        self.tk_variable_slider = ttk.Scale(self.frame, from_= min, to = max, bootstyle=DARK, command=self.update_tk_variable, variable = self.tk_variable)
        
        self.tk_variable_increment_button = ttk.Button(self.frame, text="+1", command=self.tk_variable_increment)
        self.tk_variable_increment_button_10 = ttk.Button(self.frame, text="+10", command=self.tk_variable_increment_10)
        self.tk_variable_increment_button_100 = ttk.Button(self.frame, text="+100", command=self.tk_variable_increment_100)
        self.tk_variable_increment_button_1000 = ttk.Button(self.frame, text="+1000", command=self.tk_variable_increment_1000)
        
        self.tk_variable_decrement_button = ttk.Button(self.frame, text="-1", command=self.tk_variable_decrement)
        self.tk_variable_decrement_button_10 = ttk.Button(self.frame, text="-10", command=self.tk_variable_decrement_10)
        self.tk_variable_decrement_button_100 = ttk.Button(self.frame, text="-100", command=self.tk_variable_decrement_100)
        self.tk_variable_decrement_button_1000 = ttk.Button(self.frame, text="-1000", command=self.tk_variable_decrement_1000)    
    
    def display(self):
        
        self.frame.pack(side=TOP, padx=50, pady=5)
        
        self.tk_variable_label.pack(side=TOP, padx=50, pady=5)
        self.tk_variable_slider.pack(side=BOTTOM, padx=50, pady=5)
        
        self.tk_variable_decrement_button_1000.pack(side=LEFT, padx=5, pady=5)
        self.tk_variable_decrement_button_100.pack(side=LEFT, padx=5, pady=5)
        self.tk_variable_decrement_button_10.pack(side=LEFT, padx=5, pady=5)
        self.tk_variable_decrement_button.pack(side=LEFT, padx=5, pady=5)
        self.tk_variable_increment_button_1000.pack(side=RIGHT, padx=5, pady=5)
        self.tk_variable_increment_button_100.pack(side=RIGHT, padx=5, pady=5)
        self.tk_variable_increment_button_10.pack(side=RIGHT, padx=5, pady=5)
        self.tk_variable_increment_button.pack(side=RIGHT, padx=5, pady=5)
    
    def undisplay(self):
        self.frame.pack_forget()
 
    def update_tk_variable(self, *args):
        self.tk_variable_str.set(f"{self.name} {self.tk_variable.get()}")    
    
    def tk_variable_increment(self):
        self.tk_variable.set(self.tk_variable.get() + 1)
        if self.tk_variable.get() > self.max:
                self.tk_variable.set(self.max)
        self.update_tk_variable(self)
    def tk_variable_increment_10(self):
        self. tk_variable.set(self. tk_variable.get() + 10)
        if self. tk_variable.get() > self.max:
            self. tk_variable.set(self.max)
        self.update_tk_variable(self)
    def tk_variable_increment_100(self):
        self. tk_variable.set(self. tk_variable.get() + 100)
        if self. tk_variable.get() > self.max:
            self. tk_variable.set(self.max)
        self.update_tk_variable(self)
    def tk_variable_increment_1000(self):
        self. tk_variable.set(self. tk_variable.get() + 1000)
        if self.tk_variable.get() > self.max:
            self.tk_variable.set(self.max)
        self.update_tk_variable(self)
            
    def tk_variable_decrement(self):
        self. tk_variable.set(self. tk_variable.get() - 1)
        if self. tk_variable.get() < 1:
            self. tk_variable.set(1)
        self.update_tk_variable(self)
    def tk_variable_decrement_10(self):
        self. tk_variable.set(self. tk_variable.get() - 10)
        if self. tk_variable.get() < 1:
            self. tk_variable.set(1)
        self.update_tk_variable(self)
    def tk_variable_decrement_100(self):
        self. tk_variable.set(self. tk_variable.get() - 100)
        if self. tk_variable.get() < 1:
            self. tk_variable.set(1)
        self.update_tk_variable(self)
    def tk_variable_decrement_1000(self):
        self. tk_variable.set(self.tk_variable.get() - 1000)
        if self. tk_variable.get() < 1:
            self. tk_variable.set(1)
        self.update_tk_variable(self)
                             
class filter:
    
    def __init__(self, root):
        
        #root
        self.root = root
        
        #Frame
        self.interger_precision = tk_scale_variable(root, "Interger Precision (0 = Use full precision floating Point)", 0, max_interger_precision)
        self.sample_rate = tk_scale_variable(root, "Sample Rate", 1, max_sample_rate)
        self.number_of_taps = tk_scale_variable(root, "N taps", 1, max_taps)

        #Window Function
        self.window_function = ""
           
    def display_base(self):
        
        self.interger_precision.display()
        self.sample_rate.display()
        self.number_of_taps.display()
    
    def undisplay_base(self):
        
        self.interger_precision.undisplay()
        self.sample_rate.undisplay()
        self.number_of_taps.undisplay()
        
class low_pass_filter(filter):
    
    def __init__(self, root):
        
        filter.__init__(self, root)
        
        #Low pass
        self.cutoff = tk_scale_variable(root, "Cutoff frequency fL [Hz]:", 1, max_cutoff)
              
    def display(self):
        
        filter.display_base(self)
        self.cutoff.display()
        
    def undisplay(self):
        filter.undisplay_base(self)
        self.cutoff.undisplay()
              
    def run(self):
        
        # Sampling rate.
        fS = self.sample_rate.tk_variable.get()  
        
        # Cutoff frequency.
        fL = self.cutoff.tk_variable.get()
        
        # Filter length (Odd number)
        N = self.number_of_taps.tk_variable.get()  
        if N % 2 == 0:
            N += 1
            
        # Convert coeffecients to interger 
        Precision = self.interger_precision.tk_variable.get()
        
        # Compute sinc filter.
        self.h = np.sinc(2 * fL / fS * (np.arange(N) - (N - 1) / 2))

        # Apply window.
        if self.window_function == "Rectangular":
            pass
        elif self.window_function == "Blackman":
            self.h *= np.blackman(N)
        elif self.window_function == "Hamming":
            self.h *= np.hamming(N)
             
        # Normalize to get unity gain.
        self.h /= np.sum(self.h)
        
        #Set precision for interger coeffecients
        if Precision > 0:
            self.h *= Precision
            self.h = self.h.astype(int)
        
        #View coeffecietns
        print(self.h)
        
        #Create Impulse singal
        signal = np.zeros(N)
        signal[-1] = 1

        #Apply Filter + Get FFT Magnitude (normalise again if Presision is set)
        s = np.convolve(signal, self.h)
        if Precision > 0:
            s /= Precision
        s_fft = abs(np.fft.fft(s))
        
        #Plot all
        x1 = np.arange(0, N * 2 - 1)
        x2 = np.arange(1, N) * fS/(2 * N)
        
        fig, axs = plt.subplots(3, 1, figsize=(10, 15))
        fig.canvas.manager.set_window_title(f'Low Pass fS = {fS}, fL = {fL}, N = {N}, Precision = {Precision}')
        
        axs[0].plot(x1, s)
        axs[0].set_title("Impulse response")
        axs[0].set_xlabel('t(s)', fontsize=10)
        axs[0].set_ylabel('A', fontsize=10)
        
        axs[1].plot(x2, s_fft[1:N])
        axs[1].set_title("FFT (G)")
        axs[1].set_xlabel('f(Hz)', fontsize=10)
        axs[1].set_ylabel('A', fontsize=10)
        
        axs[2].plot(x2, 10 * np.log(s_fft[1:N]))
        axs[2].set_title("FFT (dB)")
        axs[2].set_xlabel('f(Hz)', fontsize=10)
        axs[2].set_ylabel('dB', fontsize=10)
        
        # Show the plots
        plt.subplots_adjust(hspace=0.5)
        plt.show()

class high_pass_filter(filter):
    
    def __init__(self, root):
        
        filter.__init__(self, root)
        
        #high pass cutoff
        self.cutoff = tk_scale_variable(root, "Cutoff frequency fH [Hz]:", 1, max_cutoff)
              
    def display(self):
        
        filter.display_base(self)
        self.cutoff.display()
        
    def undisplay(self):
        filter.undisplay_base(self)
        self.cutoff.undisplay()
              
    def run(self):
        
        # Sampling rate.
        fS = self.sample_rate.tk_variable.get()  
        
        # Cutoff frequency.
        fH = self.cutoff.tk_variable.get()
        
        # Filter length (Odd number)
        N = self.number_of_taps.tk_variable.get()  
        if N % 2 == 0:
            N += 1
            
        # Convert coeffecients to interger 
        Precision = self.interger_precision.tk_variable.get()
        
        # Compute sinc filter.
        self.h = np.sinc(2 * fH / fS * (np.arange(N) - (N - 1) / 2))

        # Apply window.
        if self.window_function == "Rectangular":
            pass
        elif self.window_function == "Blackman":
            self.h *= np.blackman(N)
        elif self.window_function == "Hamming":
            self.h *= np.hamming(N)
             
        # Normalize to get unity gain.
        self.h /= np.sum(self.h)
        
        #Create highpass from lowpass
        self.h = -self.h
        self.h[(N - 1) // 2] += 1
        
        #Set precision for interger coeffecients
        if Precision > 0:
            self.h *= Precision
            self.h = self.h.astype(int)
        
        #View coeffecietns
        print(self.h)
        
        #Create Impulse singal
        signal = np.zeros(N)
        signal[-1] = 1

        #Apply Filter + Get FFT Magnitude (normalise again if Presision is set)
        s = np.convolve(signal, self.h)
        if Precision > 0:
            s /= Precision
        s_fft = abs(np.fft.fft(s))
        
        #Plot all
        x1 = np.arange(0, N * 2 - 1)
        x2 = np.arange(1, N) * fS/(2 * N)
        
        fig, axs = plt.subplots(3, 1, figsize=(10, 15))
        fig.canvas.manager.set_window_title(f'High Pass fS = {fS}, fL = {fH}, N = {N}, Precision = {Precision}')
        
        axs[0].plot(x1, s)
        axs[0].set_title("Impulse response")
        axs[0].set_xlabel('t(s)', fontsize=10)
        axs[0].set_ylabel('A', fontsize=10)
        
        axs[1].plot(x2, s_fft[1:N])
        axs[1].set_title("FFT (G)")
        axs[1].set_xlabel('f(Hz)', fontsize=10)
        axs[1].set_ylabel('A', fontsize=10)
        
        axs[2].plot(x2, 10 * np.log(s_fft[1:N]))
        axs[2].set_title("FFT (dB)")
        axs[2].set_xlabel('f(Hz)', fontsize=10)
        axs[2].set_ylabel('dB', fontsize=10)
        
        # Show the plots
        plt.subplots_adjust(hspace=0.5)
        plt.show()
               
class band_pass_filter(filter):
    
    def __init__(self, root):
        
        filter.__init__(self, root)
        
        #Frames
        self.cutoff_frame = tk.Frame(self.root)
        self.cutoff_frame.pack()
        
        #high pass and low pass cutoff
        self.cutoff_low = tk_scale_variable(root, "Cutoff frequency fL [Hz]:", 1, max_cutoff)
        self.cutoff_high = tk_scale_variable(root, "Cutoff frequency fH [Hz]:", 1, max_cutoff)  
        
    def display(self):
        
        filter.display_base(self)
        self.cutoff_low.display()
        self.cutoff_high.display()
        
    def undisplay(self):
        filter.undisplay_base(self)
        self.cutoff_low.undisplay()
        self.cutoff_high.undisplay()
              
    def run(self):
        
        # Sampling rate.
        fS = self.sample_rate.tk_variable.get()  
        
        # Cutoff frequencies.
        fL = self.cutoff_low.tk_variable.get()
        fH = self.cutoff_high.tk_variable.get()
        
        # Filter length (Must be off)
        N = self.number_of_taps.tk_variable.get()
        if N % 2 == 0:
            N += 1
        half_N = int((N + 1)/2)
        
        # Convert coeffecients to interger 
        Precision = self.interger_precision.tk_variable.get()
        
        # High pass filter with fL, Low pass filter with fH
        self.lp = np.sinc(2 * fH / fS * (np.arange(half_N) - (half_N - 1) / 2))
        self.hp = np.sinc(2 * fL / fS * (np.arange(half_N) - (half_N - 1) / 2))
        
        # Apply window.
        if self.window_function == "Rectangular":
            pass
        elif self.window_function == "Blackman":
            self.lp *= np.blackman(half_N)
            self.hp *= np.blackman(half_N)
        elif self.window_function == "Hamming":
            self.lp *= np.hamming(half_N)
            self.hp *= np.hamming(half_N)
             
        # Normalize to get unity gain.
        self.lp /= np.sum(self.lp)
        self.hp /= np.sum(self.hp)
            
        #Create highpass from lowpass
        self.hp = -self.hp
        self.hp[(half_N - 1) // 2] += 1
        
        #Convolve together to get complete filter
        self.h = np.convolve(self.lp, self.hp)
        
        #Set precision for interger coeffecients
        if Precision > 0:
            self.h *= Precision
            self.h = self.h.astype(int)
            
        #View coeffecietns
        print(self.h)
        
        #Create Impulse singal
        signal = np.zeros(N)
        signal[-1] = 1

        #Apply Filter + Get FFT Magnitude (normalise again if Presision is set)
        s = np.convolve(signal, self.h)
        if Precision > 0:
            s /= Precision
        s_fft = abs(np.fft.fft(s))
        
        #Plot all
        x1 = np.arange(0, N * 2 - 1)
        x2 = np.arange(1, N) * fS/(2 * N)
        
        fig, axs = plt.subplots(3, 1, figsize=(10, 15))
        fig.canvas.manager.set_window_title(f'Band Pass fS = {fS}, fL = {fL}, fH = {fH}, N = {N}, Precision = {Precision}')
        
        axs[0].plot(x1, s)
        axs[0].set_title("Impulse response")
        axs[0].set_xlabel('t(s)', fontsize=10)
        axs[0].set_ylabel('A', fontsize=10)
        
        axs[1].plot(x2, s_fft[1:N])
        axs[1].set_title("FFT (G)")
        axs[1].set_xlabel('f(Hz)', fontsize=10)
        axs[1].set_ylabel('A', fontsize=10)
        
        axs[2].plot(x2, 10 * np.log(s_fft[1:N]))
        axs[2].set_title("FFT (dB)")
        axs[2].set_xlabel('f(Hz)', fontsize=10)
        axs[2].set_ylabel('dB', fontsize=10)
        
        # Show the plots
        plt.subplots_adjust(hspace=0.5)
        plt.show()

class band_stop_filter(filter):
    
    def __init__(self, root):
        
        filter.__init__(self, root)
        
        #Frames
        self.cutoff_frame = tk.Frame(self.root)
        self.cutoff_frame.pack()
        
        #high pass and low pass cutoff
        self.cutoff_low = tk_scale_variable(root, "Cutoff frequency fL [Hz]:", 1, max_cutoff)
        self.cutoff_high = tk_scale_variable(root, "Cutoff frequency fH [Hz]:", 1, max_cutoff)  
        
    def display(self):
        
        filter.display_base(self)
        self.cutoff_low.display()
        self.cutoff_high.display()
        
    def undisplay(self):
        filter.undisplay_base(self)
        self.cutoff_low.undisplay()
        self.cutoff_high.undisplay()
              
    def run(self):
        
        # Sampling rate.
        fS = self.sample_rate.tk_variable.get()  
        
        # Cutoff frequencies.
        fL = self.cutoff_low.tk_variable.get()
        fH = self.cutoff_high.tk_variable.get()
        
        # Filter length (Must be off)
        N = self.number_of_taps.tk_variable.get()
        if N % 2 == 0:
            N += 1
        
        # Convert coeffecients to interger 
        Precision = self.interger_precision.tk_variable.get()
        
        # High pass filter with fL, Low pass filter with fH
        self.lp = np.sinc(2 * fL / fS * (np.arange(N) - (N - 1) / 2))
        self.hp = np.sinc(2 * fH / fS * (np.arange(N) - (N - 1) / 2))
        
        # Apply window.
        if self.window_function == "Rectangular":
            pass
        elif self.window_function == "Blackman":
            self.lp *= np.blackman(N)
            self.hp *= np.blackman(N)
        elif self.window_function == "Hamming":
            self.lp *= np.hamming(N)
            self.hp *= np.hamming(N)
             
        # Normalize to get unity gain.
        self.lp /= np.sum(self.lp)
        self.hp /= np.sum(self.hp)
 
        #Create highpass from lowpass
        self.hp = -self.hp
        self.hp[(N - 1) // 2] += 1
        
        # Add both filters.
        self.h = self.hp + self.lp
        
        #Set precision for interger coeffecients
        if Precision > 0:
            self.h *= Precision
            self.h = self.h.astype(int)
            
        #View coeffecietns
        print(self.h)
        
        #Create Impulse singal
        signal = np.zeros(N)
        signal[-1] = 1

        #Apply Filter + Get FFT Magnitude (normalise again if Presision is set)
        s = np.convolve(signal, self.h)
        if Precision > 0:
            s /= Precision
        s_fft = abs(np.fft.fft(s))
        
        #Plot all
        x1 = np.arange(0, N * 2 - 1)
        x2 = np.arange(1, N) * fS/(2 * N)
        
        fig, axs = plt.subplots(3, 1, figsize=(10, 15))
        fig.canvas.manager.set_window_title(f'Band Stop fS = {fS}, fL = {fL}, fH = {fH}, N = {N}, Precision = {Precision}')
        
        axs[0].plot(x1, s)
        axs[0].set_title("Impulse response")
        axs[0].set_xlabel('t(s)', fontsize=10)
        axs[0].set_ylabel('A', fontsize=10)
        
        axs[1].plot(x2, s_fft[1:N])
        axs[1].set_title("FFT (G)")
        axs[1].set_xlabel('f(Hz)', fontsize=10)
        axs[1].set_ylabel('A', fontsize=10)
        
        axs[2].plot(x2, 10 * np.log(s_fft[1:N]))
        axs[2].set_title("FFT (dB)")
        axs[2].set_xlabel('f(Hz)', fontsize=10)
        axs[2].set_ylabel('dB', fontsize=10)
        
        # Show the plots
        plt.subplots_adjust(hspace=0.5)
        plt.show()
       
class window:
    def __init__(self):
        
        #Root
        self.root = tk.Tk()
        
        #Frames
        self.title_frame = tk.Frame(self.root)  
        self.filter_select_frame = tk.Frame(self.root)
        self.title_frame.pack()
        self.filter_select_frame.pack()  
        
        #Filters
        self.low_pass_filter = low_pass_filter(self.root)
        self.high_pass_filter = high_pass_filter(self.root)
        self.band_pass_filter = band_pass_filter(self.root)
        self.band_stop_filter = band_stop_filter(self.root)
        
        #Window Function Select
        self.window_function = tk.StringVar()
        self.window_function.set('Rectangular')
        
        #Filter Select (Default = low pass)
        self.filter_select_var = tk.IntVar()
        self.filter_select_var.set(1)

        
    def exit(self):
        exit()
        
    def filter_run(self): 
        if self.filter_select_var.get() == 1:
            self.low_pass_filter.window_function = self.window_function.get()
            self.low_pass_filter.run() 
        elif self.filter_select_var.get() == 2:
            self.high_pass_filter.window_function = self.window_function.get()
            self.high_pass_filter.run() 
        elif self.filter_select_var.get() == 3:
            self.band_pass_filter.window_function = self.window_function.get()
            self.band_pass_filter.run() 
        elif self.filter_select_var.get() == 4:
            self.band_stop_filter.window_function = self.window_function.get()
            self.band_stop_filter.run() 
            
    def select_low_pass(self):
        self.low_pass_filter.display()
        self.high_pass_filter.undisplay()
        self.band_pass_filter.undisplay()
        self.band_stop_filter.undisplay()
    def select_high_pass(self):
        self.high_pass_filter.display()
        self.low_pass_filter.undisplay()
        self.band_pass_filter.undisplay()
        self.band_stop_filter.undisplay()
    def select_band_pass(self):
        self.band_pass_filter.display()
        self.low_pass_filter.undisplay()
        self.high_pass_filter.undisplay()
        self.band_stop_filter.undisplay() 
    def select_band_stop(self):
        self.band_stop_filter.display()
        self.low_pass_filter.undisplay()
        self.high_pass_filter.undisplay()
        self.band_pass_filter.undisplay()
        
            
    def run(self):
        
        label = ttk.Label(self.title_frame, text="Select Window Type")
        label.pack(side=TOP, padx=5, pady=5)
        
        self.combobox = ttk.Combobox(self.title_frame, bootstyle="primary", width = 27, textvariable = self.window_function)
        self.combobox['values'] = ('Rectangular', 'Blackman', 'Hamming')
        self.combobox.pack(side=BOTTOM, padx=5, pady=5)
        
        ttk.Button(self.filter_select_frame, text='Run', bootstyle=PRIMARY, command = self.filter_run).pack(side=BOTTOM, padx=5, pady=5)
        ttk.Radiobutton(self.filter_select_frame, bootstyle="primary-toolbutton", text="Low Pass", variable=self.filter_select_var, value=1, command=self.select_low_pass).pack(side=LEFT, padx=5, pady=5,)
        ttk.Radiobutton(self.filter_select_frame, bootstyle="primary-toolbutton", text="High Pass", variable=self.filter_select_var, value=2, command=self.select_high_pass).pack(side=LEFT, padx=5, pady=5)
        ttk.Radiobutton(self.filter_select_frame, bootstyle="primary-toolbutton", text="Band Pass", variable=self.filter_select_var, value=3, command=self.select_band_pass).pack(side=LEFT, padx=5, pady=5)
        ttk.Radiobutton(self.filter_select_frame, bootstyle="primary-toolbutton", text="Band Stop", variable=self.filter_select_var, value=4, command=self.select_band_stop).pack(side=LEFT, padx=5, pady=5)
        self.select_low_pass()
        
        self.root.mainloop()
           
win = window()
win.run()
