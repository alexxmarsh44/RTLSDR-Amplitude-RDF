import matplotlib.pyplot as plt
import math
from matplotlib import animation
from scipy import signal as sc
import tkinter as Tk
import numpy as np
from scipy.signal import butter, lfilter, freqz
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.signal import kaiserord, lfilter, firwin, freqz
import sys
from rtlsdr import RtlSdr
import time
from threading import Thread
from multiprocessing import Process
from tkinter import filedialog as fd
from matplotlib.widgets import Slider, Button
np.set_printoptions(threshold=sys.maxsize)

serial_numbers = RtlSdr.get_device_serial_addresses()
print(serial_numbers)
#device_index = RtlSdr.get_device_index_by_serial('00000001')
sdr = RtlSdr(serial_number= '00000001')
sdr2 = RtlSdr(serial_number= '00000002')
sdr3 = RtlSdr(serial_number= '00000003')
sdr4 = RtlSdr(serial_number= '00000004')
# configure device
sdr.sample_rate = 3.048e6  # Hz   # PPM
sdr.gain = 45
sdr2.sample_rate = 3.048e6  # Hz  # PPM
sdr2.gain = 45
sdr3.sample_rate = 3.048e6  # Hz   # PPM
sdr3.gain = 45
sdr4.sample_rate = 3.048e6  # Hz # PPM
sdr4.gain = 45

plt.style.use('dark_background')
fig = plt.figure()
fig.set_size_inches(10, 6)
root = Tk.Tk()
root.winfo_toplevel().title("Radar")
label = Tk.Label(root, text="RTL Radar").grid(column=0, row=0)
img = Tk.PhotoImage(file='Screenshot_2022-11-01_21-04-52.ico')
#root.tk.call('wm', 'iconphoto', root._w, img)
width= root.winfo_screenwidth()
height= root.winfo_screenheight()
ax1 = fig.add_subplot(3, 1, 1)
ax3 = fig.add_subplot(3, 1, 3,  projection='polar')
ax2 = fig.add_subplot(3, 1, 2)
var1 = Tk.IntVar()

class RDF:
    def butter_lowpass(cutoff, fs, order=5):
        return butter(order, cutoff, fs=fs, btype='low', analog=False)

    def butter_lowpass_filter(data, cutoff, fs, order=5):
        b, a = butter_lowpass(cutoff, fs, order=order)
        y = lfilter(b, a, data)
        return y

    def normalize(arr, t_min, t_max):
        norm_arr = []
        diff = t_max - t_min
        diff_arr = max(arr) - min(arr)
        for i in arr:
            temp = (((i - min(arr)) * diff) / diff_arr) + t_min
            norm_arr.append(temp)
        return norm_arr

    i1 = 1

    def roto():
        global i1
        step = 4
        i1 += step
        if i1 >= 360:
            i1 = 1
        else:
            pass
        return np.deg2rad(i1)

    def antsel(ant):
        if ant == None:
            degree = 45
            return degree
        if ant == 0:
            degree = 45
            return degree
        if ant == 1:
            degree = 135
            return degree
        if ant == 2:
            degree = 225
            return degree
        if ant == 3:
            degree = 315
            return degree
        else:
            degree = 45
            return degree

    def samples22():

        if var1.get() == 0:
            samples = sdr.read_samples(1024)
            samples2 = sdr2.read_samples(1024)
            samples3 = sdr3.read_samples(1024)
            samples4 = sdr4.read_samples(1024)
            return samples, samples2, samples3, samples4
        if var1.get() == 1:
            samples = sdr.read_samples(1024)
            samples2 = sdr2.read_samples(1024)
            samples3 = sdr3.read_samples(1024)
            samples4 = sdr4.read_samples(1024)
            samples5 = butter_lowpass_filter(samples, 0.03e6, 2.4e6, 10)
            samples6 = butter_lowpass_filter(samples2, .03e6, 2.4e6, 10)
            samples7 = butter_lowpass_filter(samples3, 0.03e6, 2.4e6, 10)
            samples8 = butter_lowpass_filter(samples4, 0.03e6, 2.4e6, 10)
            return samples5, samples6, samples7, samples8
file1 = open("Record", "w")
file2 = open("Record2", "w")
file3 = open("Record3", "w")
file4 = open("Record4", "w")

var2 = Tk.IntVar()

class GUI:
    def freqdown1():
        if sdr.center_freq > 0.1e6:
            sdr.center_freq = sdr.center_freq - .01e6
            sdr2.center_freq = sdr2.center_freq - .01e6
            sdr3.center_freq = sdr3.center_freq - .01e6
            sdr4.center_freq = sdr4.center_freq - .01e6
        else:
            pass
        return print(sdr.center_freq)

    def freqdown2():
        if sdr.center_freq > 24e6:
            sdr.center_freq = sdr.center_freq - 25e6
            sdr2.center_freq = sdr2.center_freq - 25e6
            sdr3.center_freq = sdr3.center_freq - 25e6
            sdr4.center_freq = sdr4.center_freq - 25e6
        else:
            pass
        return print(sdr.center_freq)

    def frequp1():
        if sdr.center_freq < 1700e6:
            sdr.center_freq = sdr.center_freq + .01e6
            sdr2.center_freq = sdr2.center_freq + .01e6
            sdr3.center_freq = sdr3.center_freq + .01e6
            sdr4.center_freq = sdr4.center_freq + .01e6
        else:
            pass
        return print(sdr.center_freq)

    def frequp2():
        if sdr.center_freq < 1700e6:
            sdr.center_freq = sdr.center_freq + 25e6
            sdr2.center_freq = sdr2.center_freq + 25e6
            sdr3.center_freq = sdr3.center_freq + 25e6
            sdr4.center_freq = sdr4.center_freq + 25e6
        else:
            pass
        return print(sdr.center_freq)

    def freqdown3():
        if sdr.center_freq > 1e6:
            sdr.center_freq = sdr.center_freq - .1e6
            sdr2.center_freq = sdr2.center_freq - .1e6
            sdr3.center_freq = sdr3.center_freq - .1e6
            sdr4.center_freq = sdr4.center_freq - .1e6
        else:
            pass
        return print(sdr.center_freq)

    def frequp3():
        if sdr.center_freq < 1700e6:
            sdr.center_freq = sdr.center_freq + .1e6
            sdr2.center_freq = sdr2.center_freq + .1e6
            sdr3.center_freq = sdr3.center_freq + .1e6
            sdr4.center_freq = sdr4.center_freq + .1e6
        else:
            pass
        return print(sdr.center_freq)

    def slider_changed(event):

        sdr.center_freq = slider.get()
        sdr2.center_freq = slider.get()
        sdr3.center_freq = slider.get()
        sdr4.center_freq = slider.get()
        print(slider.get())



    def recording():
        var2.set(1)
        frequp1["state"] = "disabled"
        frequp2["state"] = "disabled"
        frequp3["state"] = "disabled"
        freqdown1["state"] = "disabled"
        freqdown2["state"] = "disabled"
        freqdown3["state"] = "disabled"
        slider["state"] = "disabled"
        if var2.get() == 1:
            while True:
                file1 = open("Record", "a")
                file2 = open("Record2", "a")
                file3 = open("Record3", "a")
                file4 = open("Record4", "a")
                samples = sdr.read_samples(1024)
                samples = samples.astype(np.complex64)  # Convert to 64
                print(type(samples))  # Verify it's 64
                samples.tofile(file1)
                samples2 = sdr2.read_samples(1024)
                samples2 = samples2.astype(np.complex64)  # Convert to 64
                print(type(samples2))  # Verify it's 64
                samples2.tofile(file2)
                samples3 = sdr3.read_samples(1024)
                samples3 = samples3.astype(np.complex64)  # Convert to 64
                print(type(samples3))  # Verify it's 64
                samples3.tofile(file3)
                samples4 = sdr4.read_samples(1024)
                samples4 = samples4.astype(np.complex64)  # Convert to 64
                print(type(samples4))  # Verify it's 64
                samples4.tofile(file4)
                ax1.text(100, 100, 'Record', color='red')
                if var2.get() == 0:
                    break
                file1.close()
                file2.close()
                file3.close()
                file4.close()

    def stoprec():
        var2.set(0)
        frequp1["state"] = "normal"
        frequp2["state"] = "normal"
        frequp3["state"] = "normal"
        freqdown1["state"] = "normal"
        freqdown2["state"] = "normal"
        freqdown3["state"] = "normal"
        slider["state"] = "normal"

    def store():
        channel = open("Channel", "w")
        channel.write(str(sdr.center_freq))
        channel.close()

    def tune1():
        channel = open("Channel", "r")
        chan = channel.readline()
        sdr.center_freq = chan
        sdr2.center_freq = chan
        sdr3.center_freq = chan
        sdr4.center_freq = chan
        channel.close()

    def store2():
        channel = open("Channel2", "w")
        channel.write(str(sdr.center_freq))
        channel.close()

    def tune22():
        channel = open("Channel2", "r")
        chan = channel.readline()
        sdr.center_freq = chan
        sdr2.center_freq = chan
        sdr3.center_freq = chan
        sdr4.center_freq = chan
        channel.close()

    def store3():
        channel = open("Channel3", "w")
        channel.write(str(sdr.center_freq))
        channel.close()

    def tune33():
        channel = open("Channel3", "r")
        chan = channel.readline()
        sdr.center_freq = chan
        sdr2.center_freq = chan
        sdr3.center_freq = chan
        sdr4.center_freq = chan
        channel.close()

    def threading():
        # Call work function
        record.config(text='Record', bg="red", activebackground="red")
        t1 = Thread(target=GUI.recording)
        t1.start()

    def threading2():
        # Call work function
        record.config(text='Record', bg="grey", activebackground="grey")
        recordstop.config(text='Stop', bg="grey", activebackground="grey")
        t1 = Thread(target=GUI.stoprec)
        t1.start()

    def threading3():
        # Call work function
        savestation.config(text=str(sdr.center_freq / 1e6), bg="grey", activebackground="orange")
        t1 = Thread(target=GUI.store)
        t1.start()

    def threading4():
        # Call work function
        tune.config(text='Tune', bg="grey", activebackground="orange")
        t1 = Thread(target=GUI.tune1)
        t1.start()

    def threading5():
        # Call work function
        savestation2.config(text=str(sdr.center_freq / 1e6), bg="grey", activebackground="orange")
        t1 = Thread(target=GUI.store2)
        t1.start()

    def threading6():
        # Call work function
        tune2.config(text='Tune', bg="grey", activebackground="orange")
        t1 = Thread(target=GUI.tune22)
        t1.start()

    def threading7():
        # Call work function
        savestation3.config(text=str(sdr.center_freq / 1e6), bg="grey", activebackground="orange")
        t1 = Thread(target=GUI.store3)
        t1.start()

    def threading8():
        # Call work function
        tune3.config(text='Tune', bg="grey", activebackground="orange")
        t1 = Thread(target=GUI.tune33)
        t1.start()



    def playback34():
        filename = 'Record'
        samples = np.fromfile(filename, np.complex64)
        samples2 = np.fromfile(filename + '2', np.complex64)
        samples3 = np.fromfile(filename + '3', np.complex64)
        samples4 = np.fromfile(filename + '4', np.complex64)
        strin = np.shape(samples)
        strin = str(strin)[1:-2]
        strin = int(strin)
        for i in range(1, strin, 1024):
            return samples[i:i + 1024], samples2[i:i + 1024], samples3[i:i + 1024], samples4[i:i + 1024]

    def threading9():
        # Call work function

        t1 = Thread(target=GUI.playback34)
        t1.start()

    def varsetter():
        var4.set(1)

    def play1(sam, sam2, sam3, sam4):
        ax1.psd(sam, Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                sides=None, data=None, color='yellow')
        ax1.psd(sam2, Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                sides=None, data=None, color='green')
        ax1.psd(sam3, Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                sides=None, data=None, color='orange')
        ax1.psd(sam4, Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                sides=None, data=None, color='blue')
        ax1.axhspan(-100, 40, xmin=0.5, xmax=0.5, color='white')
        ax1.set_xlabel('')
        # ax1.magnitude_spectrum(samples2)
        ax1.set_ylim(-50, 0)
        ax3.grid(visible=None, which='major', axis='both', color='black')
        ax3.set_theta_zero_location('N')
        ax3.set_theta_direction(-1)
        ax3.set_facecolor('green')
        ax3.set_xticklabels([])
        ax3.set_yticklabels([])
        t = np.arange(0, 1024, 1)
        ax2.plot(sam)
        ax2.plot(sam2)
        ax2.plot(sam3)
        ax2.plot(sam4)
        ax2.set_ylim(-2, 2)


    def recplay(sam, sam2, sam3, sam4):
        strin = np.shape(sam)
        strin = str(strin)[1:-2]
        strin = int(strin)
        for i in range(1, strin, 1024):
            ax1.clear()
            ax2.clear()
            ax1.psd(sam[i:i + 1024], Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                    sides=None, data=None, color='yellow')
            ax1.psd(sam2[i:i + 1024], Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                    sides=None, data=None, color='green')
            ax1.psd(sam3[i:i + 1024], Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                    sides=None, data=None, color='orange')
            ax1.psd(sam4[i:i + 1024], Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                    sides=None, data=None, color='blue')
            # ax1.magnitude_spectrum(gg1, Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, window=np.blackman(1024), pad_to=None,
            #                       sides=None, scale='dB', data=None, color='blue')
            # ax1.magnitude_spectrum(gg2, Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, window=np.blackman(1024), pad_to=None,
            #                       sides=None, scale='dB', data=None, color='green')
            # ax1.magnitude_spectrum(gg3, Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, window=np.blackman(1024), pad_to=None,
            #                       sides=None, scale='dB', data=None, color='red')
            # ax1.magnitude_spectrum(gg4, Fs=2.4e6 / 1e6, Fc=sdr.center_freq / 1e6, window=np.blackman(1024), pad_to=None,
            #                       sides=None, scale='dB', data=None, color='orange')
            ax1.axhspan(-100, 40, xmin=0.5, xmax=0.5, color='white')
            ax1.set_xlabel('')
            # \ax1.magnitude_spectrum(samples2)
            ax1.set_ylim(-50, 0)
            ax3.grid(visible=None, which='major', axis='both', color='black')
            ax3.set_theta_zero_location('N')
            ax3.set_theta_direction(-1)
            ax3.set_facecolor('green')
            ax3.set_xticklabels([])
            ax3.set_yticklabels([])
            t = np.arange(0, 1024, 1)
            ax2.plot(sam[i:i + 1024])
            ax2.set_ylim(-2, 2)


    def varsetter2():

        var4.set(0)


var3 = Tk.IntVar()
var4 = Tk.IntVar()



def animate(i):
    ax1.clear()
    ax2.clear()
    ax3.clear()




    samples23 = RDF.samples22()
    samples = samples23[0]
    samples2 = samples23[1]
    samples3 = samples23[2]
    samples4 = samples23[3]





    ###array antenna pairs
    northeast = np.array(abs(np.mean(samples)), abs(np.mean(samples2)))
    southeast = np.array(abs(np.mean(samples2)), abs(np.mean(samples3)))
    southwest = np.array(abs(np.mean(samples3)), abs(np.mean(samples4)))
    northwest = np.array(abs(np.mean(samples4)), abs(np.mean(samples)))
    ##antenna array
    antennaarr = np.array([northeast, southeast, southwest, northwest])
    maxarr = np.argmax(antennaarr)
    maxarr1 = maxarr or 0
    if maxarr is None:
        maxarr = 0
    ###Filter
    b, a = sc.butter(3, 0.05)
    # b, a = sc.cheby2(4, 40, 100, 'low', analog=True)
    filt = sc.lfilter(b, a, samples, axis=- 1, zi=None)
    filt2 = sc.lfilter(b, a, samples2, axis=- 1, zi=None)
    arr1 = antennaarr[maxarr]
    ###DOA Algorithm
    R = 10 * np.log(abs(np.mean(filt)) / abs(np.mean(filt2)))
    DOA = (31 ** 2 * R) / (20 * 0.69 * 90 * 0.4342944819)
    antcor = RDF.antsel(maxarr)
    antcor1 = antcor or 45
    def algorith(sel):
        if sel == 0:
            R = 10 * np.log(abs(np.mean(samples)) / abs(np.mean(samples2)))
            DOA = (31 ** 2 * R) / (20 * 0.69 * 90 * 0.4342944819)
            return DOA
        if sel == 1:
            R = 10 * np.log(abs(np.mean(samples2)) / abs(np.mean(samples3)))
            DOA1 = (31 ** 2 * R) / (20 * 0.69 * 90 * 0.4342944819)
            return DOA1
        if sel == 2:
            R = 10 * np.log(abs(np.mean(samples3)) / abs(np.mean(samples4)))
            DOA2 = (31 ** 2 * R) / (20 * 0.69 * 90 * 0.4342944819)
            return DOA2
        if sel == 3:
            R = 10 * np.log(abs(np.mean(samples4)) / abs(np.mean(samples)))
            DOA3 = (31 ** 2 * R) / (20 * 0.69 * 90 * 0.4342944819)
            return DOA3
    alg1 = algorith(maxarr)
    algorith1 = alg1 or DOA
    DOA1 = np.deg2rad(antcor1 - algorith1)
    ####Normalizing
    dis = np.array(RDF.normalize(samples, 0, 1))
    ###Filter
    b, a = sc.butter(3, 0.05)

    ax1.set_facecolor('black')
    ###SPectrum Plot
    ax3.vlines(DOA1, 0, 10, color='white')
    if var4.get() == 0:
        ax1.magnitude_spectrum(samples, Fs=3048000.222855 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                sides=None, data=None, color='red')
        ax1.magnitude_spectrum(samples2, Fs=3048000.222855 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                sides=None, data=None, color='blue')
        ax1.magnitude_spectrum(samples3, Fs=3048000.222855 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                sides=None, data=None, color='white')
        ax1.magnitude_spectrum(samples4, Fs=3048000.222855 / 1e6, Fc=sdr.center_freq / 1e6, pad_to=None,
                sides=None, data=None, color='grey')
        #ax1.axhspan(-100, 40, xmin=0.5, xmax=0.5, color='white')
        ax1.set_xlabel('')
        # ax1.magnitude_spectrum(samples2)
        #ax1.set_ylim(0.000, 0.398)
        ax3.grid(visible=None, which='major', axis='both', color='black')
        ax3.set_theta_zero_location('N')
        ax3.set_theta_direction(-1)
        ax3.set_facecolor('red')
        ax3.set_xticklabels([])
        ax3.set_yticklabels([])
        t = np.arange(0, 1024, 1)
        ax2.plot(samples, color='blue')
        ax2.plot(samples2, color='red')
        ax2.plot(samples3, color='white')
        ax2.plot(samples4, color='grey')
        ax2.set_ylim(-0.23, 0.23)
    else:
        pass
    if var4.get() == 1:
        filename = 'Record'
        samples = np.fromfile(filename, np.complex64)
        samples2 = np.fromfile(filename + '2', np.complex64)
        samples3 = np.fromfile(filename + '3', np.complex64)
        samples4 = np.fromfile(filename + '4', np.complex64)

        t1 = Process(target=GUI.recplay(samples,samples2,samples3,samples4))
        t1.start()
    else:
        pass

current_value = Tk.DoubleVar()
slider = Tk.Scale(
    root,
    from_=1e6,
    to=1766e6,
    orient='horizontal',
    variable=current_value,
    command=GUI.slider_changed,
    sliderlength=10,
    length=700

)
slider.grid(column=0, row=5)

frequp1 = Tk.Button(root, text='0.01mhz->', command=GUI.frequp1)
frequp1.grid(column=2, row=1)

freqdown1 = Tk.Button(root, text='<-0.01mhz', command=GUI.freqdown1)
freqdown1.grid(column=1, row=1)

freqdown2 = Tk.Button(root, text='<-25mhz', command=GUI.freqdown2)
freqdown2.grid(column=1, row=3)

frequp2 = Tk.Button(root, text='25mhz->', command=GUI.frequp2)
frequp2.grid(column=2, row=3)

freqdown3 = Tk.Button(root, text='<-0.1mhz', command=GUI.freqdown3)
freqdown3.grid(column=1, row=2)

frequp3 = Tk.Button(root, text='0.1mhz->', command=GUI.frequp3)
frequp3.grid(column=2, row=2)

c1 = Tk.Checkbutton(root, text='Butter Filter',variable=var1, onvalue=1, offvalue=0, command=RDF.samples22)
c1.grid(column=2, row=4)

record = Tk.Button(root, text='Record', command=GUI.threading, activebackground="red", state="normal", relief="raised")
record.grid(column=2, row=5)

recordstop = Tk.Button(root, text='Stop', command=GUI.threading2)
recordstop.grid(column=1, row=5)

savestation = Tk.Button(root, text='Store', command=GUI.threading3, activebackground="orange", state="normal", relief="raised")
savestation.grid(column=2, row=6)

tune = Tk.Button(root, text='Tune', command=GUI.threading4, activebackground="grey", state="normal", relief="raised")
tune.grid(column=1, row=6)

savestation2 = Tk.Button(root, text='Store', command=GUI.threading5, activebackground="orange", state="normal", relief="raised")
savestation2.grid(column=2, row=7)

tune2 = Tk.Button(root, text='Tune', command=GUI.threading6, activebackground="grey", state="normal", relief="raised")
tune2.grid(column=1, row=7)

savestation3 = Tk.Button(root, text='Store', command=GUI.threading7, activebackground="orange", state="normal", relief="raised")
savestation3.grid(column=2, row=8)

stop = Tk.Button(root, text='Stop', command=GUI.varsetter2, activebackground="orange", state="normal", relief="raised")
stop.grid(column=2, row=9)

play = Tk.Button(root, text='Play', command=GUI.varsetter, activebackground="orange", state="normal", relief="raised")
play.grid(column=1, row=9)

tune3 = Tk.Button(root, text='Tune', command=GUI.threading8, activebackground="grey", state="normal", relief="raised")
tune3.grid(column=1, row=8)

#playback11 = Tk.Button(root, text='Playback', command=playback34)
#playback11.grid(column=1, row=6)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(column=0, row=1)
ani = animation.FuncAnimation(fig, animate, interval=1)

Tk.mainloop()

sdr.stop()
sdr2.stop()
sdr3.stop()
sdr4.stop()

sdr.close()
sdr2.close()
sdr3.close()
sdr4.close()[4]