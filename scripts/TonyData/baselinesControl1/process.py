import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

import pandas as pd

x = [1000, 2000, 4000]
exclude = [[0,15], [], [0,1]]
filenames = ["20200907 - Atrial EHM - Control - rec2 - 1 Hz.txt",\
        "20200907 - Atrial EHM - Control - rec2 - 05 Hz.txt",\
        "20200907 - Atrial EHM - Control - rec2 - 025 Hz.txt"]

paint_interval = [[517, 547], [716, 756], [2272, 2308]]

time_shift = [-303,-500,-2050]

final_series = []

for i in range(len(x)):
    data = pd.read_csv(filenames[i], sep='\s+')
    data = pd.DataFrame(data)
    plt.plot(data['t'], data['v'], label=str(x[i]) + " ms")
    plt.legend()
    
    
    
    fs = 1000
    n = 101
   # fil = signal.firwin(n, cutoff = [40, 60], window = 'blackmanharris', pass_zero = False, fs=1000)


   # fil = signal.firwin(n, cutoff = 45, window = "hamming", pass_zero='lowpass', fs=fs)
    fil = signal.firwin(n, cutoff = 45, window = "flattop", pass_zero='lowpass', fs=fs)
    filtered_v = signal.lfilter(fil, 1.0, data['v'])
    filtered_v = np.roll(filtered_v, -n//2)
   
   
    nyq_rate = fs / 2.0
    '''
    # The desired width of the transition from pass to stop,
    # relative to the Nyquist rate.  We'll design the filter
    # with a 5 Hz transition width.
    width = 5.0/nyq_rate
    # The desired attenuation in the stop band, in dB.
    ripple_db = 60.0
    # Compute the order and Kaiser parameter for the FIR filter.
    N, beta = signal.kaiserord(ripple_db, width)
    # The cutoff frequency of the filter.
    cutoff_hz = 70
    # Use firwin with a Kaiser window to create a lowpass FIR filter.
    fil = signal.firwin(N, cutoff_hz/nyq_rate, window=('kaiser', beta))   
    # Use lfilter to filter x with the FIR filter.
    filtered_v = signal.lfilter(fil, 1.0, data['v'])
    filtered_v = np.roll(filtered_v, -N//2)
    '''
    
    plt.plot(data['t'], filtered_v)
    plt.show()
    
    
    w, h = signal.freqz(fil)
    w = w / np.pi * fs/2
    fig, ax1 = plt.subplots()
    ax1.set_title('Digital filter frequency response')
    ax1.plot(w, 20 * np.log10(abs(h)), 'b')
    ax1.set_ylabel('Amplitude [dB]', color='b')
    ax1.set_xlabel('Frequency [rad/sample]')
    ax2 = ax1.twinx()
    angles = np.unwrap(np.angle(h))
    ax2.plot(w, angles, 'g')
    ax2.set_ylabel('Angle (radians)', color='g')
    ax2.grid()
    ax2.axis('tight')
    plt.show()
    
'''
    

    period = x[i]
    num = len(data['v'])//period
    b = data['v'].values[0: period * num ].reshape(period, num, order='F')
    a = pd.DataFrame(b)
    
    if len(exclude[i]) != 0:
    	ex = a[exclude[i]]
    	#plt.plot(ex, color="black")
    	#plt.legend(ex.columns.values, loc='upper left')
    	a = a.drop(exclude[i], axis=1)
    #plt.plot(a)
    #plt.legend(a.columns.values, loc='upper left')
    #plt.show()
    #plt.plot(a.mean(axis=1), color = "brown", linewidth=1, alpha=1)
    #plt.show()

    mean_v = a.mean(axis=1)
    start = paint_interval[i][0]
    end = paint_interval[i][1]
    for j in range(start, end + 1):
        mean_v[j] = (j - end) / (start - end) * mean_v[start] + (start - j) /\
        (start - end) * mean_v[end]
    mean_v = pd.Series(np.roll(mean_v.values, time_shift[i]), index=mean_v.index)
    #plt.plot(mean_v)
    #plt.show()
    final_series.append(mean_v)
    
    z = a[a.columns[0]].values
    for j in range(start, end + 1):
        z[j] = (j - end) / (start - end) * z[start] + (start - j) /\
        (start - end) * z[end]
    
    f, Pwelch_spec = signal.welch(z, 1000, scaling='spectrum')
    plt.semilogy(f, Pwelch_spec)
    #f, Pper_spec = signal.periodogram(z, 1000, 'flattop', scaling='spectrum')
    #plt.semilogy(f, Pper_spec)
    
    plt.xlabel('frequency [Hz]')
    plt.ylabel('PSD')
    plt.grid()
    plt.show()

for s in final_series:
    plt.plot(s, label=str(len(s)) + " ms")
    plt.legend()
plt.show()
'''
#data = pd.read_csv("control_1s_atrial_EHM.txt")
#data = pd.DataFrame(data)["m1avg"]

#data.to_csv("control_1s_atrial_EHM_AVG.txt", index=False)
