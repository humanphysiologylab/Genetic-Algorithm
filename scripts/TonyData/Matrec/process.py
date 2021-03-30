import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import pandas as pd
from dataclasses import dataclass
import scipy.io

@dataclass
class Data:
    periods: list
    exclude: list
    prefix: str
    filenames: list
    paint_intervals: list # s
    paint_intervals2: list #s
    shifts: list

def kaiser_filter(samplerate, cutoff_freq_hz):
    nyq_rate = samplerate / 2.0
    # The desired width of the transition from pass to stop,
    # relative to the Nyquist rate.  We'll design the filter
    # with a 5 Hz transition width.
    width = 5.0 / nyq_rate
    # The desired attenuation in the stop band, in dB.
    ripple_db = 60.0
    # Compute the order and Kaiser parameter for the FIR filter.
    N, beta = signal.kaiserord(ripple_db, width)
    # Use firwin with a Kaiser window to create a lowpass FIR filter.
    fil = signal.firwin(N, cutoff_freq_hz/nyq_rate, window=('kaiser', beta))
    return N, fil

def analyze_filter(fil, fs):
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

def freq_analysis(functions, samplerate):
    for f in functions:
        f, Pwelch_spec = signal.welch(f, samplerate, nperseg=2048, scaling='spectrum')
        plt.semilogy(f, Pwelch_spec)
        #f, Pper_spec = signal.periodogram(v_data, samplerate, scaling='spectrum')
        #plt.semilogy(f, Pper_spec)
        
    plt.xlabel('frequency [Hz]')
    plt.ylabel('PSD')
    plt.grid()
    plt.show()

def custom_comb_filtering(w0, samplerate, si):
    res = si.copy()
    Q = 1000
    for w in np.arange(w0, samplerate/2 + 0.001, 2*w0):
        if w == 150:
            b, a = signal.iirnotch(w, 100, fs=samplerate)
        else:
            b, a = signal.iirnotch(w, int(Q * (w0 / w)**0.5), fs=samplerate)
        res = signal.filtfilt(b, a, res)
    return res


def meanGroup(time, functions, period, samplerate, exclude):
    '''
    functions to plot
    period, s
    samplerate, Hz
    '''
    
    for i, f in enumerate(functions):
        samplesize = int(period * samplerate)
        num = len(f)//samplesize
        b = f[0: samplesize * num].reshape(samplesize, num, order='F')
        a = pd.DataFrame(b)

        #exclude maybe
        if len(exclude) != 0:
            ex = a[exclude]
            plt.plot(time[:samplesize], ex, color="black")
            plt.legend(ex.columns.values, loc='upper left')
            a = a.drop(exclude, axis=1)

        mean_v = a.mean(axis=1)
        plt.plot(time[:samplesize], a)
        plt.legend(a.columns.values, loc='upper left')
        plt.plot(time[:samplesize], mean_v, label=str(i))

    plt.show()

def getMean(time, f, period, samplerate, exclude):
    '''
    period, s
    samplerate, Hz
    '''
    
    samplesize = int(period * samplerate)
    num = len(f)//samplesize
    b = f[0: samplesize * num].reshape(samplesize, num, order='F')
    a = pd.DataFrame(b)

    #exclude maybe
    if len(exclude) != 0:
        ex = a[exclude]
        a = a.drop(exclude, axis=1)

    mean_v = a.mean(axis=1)
    return mean_v

def paint(v_data, startx, endx, period, samplerate):
    '''
    all in s
    '''
    for k in range(len(v_data) // int(samplerate * period)):
        r_start = int((period * k + startx) * samplerate)
        r_end = int((period * k + endx) * samplerate)
        for j in range(r_start, r_end + 1):
            v_data[j] = (j - r_end) / (r_start - r_end) * v_data[r_start] +\
                (r_start - j) / (r_start - r_end) * v_data[r_end]

def process(data):
    filenames = [data.prefix + "/" + s for s in data.filenames]
    periods = data.periods
    exclude = data.exclude
    paint_intervals = data.paint_intervals
    paint_intervals2 = data.paint_intervals2
    shifts = data.shifts
    

    final_series = []
    for i in range(len(periods)):
        #v_data = pd.read_csv(filenames[i], sep='\s+')
        #v_data = pd.DataFrame(v_data)
        
        mat = scipy.io.loadmat(filenames[i])

        samplerate = int(mat['samplerate'][0][0])
        assert samplerate == 40000
    
        v_data = np.array(mat['data'][0], dtype='float64')
        
        time_data = np.arange(0, len(v_data)) / samplerate # in s
        #plt.plot(time_data, v_data, label=str(periods[i]) + " ms")
        #plt.vlines(np.arange(0.97,30,1), -80, 50, color='r', zorder=3  )
        #plt.show()
        #paint
        startx = paint_intervals[i][0]
        endx = paint_intervals[i][1]
        period = periods[i]/1000 # s

        paint(v_data, startx, endx, period, samplerate)

        #downsample
        oldsamplerate = samplerate
        samplerate = 1000
        v_data = scipy.signal.decimate(v_data, oldsamplerate//samplerate,
        ftype='fir')
        time_data = np.arange(0, len(v_data)) / samplerate # in s
       # plt.plot(time_data, v_data, label=str(periods[i]) + " ms, downsampled")
       # plt.legend()
       # plt.show()
        




        N = int(samplerate/10)
        N += N%2 + 1 #make N odd
        window = 'blackmanharris'
        #window = 'hamming'
        #window = 'flattop'
        fil = signal.firwin(N, cutoff = 45, window=window,
                pass_zero = 'lowpass', fs=samplerate)

        filtered_v = signal.lfilter(fil, 1.0, v_data)
        filtered_v = np.roll(filtered_v, -N//2)


        b, a = signal.iirnotch(50, 30, fs=samplerate)
        #b, a = signal.iircomb(w0 = 50, Q = 100, ftype='notch', fs=samplerate)
        filtered_v2 = signal.filtfilt(b, a, v_data)
        #b, a = signal.iirpeak(w0 = 0, Q = 30, fs=samplerate)
        #base = signal.filtfilt(b, a, v_data)
        #filtered_v += base

        
        filtered_v = custom_comb_filtering(50, samplerate, v_data)
        #filtered_v = v_data





        #filtered_v = signal.filtfilt(b, a, filtered_v)
        #plt.plot(time_data, filtered_v2, label='filteredv2', alpha=0.7)
       # plt.plot(time_data, filtered_v, label='filtered', alpha=0.7)
       # plt.plot(time_data, v_data, label='orig', alpha=0.7)
       # plt.legend()
       # plt.show()
       # freq_analysis([v_data, filtered_v, filtered_v2], samplerate)
        #freq_analysis([v_data, filtered_v], samplerate) 
        #plt.plot(time_data, filtered_v)
        #plt.show()
        

        #meanGroup(time_data, [filtered_v], period, samplerate, exclude[i])
        mean_v = getMean(time_data, filtered_v, period, samplerate, exclude[i])


        #additional filtering
        #if filenames[i] =="1/20200916 - Atrial EHM - 1 Hz pacing - rec8 - 1 Hz.mat":
        #    b, a = signal.iirnotch(50, 15, fs=samplerate)
            #b, a = signal.iircomb(w0 = 50, Q = 100, ftype='notch', fs=samplerate)
        #    mean_v = signal.filtfilt(b, a, mean_v)
        #    print("OK")
            #mean_v = custom_comb_filtering(50, samplerate, mean_v)


        paint(mean_v, paint_intervals2[i][0], paint_intervals2[i][1], period, samplerate)
#        plt.plot(time_data[:int(period * samplerate)], mean_v,
#                label=filenames[i])
#        plt.legend()
#        plt.show()
#        freq_analysis([v_data, filtered_v, mean_v], samplerate) 
        
        
        mean_v = pd.Series(np.roll(mean_v, int(shifts[i]/1000*samplerate)))#, index=mean_v.index)
        final_series.append(mean_v)

    for i, s in enumerate(final_series):
        plt.plot(s, label=str(periods[i]) + " ms")
        plt.legend()
    plt.show()
  #  for i, s in enumerate(final_series):
  #      pd.Series(s).to_csv(filenames[i].replace('mat', 'txt'), index=False)
    #data = pd.read_csv("control_1s_atrial_EHM.txt")
    #data = pd.DataFrame(data)["m1avg"]

    #data.to_csv("control_1s_atrial_EHM_AVG.txt", index=False)

hz1_data = Data(
        [1000, 2000, 4000],
        [[0],[0],[0,1]],
        "1",
        ["20200916 - Atrial EHM - 1 Hz pacing - rec8 - 1 Hz.mat",
         "20200916 - Atrial EHM - 1 Hz pacing - rec8 - 05 Hz.mat",
         "20200916 - Atrial EHM - 1 Hz pacing - rec8 - 025 Hz.mat"
        ],
        [[0.05, 0.067], [1.15, 1.167], [0.9, 0.915]],
        [[0.05, 0.067], [1.15, 1.167], [0.9, 0.9166]],
        [100,-1000,-745]
)
hz3_data = Data(
        [1000, 500, 2000],
        [[0,1,2],[0,1],[0,1]],
        "3",
        ["20200903 - Atrial EHM - 3 Hz pacing - rec18 - 1 Hz.mat",
         "20200903 - Atrial EHM - 3 Hz pacing - rec18 - 2 Hz.mat",
         "20200903 - Atrial EHM - 3 Hz pacing - rec18 - 05 Hz.mat"
        ],
        [[0.425, 0.445], [0.05, 0.07], [0.3, 0.318]],
        [[0.425, 0.47], [0.05, 0.078], [0.3, 0.343]],
        [-350,0,-220]
)
control2_data = Data(
    [1000, 2000, 4000],
    [[0,1,15],[0],[0,1]],#[[0,15], [], [0,1]]
    "c2",
    ["20200907 - Atrial EHM - Control - rec2 - 1 Hz.mat",
    "20200907 - Atrial EHM - Control - rec2 - 05 Hz.mat",
    "20200907 - Atrial EHM - Control - rec2 - 025 Hz.mat"],
    [[0.525, 0.54], [0.72, 0.735], [2.27, 2.285]],
    [[0.517, 0.549], [0.72, 0.748], [2.27, 2.297]],
    [-303,-500,-2050]
)
'''
control1_data = Data(
    [1000, 500, 333.3333333333],
    [],#[[0,15], [], [0,1]]
    "c1",
    ["20200420 - Atrial EHM - Control - rec10 - 1 Hz.mat",
    "20200420 - Atrial EHM - Control - rec10 - 2 Hz.mat",
    "20200420 - Atrial EHM - Control - rec10 - 3 Hz.mat"],
    [[0.85, 0.86], [0.235, 0.25], [0.34, 0.35]],
    [[0.85, 0.86], [0.235, 0.25], [0.34, 0.35]],
    [-303,-500,-2050]
)
'''
#process(control1_data)
process(hz3_data)
process(hz1_data)
process(control2_data)
