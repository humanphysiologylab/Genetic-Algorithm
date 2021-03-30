import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

import pandas as pd

x = [1000, 2000, 4000]
exclude = [[0,1,2], [0, 1], [0,1]]
filenames = ["20200916 - Atrial EHM - 1 Hz pacing - rec8 - 1 Hz.txt",\
        "20200916 - Atrial EHM - 1 Hz pacing - rec8 - 05 Hz.txt",\
        "20200916 - Atrial EHM - 1 Hz pacing - rec8 - 025 Hz.txt"]

paint_interval = [[46, 68], [1147, 1166], [898, 916]]

time_shift = [95,-1000,-750]

final_series = []

for i in range(len(x)):
    data = pd.read_csv(filenames[i], sep='\s+')
    data = pd.DataFrame(data)
    #plt.plot(data['t'], data['v'], label=str(x[i]) + " ms")
    #plt.legend()
    #plt.show()
    
      
    period = x[i]
    num = len(data['v'])//period
    b = data['v'].values[0: period * num ].reshape(period, num, order='F')
    a = pd.DataFrame(b)
    
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
    
    #f, Pwelch_spec = signal.welch(z, 10000, scaling='spectrum')
    #plt.semilogy(f, Pwelch_spec)
    #z = np.roll(z, time_shift[i])[-500:]
    f, Pper_spec = signal.periodogram(z, 1000)#, 'flattop', scaling='spectrum')
    plt.semilogy(f, Pper_spec)
    
    plt.xlabel('frequency [Hz]')
    plt.ylabel('PSD')
    plt.grid()
    plt.show()



for s in final_series:
    plt.plot(s, label=str(len(s)) + " ms")
    plt.legend()
plt.show()
#data = pd.read_csv("control_1s_atrial_EHM.txt")
#data = pd.DataFrame(data)["m1avg"]

#data.to_csv("control_1s_atrial_EHM_AVG.txt", index=False)