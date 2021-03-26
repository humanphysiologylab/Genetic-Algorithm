import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


data = []
for i in range(1, len(sys.argv)):
    d = pd.read_csv(sys.argv[i], sep="\s+")
    data.append(pd.DataFrame(d))

files = {
'3/20200903 - Atrial EHM - 3 Hz pacing - rec18 - 05 Hz.txt',
'3/20200903 - Atrial EHM - 3 Hz pacing - rec18 - 1 Hz.txt',
'3/20200903 - Atrial EHM - 3 Hz pacing - rec18 - 2 Hz.txt'
}

for column in data[0]:
    if column == "time":
        continue
    if column == "V":
        for f in files:
            d = pd.read_csv(f, sep='\s+', header=None)
            d = pd.DataFrame(d)
            plt.plot(np.array(range(0, len(d[0]))) / 1000, d[0], label='exp')
    for i in range(len(data)):
        #plt.plot((data[i]["time"]) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])

        if data[i]["time"][0] > 100:
           plt.plot((data[i]["time"] - data[i]["time"][0] ) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])
        else:
            plt.plot((data[i]["time"]) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])
    plt.legend()
    plt.xlabel("time, s")
    plt.show()
