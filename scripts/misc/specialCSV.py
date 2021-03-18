import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

data = []
d = pd.read_csv(sys.argv[1], sep="\s+", header=None)
d = pd.DataFrame(d)

for i in range(2, len(sys.argv)):
    de = pd.read_csv(sys.argv[i], sep="\s+")
    data.append(pd.DataFrame(de))


for column in data[0]:
    if column == "time":
        continue

    for i in range(len(data)):
        #plt.plot((data[i]["time"]) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])

        if data[i]["time"][0] > 100:
           plt.plot((data[i]["time"] - data[i]["time"][0] ) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])
        else:
            plt.plot((data[i]["time"]) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])
    if column == "V":
        plt.plot(np.arange(0, len(d[0])/1000, 0.001), d[0], label="ga input")

    plt.legend()
    plt.xlabel("time, s")
    plt.show()
