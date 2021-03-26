import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


fields = {
        "Nai",
        "Ki",
        "Ca_SR",
        "Cai"
        }

data = []
for i in range(1, len(sys.argv)):
    d = pd.read_csv(sys.argv[i], sep="\s+")
    data.append(pd.DataFrame(d))



s = 1
for column in data[0]:
    if column not in fields:
        continue
    plt.subplot(2,2, s)
    s+=1
    for i in range(len(data)):
        #plt.plot((data[i]["time"]) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])
        
        if data[i]["time"][0] > 100:
           plt.plot((data[i]["time"] - data[i]["time"][0] ) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])
        else:
            plt.plot((data[i]["time"]) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])
    plt.legend()
    plt.xlabel("time, s")
    plt.tight_layout()
plt.show()


