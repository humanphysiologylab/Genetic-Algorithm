import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

from matplotlib.backends.backend_pdf import PdfPages


data = []
for i in range(1, len(sys.argv)):
    d = pd.read_csv(sys.argv[i], sep="\s+")
    data.append(pd.DataFrame(d))


num_columns = len(data[0])

s = 1
with PdfPages('multipage_pdf.pdf') as pdf:
    for column in data[0]:
        if column == 'time':
            continue
        plt.figure(figsize=(3, 3))
	    #plt.subplot(1, 2, s)#num_columns//2 + 1, 2, s)
        #s += 1
        for i in range(len(data)):
        #plt.plot((data[i]["time"]) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])
            if data[i]["time"][0] > 100:
                plt.plot((data[i]["time"] - data[i]["time"][0] ) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])
            else:
                plt.plot((data[i]["time"]) / 1000, data[i][column], label=column + "_" + sys.argv[1+i])
            plt.legend(prop={'size': 6})
            plt.xlabel("time, s")
            plt.tight_layout()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
