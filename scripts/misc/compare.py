import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


x = ["1000.txt", "500.txt", "2000.txt"]
for i in x:
    d = pd.read_csv("longest/" + i, sep="\s+", header=None)
    plt.plot(d, label="10000 s")
 
    d = pd.read_csv("new_bs/" + i, sep="\s+", header=None)
    plt.plot(d, label="100 s")
    d = pd.read_csv("longer_bs/" + i, sep="\s+", header=None)
    plt.plot(d, label="1000 s")

    plt.legend()
    plt.show()
