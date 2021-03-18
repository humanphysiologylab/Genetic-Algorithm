import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


save_dir = sys.argv[-1]
p = [500, 1000, 2000]
for i in range(1, len(sys.argv)-1):
    d = pd.read_csv(sys.argv[i], sep="\s+")
    d = pd.DataFrame(d)
    d["V"][-p[i-1]: -1].to_csv(save_dir + "/" +  str(p[i-1]) + ".txt", sep=" ", index=False)
    plt.plot(d["V"][-p[i-1]: -1], label = str(p[i-1]))
plt.legend()
plt.show()
