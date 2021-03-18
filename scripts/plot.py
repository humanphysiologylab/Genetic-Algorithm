import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

files = ["ap0.txt",
        'ap1.txt',
        'ap2.txt',
        '1/20200916 - Atrial EHM - 1 Hz pacing - rec8 - 1 Hz.txt',
        '1/20200916 - Atrial EHM - 1 Hz pacing - rec8 - 05 Hz.txt',
        '1/20200916 - Atrial EHM - 1 Hz pacing - rec8 - 025 Hz.txt'
        ]

for f in files:
    data = pd.read_csv(f, sep='\s+', header=None)
    data = pd.DataFrame(data)
    xx = data[0]
    plt.plot(np.array(range(0, len(xx))), xx, label='exp' if f[:2] == '1/' else 'ga')

plt.legend()
plt.show()

data = pd.read_csv("convergence_GA.txt", sep='\s+', header=None)
data = pd.DataFrame(data)
plt.plot(data[0], data[1], label="error" )
plt.yscale('log')
plt.legend()
plt.show()
