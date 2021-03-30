import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

from os import listdir
from os.path import isfile, join

import re
import json


def config_line(dirname):
    
    files = [join(dirname,f) for f in listdir(dirname)\
        if isfile(join(dirname, f)) and "table_" in f]

    p = re.compile("^table_(.*)_(.*)\.csv$")


    with open(sys.argv[1]) as f:
        config = json.load(f)

    glob = config['global']

    frame = pd.DataFrame()

    for x in glob:
        if "value" in x:
            continue
        filename = dirname + "/table_" + x['name'] + "_global.csv"
        files.remove(filename)
        try:
            with open(filename) as f:
                #save best organism value
                data = pd.read_csv(f, sep="\s+")
                data = pd.DataFrame(data)
                frame[x['name']] = [data["best_organism"].iloc[-1]]
        except IOError:
            print("cannot process " + filename)
            break

    baselines = config['baselines']

    for bl in baselines:
        for x in bl['params']:
            if "value" in x:
                continue
            filename = dirname + "/table_" + x['name'] +"_"+ bl['name'] + '.csv'
            files.remove(filename)
            try:
                with open(filename) as f:
                    #save best organism value
                    data = pd.read_csv(f, sep="\s+")
                    data = pd.DataFrame(data)
                    frame[x['name'] + "_" + bl['name']] = [data["best_organism"].iloc[-1]]

            except IOError:
                print("cannot process " + filename)
                break
    
    return frame 

dirs = ['r1', 'r2', 'r3', 'r4', 'r5']

frame = config_line(dirs[0])
for i in range(1, len(dirs)):
    frame = frame.append(config_line(dirs[i]))

print(frame)
plt.figure()

bp = frame.boxplot()
plt.xticks(rotation=90)
plt.show()
