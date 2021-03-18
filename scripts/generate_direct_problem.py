import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

from os import listdir
from os.path import isfile, join

import re
import json

files = [f for f in listdir(".")\
        if isfile(join(".", f)) and "table_" in f]

p = re.compile("^table_(.*)_(.*)\.csv$")


with open(sys.argv[1]) as f:
    config = json.load(f)

#we prepare new config for direct problem
config['script'] = "Direct Problem"

glob = config['global']

for x in glob:
    if "value" in x:
        continue
    filename = "table_" + x['name'] + "_global.csv"
    try:
        with open(filename) as f:
            #save best organism value
            name = x['name']
            x.clear()
            x['name'] = name
            data = pd.read_csv(f, sep="\s+")
            data = pd.DataFrame(data)
            x['value'] = data["best_organism"].iloc[-1]

    except IOError:
        print("cannot process " + filename)
        break

baselines = config['baselines']

for bl in baselines:
    for x in bl['params']:
        if "value" in x:
            continue
        filename = "table_" + x['name'] +"_"+ bl['name'] + '.csv'
        try:
            with open(filename) as f:
                #save best organism value
                name = x['name']
                x.clear()
                x['name'] = name
                data = pd.read_csv(f, sep="\s+")
                data = pd.DataFrame(data)
                x['value'] = data["best_organism"].iloc[-1]

        except IOError:
            print("cannot process " + filename)
            break

#finally, save new config file
with open("direct_config.json", 'w') as newjson:
    json.dump(config, newjson, indent=4)


'''
for x in files:
    search = re.match(p, x)
    
    if search.group(2) == "global":
        config['global'].
    


    if search:
        print(search.group(1), search.group(2))




data = pd.read_csv(sys.argv[1], sep="\s+")
data = pd.DataFrame(data)


for column in data:
    if column == "time":
        continue
    plt.plot(data["time"] / 1000, data[column], label=column)
    plt.legend()
    plt.xlabel("time, s")
    plt.show()
'''
