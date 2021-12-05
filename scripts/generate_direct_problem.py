import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

from os import listdir
from os.path import isfile, join

import re
import json

files_path = sys.argv[2]

files = [f for f in listdir(files_path)\
        if isfile(join(files_path, f)) and "table_" in f]

print(files)

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
    files.remove(filename)
    try:
        with open(join(files_path, filename)) as f:
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
        files.remove(filename)
        try:
            with open(join(files_path, filename)) as f:
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

# the rest of the files is for RESET_STATE = 0 case
# fields are not explicitly stated in config
for x in files:
    search = re.match(p, x)
    assert search
    print(search.group(1), search.group(2))
    bs = next(b for b in baselines if b['name'] == search.group(2))
    data = pd.read_csv(join(files_path, x), sep="\s+")
    data = pd.DataFrame(data)
    val = data["best_organism"].iloc[-1]
    bs['params'].append({'name': search.group(1), 'value': val})


#finally, save new config file
with open("direct_config.json", 'w') as newjson:
    json.dump(config, newjson, indent=4)

