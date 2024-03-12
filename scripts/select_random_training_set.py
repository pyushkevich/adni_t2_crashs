import glob
import os
import json
import pathlib
from sklearn.model_selection import KFold
import pandas as pd
import argparse
import random

parser = argparse.ArgumentParser(prog='upsample_net')
parser.add_argument('--seed', type=int, default=None)
parser.add_argument('-n', type=int, default=5)
parser.add_argument('manifest', type=argparse.FileType('rt'))
parser.add_argument('output', type=argparse.FileType('wt'))
args = parser.parse_args()

# Read all the ids
mf = pd.read_csv(args.manifest)
print(mf)
id_all = mf.iloc[:,0].tolist()
id_sided = [ (k,'left') for k in id_all ] + [ (k,'right') for k in id_all ]

# Get the ids we want
r = random.Random(args.seed)
r.shuffle(id_sided)
id_sided_sel = id_sided[0:args.n]

# Create a dataframe to save
res = pd.DataFrame({'id': [k[0] for k in id_sided_sel], 'side': [k[1] for k in id_sided_sel]})
res.to_csv(args.output, header=False, index=False)

