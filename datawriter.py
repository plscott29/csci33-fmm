import pandas as pd
import numpy as np
import sys

# example to run: python datawriter.py 1000 data3.csv

N = int(sys.argv[1])
filename = sys.argv[2]

x = np.float32(np.random.random_sample((N, 2)))
q = np.random.randint(-5, 6, (N, 1))
data = np.concatenate((x, q), axis=1)
df = pd.DataFrame(data, columns=['x', 'y', 'q'])
df['q'] = df['q'].astype(int)
df['x'] = df['x'].round(8)
df['y'] = df['y'].round(8)
df.to_csv(filename, header = False, index=False)