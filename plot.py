#!/usr/bin/env python3

from sys import argv, exit
import pandas as pd
import matplotlib.pyplot as plt

if len(argv) != 2:
    exit("Usage: ./plot filename.dat")

data = pd.read_csv(argv[1], sep='\t', names=['f', 'r', 'i'])

fig, axs = plt.subplots(2, 1)
axs[0].plot(data.f, data.r)
axs[0].grid()
axs[0].set_xscale('log')
axs[1].plot(data.f, -data.i)
axs[1].grid()
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[0].set_xlabel('f, Hz')
axs[1].set_xlabel('f, Hz')
axs[0].set_title(r'$\Re{Z}$')
axs[1].set_title(r'$-\Im{Z}$')
plt.tight_layout()
plt.show()
