from functions import *

path = '/Volumes/BckSmoreau6/2025-NACA0015'
filename = 'NACA0015_2deg_sweep-1.h5'

chanels = [13,14,15,16,25,26]
pr_i = 1
pr_n = 47
aoa = 2
ylim = [-20, 40]

psd_generation(pr_i, pr_n, aoa, ylim, chanels, path, filename)
# time_trace_generation(pr_i, pr_n, aoa, chanels, path, filename)