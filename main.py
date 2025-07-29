from functions import *

path = '/Volumes/BckSmoreau6/2025-NACA0015/Exp-NACA0015-2'
path_cp = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/Others/Mine/2025-ISAE-Exp/NACA0015/Cp/Exp-NACA0015-2'
filename = 'NACA0015_sweep-1.h5'
filename_cp = ['NACA0015_SVA_sweep-1.csv','NACA0015_SVB_sweep-1.csv',15,14]

chord_coord = np.array([0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.74, 0.8, 0.9, 0.85, 0.775, 0.7, 0.6, 0.5, 
                        0.4, 0.3, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.0125])
chanels = [13,14,15,16,25,26]
pr_i = 2
pr_n = 37
aoa = 0
ylim = [-20, 40]

# psd_generation(pr_i, pr_n, aoa, ylim, chanels, path, filename)
cp_generation(pr_i, pr_n, aoa, 0.37, 0.16, 8, chord_coord, path_cp, filename_cp)
# ladder_psd_generation(pr_i, pr_n, 0, chanels, path, 'NACA0015_sweep-1.h5')
# ladder_psd_generation(pr_i, pr_n, 1, chanels, path, 'NACA0015_1deg_sweep-1.h5')
# ladder_psd_generation(pr_i, pr_n, 2, chanels, path, 'NACA0015_2deg_sweep-1.h5')
# time_trace_generation(pr_i, pr_n, aoa, chanels, path, filename)