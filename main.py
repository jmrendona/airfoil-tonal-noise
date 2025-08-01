from functions import *

path = '/Volumes/BckSmoreau6/NACA0012_localisation-2025/Exp-NACA0012-0deg'
path_cp = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/Others/Mine/2025-ISAE-Exp/NACA0015/Cp/Exp-NACA0015-2'
filename = 'NACA0015_sweep-1.h5'
# filename_cp = ['NACA0015_6deg_SVA_sweep2-1.csv','NACA0015_6deg_SVB_sweep2-1.csv',15,14]
filename_cp = ['NACA0012_SVA_sweep-1.csv','NACA0012_SVB_sweep-1.csv',15,16]

# chord_coord = np.array([0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.74, 0.8, 0.9, 0.85, 0.775, 0.7, 0.6, 0.5, 
#                         0.4, 0.3, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.0125])
naca12_coord = np.array([0.009230769,	0.025384615,	0.043076923,	0.061538462,	0.092307692,	0.123076923, 0.153846154,
                         0.192307692,	0.230769231,	0.269230769,	0.307692308,	0.346153846, 0.423076923, 0.5, 0.576923077,
                         0, 0.061538462, 0.807692308, 0.884615385, 0.95, 1, 0.95, 0.8, 0.65, 0.5, 0.35, 0.025, 0.19, 0.12, 0.06, 0.025])
chanels = [13,14,15,16,25,26]
pr_i = 1
pr_n = 1
aoa = 0
ylim = [-20, 40]

# psd_generation(pr_i, pr_n, aoa, 2, ylim, chanels, path, filename)
cp_generation(pr_i, pr_n, aoa, 0.37, 0.16, 8, 1, naca12_coord, path, filename_cp)
# ladder_psd_generation(pr_i, pr_n, 0, 2, chanels, path, 'NACA0015_sweep-1.h5')
# ladder_psd_generation(pr_i, pr_n, 1, 2, chanels, path, 'NACA0015_1deg_sweep-1.h5')
# ladder_psd_generation(pr_i, pr_n, 2, 2, chanels, path, 'NACA0015_2deg_sweep-1.h5')
# time_trace_generation(pr_i, pr_n, aoa, 2, chanels, path, filename)