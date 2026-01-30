from functions import *

path = '/scratch/renj3003/cd-airfoil/data'
path_cp = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/Others/Mine/2025-ISAE-Exp/NACA0015/Cp/Exp-NACA0015-2'
filename = 'hover_AD30_200-1.h5'
# filename_cp = ['NACA0015_6deg_SVA_sweep2-1.csv','NACA0015_6deg_SVB_sweep2-1.csv',15,14]
filename_cp = ['NACA0012_SVA_sweep-1.csv','NACA0012_SVB_sweep-1.csv',15,16]

# chord_coord = np.array([0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.74, 0.8, 0.9, 0.85, 0.775, 0.7, 0.6, 0.5, 
#                         0.4, 0.3, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.0125])
naca12_coord = np.array([0.009230769,	0.025384615,	0.043076923,	0.061538462,	0.092307692,	0.123076923, 0.153846154,
                         0.192307692,	0.230769231,	0.269230769,	0.307692308,	0.346153846, 0.423076923, 0.5, 0.576923077,
                         0, 0.061538462, 0.807692308, 0.884615385, 0.95, 1, 0.95, 0.8, 0.65, 0.5, 0.35, 0.025, 0.19, 0.12, 0.06, 0.025])
polar_angles = np.arange(0, 184, 3)  # Polar angles from 0 to 180 degrees in steps of 5])
chanels = [90]
pr_i = 26
pr_n = 28
folder_name = 'CD-0deg'
ylim = [-10, 50]
ylim_2 = [20, 70]
pr_it = [25,26,27,28,29,30]

#h5_2_mat(pr_i, pr_n, 1, 91, chanels, path, 'CD_ISAE_5deg-14.h5')
#psd_generation(pr_i, pr_n, folder_name, 1, 91, ylim, chanels, path, 'CD_ISAE_5deg-14.h5')
psd_comparison(pr_it, folder_name, 1, 91, ylim, chanels, path, 'CD_ISAE_5deg-14.h5')
#spectrogram_generation(pr_i, pr_n, folder_name, 1, 91, ylim, chanels, path, 'CD_ISAE_5deg-14.h5')
#wavelet_generation(pr_i, pr_n, folder_name, 1, 121, ylim, chanels, path, 'CD_ISAE_0deg-14.h5')
# psd_generation(pr_i, 3, 'AD30_250', 1, ylim, chanels, path, 'hover_AD30_250-1.h5')
# psd_generation(pr_i, 3, 'AD30_250', 1, ylim, chanels, path, 'hover_AD30_250-1.h5')
# psd_generation(pr_i, 3, 'ADS01', 1, ylim, chanels, path, 'hover_ADS01-1.h5')
# psd_generation(pr_i, 3, 'ADS10', 1, ylim, chanels, path, 'hover_ADS10-1.h5')
# psd_generation(pr_i, 3, 'ADS16', 1, ylim, chanels, path, 'hover_ADS16-1.h5')
# psd_generation(pr_i, 3, 'ADS17', 1, ylim, chanels, path, 'hover_ADS17-1.h5')
# psd_generation(pr_i, 3, 'ADS18', 1, ylim, chanels, path, 'hover_ADS18-1.h5')
# psd_generation(pr_i, 3, 'PT03', 1, ylim, chanels, path, 'hover_PT03-1.h5')
# psd_generation(pr_i, 9, 'AD30_200_cross', 1, 8, ylim, chanels, path, 'Cross_AD30_200-1.h5')
# psd_generation(pr_i, 6, 'PT08', 1, ylim, chanels, path, 'hover_PT08-1.h5')

# directivity_generation(pr_i, pr_n, 'AD30_200', 1, ylim, chanels, polar_angles, path, 'hover_AD30_200-1.h5', [100,10000])
# directivity_generation(pr_i, 3, 'AD30_250', 1, ylim, chanels, polar_angles, path, 'hover_AD30_250-1.h5', [100,10000])
# directivity_generation(pr_i, 3, 'ADS01', 1, ylim, chanels, polar_angles, path, 'hover_ADS01-1.h5', [100,10000])
# directivity_generation(pr_i, 3, 'ADS10', 1, ylim, chanels, polar_angles, path, 'hover_ADS10-1.h5', [100,10000])
# directivity_generation(pr_i, 3, 'ADS16', 1, ylim, chanels, polar_angles, path, 'hover_ADS16-1.h5', [100,10000])
# directivity_generation(pr_i, 3, 'ADS17', 1, ylim, chanels, polar_angles, path, 'hover_ADS17-1.h5', [100,10000])
# directivity_generation(pr_i, 3, 'ADS18', 1, ylim, chanels, polar_angles, path, 'hover_ADS18-1.h5', [100,10000])
# directivity_generation(pr_i, 3, 'PT03', 1, ylim, chanels, polar_angles, path, 'hover_PT03-1.h5', [100,10000])
# directivity_generation(pr_i, 5, 'PT05', 1, ylim_2, chanels, polar_angles, path, 'hover_PT05-1.h5', [100,10000])
# directivity_generation(pr_i, 9, 'AD30_200_cross', 1, 8, ylim_2, chanels, polar_angles, path, 'Cross_AD30_200-1.h5', 'integrated', [100,10000])
# cp_generation(pr_i, pr_n, aoa, 0.37, 0.16, 0, 1, naca12_coord, path, filename_cp)
# ladder_psd_generation(pr_i, pr_n, aoa, 2, chanels, path, 'NACA0012_0deg-1.h5')
# ladder_psd_generation(pr_i, pr_n, 1, 2, chanels, path, 'NACA0015_1deg_sweep-1.h5')
# ladder_psd_generation(pr_i, pr_n, 2, 2, chanels, path, 'NACA0015_2deg_sweep-1.h5')
#time_trace_generation(pr_i, pr_n, 0, 1, chanels, path, 'CD_ISAE_0deg-14.h5')
