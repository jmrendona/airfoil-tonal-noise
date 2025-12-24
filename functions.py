import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os
from scipy.io import savemat
from scipy import stats
import h5py
import scipy.signal as sg
import pdb
import pywt
# from rich import print as rprint

def next_greater_power_of_2(x):
    return round(2**(x-1).bit_length())

def time_trace_generation(pr_i:int,pr_n:int,aoa:int,jump:int,chanels:list,path:str,filename:str):
    
    '''
    This function could be used to quickly evaluate the signals obtained from the B&K system\n
    after they are extracted to HDF5 format.\n
    Previous knowledge of the key names is needed in order to do the  check. The first key name\n
    is related with the time data and is as follows: 'DsX1-Time', where X is the proper zero-padding\n
    to match the number of signals used. If you have less then 100 it will be Ds01-Time, otherwise\n
    Ds-001 is used. The key name pressure refers to the name B&K uses for each signal. Usually it \n
    is defined as 'DsZZ-Signal YY' where XX is related with the signal number saved for the particular\n
    case with its associated zero-padding, if a proper saving is done XX=02 if less than 100 signals are\n
    recorded. YY on the other side refers to the used port on the B&K system, i.e. if the master card is used\n
    the first signal will have YY=1. 
    
    - path = string with the full path where the data is stored. Note that the image will be stored there\n
    - filename = string with the name of the file inclduing the .h5 extention\n
    - key_name_time = string with the key name for the time data\n
    - key_name_pressure = string with the key name for the mic to analyze\n
    - image_name = string with the name for the image that wil be saved\n
     
    '''
    
    plt.rcParams['agg.path.chunksize'] = 10000
    for k in chanels: 
        os.makedirs(os.path.join(path, f'images/timetrace/{aoa}deg/signal_{k}'), exist_ok=True)
    
    filename = filename.split('-')[0]
    pr_used = np.arange(pr_i,pr_n+1,jump)
    
    for i in pr_used:
        
        k=124
        for j in chanels:
            
            data = h5py.File(os.path.join(path,filename + f'-{i}.h5'),'r')
            time_data = data['Table1']['Ds001-Time'][:]
            pressure_data = data['Table1'][f'Ds{k:03d}-Signal {j}'][:]
            
            # Time trace plot
            plt.figure()
            plt.plot(time_data,pressure_data,'k--')
            plt.grid(True, which='both', ls='--')
            plt.xlabel('Time $[s]$',fontsize=12, style='italic')
            plt.ylabel('Pressure $[Pa]$',fontsize=12, style='italic')
            plt.ylim([np.min(pressure_data)*5,np.max(pressure_data)*5])
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.tight_layout()
            print(f'Saving timetrace for file {i} and signal {j} \n')
            print(40*'-')
            plt.savefig(os.path.join(path, f'images/timetrace/{aoa}deg/signal_{j}', filename + f'-{i}pr-s{j}-timetrace.png') , dpi=600)
            plt.close()
            
            k += 1
            

def psd_generation(pr_i:int,pr_n:int,folder_name:str,jump:int,offset:int,y_lim:list,chanels:list,path:str,filename:str):
    
    plt.rcParams['agg.path.chunksize'] = 10000
    
    for k in chanels: 
        os.makedirs(os.path.join(path, f'images/psd/{folder_name}'), exist_ok=True)
    
    filename = filename.split('-')[0]
    pr_used = np.arange(pr_i,pr_n+1,jump)
    
    for i in pr_used:
    
        k=offset
        for j in chanels:
            
            data = h5py.File(os.path.join(path,filename + f'-{i}.h5'),'r')
            time_data = data['Table1']['Ds001-Time'][:]
            pressure_data = data['Table1'][f'Ds{k:03d}-Signal {j}'][:]
        
            dt=time_data[1]-time_data[0]
            fs=1.0/dt
            
            fmin = 10
            n_chunk = 20
            lensg_exp = pressure_data.size
            nperseg_exp = lensg_exp/n_chunk
            nfft_exp = next_greater_power_of_2(int(nperseg_exp))
            noverlap_exp = nperseg_exp/2

            if nperseg_exp > lensg_exp:
                raise RuntimeError('Wrong value for $f_{min}$')

            [f,Pxx]=sg.welch(pressure_data,fs=fs,window='hann',nperseg=nperseg_exp,nfft=nfft_exp,scaling='density')
                
            # PSD plot
            plt.plot(f,10*np.log10(Pxx/4.0e-10),'k--')
            plt.grid(True, which='both', ls='--')
            plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
            plt.ylabel('PSD $\\left[\\frac{dB}{Hz}\\right]$',fontsize=12, style='italic')
            plt.xlim([70,20000])
            plt.ylim(y_lim)
            ax=plt.gca()
            ax.set_xscale('log')
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.tight_layout()
            print(f'Saving PSD for percetage {i} and signal {j} \n')
            print(40*'-')
            plt.savefig(os.path.join(path, f'images/psd/{folder_name}', filename + f'-{i}-s{j}-PSD.png'), dpi=600)
            plt.close()
            
            k += 1

def wavelet_generation(pr_i:int,pr_n:int,folder_name:str,jump:int,offset:int,y_lim:list,chanels:list,path:str,filename:str):
    
    plt.rcParams['agg.path.chunksize'] = 10000
    
    for k in chanels: 
        os.makedirs(os.path.join(path, f'images/wavelet/{folder_name}'), exist_ok=True)
    
    filename = filename.split('-')[0]
    pr_used = np.arange(pr_i,pr_n+1,jump)
    
    for i in pr_used:
    
        k=offset
        for j in chanels:
            
            data = h5py.File(os.path.join(path,filename + f'-{i}.h5'),'r')
            time_data = data['Table1']['Ds001-Time'][:]
            pressure_data = data['Table1'][f'Ds{k:03d}-Signal {j}'][:]
        
            dt=time_data[1]-time_data[0]
            fs=1.0/dt

            n_scales = 300
            scales = np.arange(1, n_scales)
            coeffs, freqs = pywt.cwt(pressure_data, scales, 'cmor', sampling_period=dt) 
            
            power = np.abs(coeffs)**2
            power_norm = power / np.var(pressure_data)
            extent = [time_data.min(), time_data.max(), freqs.min(), freqs.max()]

            plt.imshow(10*np.log10(power_norm + 1e-20), extent=extent, cmap='hot_r', aspect='auto', origin='lower',vmin=-25,vmax=0.1)
            plt.colorbar(label='Normalized $|W|^2$')
            plt.xlabel('Time [s]',fontsize=14, style='italic')
            plt.ylabel('Frequency [Hz]',fontsize=14, style='italic')
            plt.ylim([100,1500])
            plt.yscale('log')
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.tight_layout()
            print(f'Saving spectrogram for data with iterating {i}-{j}')
            print(40*'-')
            plt.savefig(os.path.join(path, f'images/wavelet/{folder_name}', filename + f'-{i}-s{j}-Wavelet.png'), dpi=600)
            plt.close()

            k += 1


def spectrogram_generation(pr_i:int,pr_n:int,folder_name:str,jump:int,offset:int,y_lim:list,chanels:list,path:str,filename:str):
    
    plt.rcParams['agg.path.chunksize'] = 10000
    
    for k in chanels: 
        os.makedirs(os.path.join(path, f'images/spectrogram/{folder_name}'), exist_ok=True)
    
    filename = filename.split('-')[0]
    pr_used = np.arange(pr_i,pr_n+1,jump)
    
    for i in pr_used:
    
        k=offset
        for j in chanels:
            
            data = h5py.File(os.path.join(path,filename + f'-{i}.h5'),'r')
            time_data = data['Table1']['Ds001-Time'][:]
            pressure_data = data['Table1'][f'Ds{k:03d}-Signal {j}'][:]
        
            dt=time_data[1]-time_data[0]
            fs=1.0/dt
            
            fmin = 10
            n_chunk = 4
            lensg_exp = pressure_data.size
            nperseg_exp = lensg_exp/n_chunk
            nfft_exp = next_greater_power_of_2(int(nperseg_exp))
            noverlap_exp = nperseg_exp/2

            if nperseg_exp > lensg_exp:
                raise RuntimeError('Wrong value for $f_{min}$')

            f, t, pxx = sg.spectrogram(pressure_data,fs=fs,window='hann',nperseg=int(nperseg_exp),nfft=nfft_exp,scaling='density')

            cmap = plt.get_cmap('hot_r')
            fig = plt.pcolormesh(t, f, 10*np.log10(pxx/4.0e-10), shading='gouraud',cmap=cmap,vmin=-10,vmax=30)
            plt.ylim([100,15000])
            ax = plt.gca()
            ax.set_yscale('log')
            ax.yaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            plt.ylabel('Frequency $[Hz]$',fontsize=14, style='italic')
            plt.xlabel('Time $(s)$',fontsize=14, style='italic')
            cbar = plt.colorbar()
            cbar.set_label('PSD $[dB/Hz]$',fontsize=14, fontstyle='italic')
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.tight_layout()
            print(f'Saving spectrogram for data with iterating {i}-{j}')
            print(40*'-')
            plt.savefig(os.path.join(path, f'images/spectrogram/{folder_name}', filename + f'-{i}-s{j}-Spectrogram.png'), dpi=600)
            plt.close()               
           
            k += 1    
        
def directivity_generation(pr_i:int,pr_n:int,folder_name:str,jump:int,offset:int,ylim:list,chanels:list,polar_angles:list,path:str,filename:str,case:str,*targets:list):
    
    '''
    This function is used to generate the directivity plot for the different cases.
    It is not used in the main script, but it can be useful for future analysis.
    '''

    plt.rcParams['agg.path.chunksize'] = 10000
    
    for k in chanels: 
        os.makedirs(os.path.join(path, f'images/directivity/{folder_name}'), exist_ok=True)
        
    oaspl_data = np.zeros(len(chanels))
    filename = filename.split('-')[0]
    pr_used = np.arange(pr_i,pr_n+1,jump)

    for z in range(int(len(targets[0])/2)):
        target_i = targets[0][z*2]
        target_n = targets[0][z*2+1]
        
        for i in pr_used:
        
            k=offset
            for j in chanels:
                
                data = h5py.File(os.path.join(path,filename + f'-{i}.h5'),'r')
                time_data = data['Table1']['Ds01-Time'][:]
                pressure_data = data['Table1'][f'Ds{k:02d}-Signal {j}'][:]
            
                dt=time_data[1]-time_data[0]
                fs=1.0/dt
                
                fmin = 10
                n_chunk = 20
                lensg_exp = pressure_data.size
                nperseg_exp = lensg_exp/n_chunk
                nfft_exp = next_greater_power_of_2(int(nperseg_exp))
                noverlap_exp = nperseg_exp/2

                if nperseg_exp > lensg_exp:
                    raise RuntimeError('Wrong value for $f_{min}$')

                [f,Pxx]=sg.welch(pressure_data,fs=fs,window='hann',nperseg=nperseg_exp,nfft=nfft_exp,scaling='density')
                f = f[target_i:target_n]
                Pxx = Pxx[target_i:target_n]
                if case == 'integrated':
                    oaspl_data[k-offset] = 10*np.log10(np.trapz(Pxx,f)/4.0e-10)
                elif case == 'max':
                    oaspl_data[k-offset] = np.max(10*np.log10(Pxx/4.0e-10))
                    
                k += 1
            plt.polar(polar_angles*np.pi/180,oaspl_data,'o',color='black')
            plt.grid(True, which='both', ls='--')
            plt.ylim(ylim)
            ax=plt.gca()
            plt.tight_layout()
            print(f'Saving directivity for percetage {i} and signal {j} \n')
            print(40*'-')
            plt.savefig(os.path.join(path, f'images/directivity/{folder_name}', filename + f'-{i}-directivity-oaspl.png'), dpi=600)    
            #plt.show()
            plt.close()
              
                    
def ladder_psd_generation(pr_i:int,pr_n:int,aoa:int,jump:int,chanels:list,path:str,filename:str):
    
    '''
    This function is used to generate the PSD for the ladder plot.
    It is not used in the main script, but it can be useful for future analysis.
    '''

    plt.rcParams['agg.path.chunksize'] = 10000
    
    os.makedirs(os.path.join('/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/Others/Mine/2025-ISAE-Exp/NACA0015', f'images/ladder_psd/{aoa}deg'), exist_ok=True)
    
    filename = filename.split('-')[0]
    
    k=2
    for j in chanels:
        pr_used = np.arange(pr_i,pr_n+1,jump)
        calculated_psd = []
        for i in pr_used:
            
            data = h5py.File(os.path.join(path,filename + f'-{i}.h5'),'r')
            time_data = data['Table1']['Ds1-Time'][:]
            pressure_data = data['Table1'][f'Ds{k}-Signal {j}'][:]
        
            dt=time_data[1]-time_data[0]
            fs=1.0/dt
            
            fmin = 10
            n_chunk = 20
            lensg_exp = pressure_data.size
            nperseg_exp = lensg_exp/n_chunk
            nfft_exp = next_greater_power_of_2(int(nperseg_exp))
            noverlap_exp = nperseg_exp/2

            if nperseg_exp > lensg_exp:
                raise RuntimeError('Wrong value for $f_{min}$')

            [f,Pxx]=sg.welch(pressure_data,fs=fs,window='hann',nperseg=nperseg_exp,nfft=nfft_exp,scaling='density')
            
            calculated_psd.append(10*np.log10(Pxx/4.0e-10))
        
        k += 1   
        Z = np.array(calculated_psd)
        f = np.array(f)
        pr_used, f = np.meshgrid(pr_used, f)

        plt.figure()
        plt.pcolormesh(pr_used, f, Z.T, shading='gouraud',cmap='grey_r', vmin=-20, vmax=45)
        plt.colorbar(label='PSD $\\left[\\frac{dB}{Hz}\\right]$')
        plt.grid(True, which='both', ls='--')
        plt.xlabel('Measurement N°$[]$',fontsize=12, style='italic')
        plt.ylabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
        plt.ylim([70,2000])
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        print(f'Saving ladder PSD for percetage {i} and signal {j} \n')
        print(40*'-')
        plt.savefig(os.path.join('/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/Others/Mine/2025-ISAE-Exp/NACA0015', f'images/ladder_psd/{aoa}deg', filename + f'-{i}-s{j}-ladder-PSD.png'), dpi=600)
        plt.close()
            
def cp_generation(pr_i:int,pr_n:int,aoa:int,slope:float,intercept:float,offset:int,jump:int,chord_coord:list,path:str,*filenames_chanels:list):
    
    '''
    This function is used to generate the Cp for the ladder plot.
    It is not used in the main script, but it can be useful for future analysis.
    '''

    plt.rcParams['agg.path.chunksize'] = 10000
    pref = 101325.0
    Tref = 20.0
    rho_2 = pref/(287.01*(273.15+Tref))
    
    filenames_chanels = filenames_chanels[0]
    sv_used = int(len(filenames_chanels)/2)
    filenames = filenames_chanels[0:sv_used]
    chanels = filenames_chanels[sv_used:]
    
    os.makedirs(os.path.join('/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/Others/Mine/2025-ISAE-Exp/NACA0015', f'images/cp/{aoa}deg'), exist_ok=True)

    pr_used = np.arange(pr_i,pr_n+1,jump)

    cp_data = {}
    for j in pr_used: 
        cp_data.clear()
        velocity = slope*j+offset + intercept
        pdyn_2 = (rho_2*velocity**2)/2
        for i in range(sv_used):
            filenames[i] = filenames[i].split('-')[0]
            cp_array = np.loadtxt(os.path.join(path, filenames[i] + f'-{j}.csv'), delimiter=',')[:,1:17].mean(axis=0)
            cp_array = cp_array[0:chanels[i]]/pdyn_2
            
            cp_data[(i, j)] = cp_array
    

        cps = [cp_data[(i, j)] for i in range(sv_used)]
        cps_concatenated = np.concatenate(cps)

        plt.figure(1)
        # plt.plot(Cp_exp[:,0],Cp_exp[:,1],"o",color='red',label='Reference')
        plt.plot(chord_coord,cps_concatenated,"o",color='black',label=str(np.round(velocity,decimals=2))+' m/s')
        # plt.plot(xC_scan,-Cp_scan_ref,linestyle='None',marker='o',label='Exp UTIAS $15^\\circ$ - 16 m/s')
        #plt.plot(xC_scan_2,-Cp_scan,linestyle='None',marker='o',color='k',label=f'Exp ISAE ${AoA}^\\circ$ - {velocity_name} m/s')
        plt.legend()
        plt.grid(True, which='both', ls='--')
        plt.xlabel('$x/C$',fontsize=16)
        plt.ylabel('$-C_p$',fontsize=16)
        plt.ylim([-2,2])
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        print(f'Saving Cp for percetage {j}\n')
        print(40*'-')
        plt.savefig(os.path.join('/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/Others/Mine/2025-ISAE-Exp/NACA0015',f'images/cp/{aoa}deg',f'cp-NACA0015-{aoa}deg-{j}pr') + '.png' , dpi=600)
        # plt.show()
        plt.close(1)

    # scan_1 = np.loadtxt(os.path.join(path,filename_1),delimiter=',')[:,1:17].mean(axis=0)
    
    
