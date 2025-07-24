import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.io import savemat
from scipy import stats
import h5py
import scipy.signal as sg
from rich import print as rprint

def next_greater_power_of_2(x):
    return round(2**(x-1).bit_length())

def time_trace_generation(pr_i:int,pr_n:int,aoa:int,chanels:list,path:str,filename:str):
    
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
    pr_used = np.arange(pr_i,pr_n+1,1)
    
    for i in pr_used:
        
        k=2
        for j in chanels:
            
            data = h5py.File(os.path.join(path,filename + f'-{i}.h5'),'r')
            time_data = data['Table1']['Ds1-Time'][:]
            pressure_data = data['Table1'][f'Ds{k}-Signal {j}'][:]
            
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
            

def psd_generation(pr_i:int,pr_n:int,aoa:int,y_lim:list,chanels:list,path:str,filename:str):
    
    plt.rcParams['agg.path.chunksize'] = 10000
    
    for k in chanels: 
        os.makedirs(os.path.join(path, f'images/psd/{aoa}deg/signal_{k}'), exist_ok=True)
    
    filename = filename.split('-')[0]
    pr_used = np.arange(pr_i,pr_n+1,1)
    
    for i in pr_used:
    
        k=2
        for j in chanels:
            
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
            plt.savefig(os.path.join(path, f'images/psd/{aoa}deg/signal_{j}', filename + f'-{i}-s{j}-PSD.png'), dpi=600)
            plt.close()
            
            k += 1