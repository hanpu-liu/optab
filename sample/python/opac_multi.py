import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py
import argparse

# Constants
K_BOL = 1.3806503e-16

def read_mono_hdf5(file_path):
    try:
        with h5py.File(file_path, 'r') as f:
            return {name: f[name][:] for name in f.keys()}
    except Exception as e:
        raise RuntimeError(f"Failed to read {file_path}: {e}")

def opac(dir_path, mean, savedata, syms):
    files = [os.path.join(dir_path, 'output', f) for f in os.listdir(os.path.join(dir_path, 'output')) if f.startswith("mono_") and f.endswith(".h5")]
    nmax = len(files)
    
    if nmax == 0:
        raise ValueError("No HDF5 files found.")

    print(f"Reading {nmax} layers ... ")

    datasets = [read_mono_hdf5(file) for file in files]
    # tmp2 = np.array([data['temp2'] for data in datasets])
    # tmp_tot = np.array([data['temp'][0] for data in datasets])
    # rho_tot = np.array([data['rho'][0] for data in datasets])
    # nden_tot = np.array([data['nden'][0] for data in datasets])
    # ros_tot = np.log10(np.array([data['ros'] for data in datasets]) / rho_tot)
    # pla_tot = np.log10(np.array([data['plac'] + data['plal'] for data in datasets]) / rho_tot)
    # pla2_tot = np.log10(np.array([data['plac2'] + data['plal2'] for data in datasets]) / rho_tot)

    df = pd.DataFrame({'temp': [data['temp'][0] for data in datasets], 
                       'rho': [data['rho'][0] for data in datasets], 
                       'nden': [data['nden'][0] for data in datasets], 
                       'logkappa_r': [np.log10(data['ros'] / data['rho'][0]) for data in datasets], 
                       'logkappa_p': [np.log10((data['plac'] + data['plal']) / data['rho'][0]) for data in datasets],
                       'logkappa_s': [np.log10(data['scamean'] / data['rho'][0]) for data in datasets],
                       'logkappa_e': [np.log10(data['eff'] / data['rho'][0]) for data in datasets],
                       'logkappa_a': [np.log10(data['absmean'] / data['rho'][0]) for data in datasets]})
    df = df.sort_values(by=['temp', 'rho'], ignore_index=True)

    if savedata:
        with h5py.File(os.path.join(dir_path, 'output.h5'), 'w') as f:
            f.create_dataset('temp', data=df['temp'])
            f.create_dataset('rho', data=df['rho'])
            f.create_dataset('nden', data=df['nden'])
            f.create_dataset('logkappa_r', data=np.array(df['logkappa_r'].to_list()))
            f.create_dataset('logkappa_p', data=np.array(df['logkappa_p'].to_list()))
            f.create_dataset('logkappa_s', data=np.array(df['logkappa_s'].to_list()))
            f.create_dataset('logkappa_e', data=np.array(df['logkappa_e'].to_list()))
            f.create_dataset('logkappa_a', data=np.array(df['logkappa_a'].to_list()))
        print(f"Data saved to {os.path.join(dir_path, 'output.h5')}")

    # save the data
    # if savedata:
    #     with h5py.File(os.path.join(dir_path, 'output.h5'), 'w') as f:
    #         f.create_dataset('temp', data=tmp_tot)
    #         f.create_dataset('rho', data=rho_tot)
    #         f.create_dataset('nden', data=nden_tot)
    #         f.create_dataset('logkappa_r', data=ros_tot)
    #         f.create_dataset('logkappa_p', data=pla_tot)
    #     print(f"Data saved to {os.path.join(dir_path, 'output.h5')}")

    # # Color scaling for visualization
    # vmax, vmin = 7.0, -6.0

    # # Title and label mapping
    # titles = {
    #     'ross': ('Rosseland-mean opacity', 'log $\kappa$ [cm$^2$/g]'),
    #     'pla': ('Planck-mean opacity', 'log $\kappa$ [cm$^2$/g]'),
    #     'pla2': (f'Planck-mean opacity at T_rad={int(tmp2[0])}K', 'log $\kappa$ [cm$^2$/g]')
    # }

    # if mean not in titles:
    #     raise ValueError(f"Unknown mean type: {mean}")
    # title, btitle = titles[mean]

    # # Plotting
    # plt.rcParams.update({'font.size': 16})  # Adjust the number to your preference
    
    # plt.figure(figsize=(10, 8))
    # scatter = plt.scatter(tmp_tot, rho_tot, c=
    #                       ros_tot[-1] if mean == 'ross' else
    #                       pla_tot[-1] if mean == 'pla' else
    #                       pla2_tot[-1],
    #                       s=syms, cmap='jet', norm=plt.Normalize(vmin=vmin, vmax=vmax), marker='s')
    # plt.yscale('log')
    # plt.xscale('log')
    # plt.xlabel('T [K]')
    # plt.ylabel(r'$\rho$ [g/cm$^3$]')
    # plt.title(title)
    # plt.colorbar(scatter, label=btitle)
    # plt.grid(True)

    # plt.show()

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='Process some integers.')
    # Add the arguments
    parser.add_argument('dir_path', type=str, help='The path to the directory containing HDF5 files')
    parser.add_argument('mean', type=str, choices=['ross', 'pla', 'pla2'], help='Type of mean opacity to plot')
    parser.add_argument('savedata', type=bool, help='Save the data to HDF5 file', default=False)
    parser.add_argument('syms', type=int, help='Marker size for scatter plot', default=100)

    # Execute the parse_args() method
    args = parser.parse_args()

    opac(args.dir_path, args.mean, args.savedata, args.syms)

if __name__ == '__main__':
    main()
