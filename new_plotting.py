#!/usr/bin/env python
import h5py, numpy
import argparse
import matplotlib; matplotlib.use('Agg')
import matplotlib.colors
import matplotlib.pyplot as plt
import pylab
import glob
import ntpath
import sys
from pycbc.events import newsnr

parser = argparse.ArgumentParser(usage='')
parser.add_argument("--reference-directory", type=str, default=None)
parser.add_argument("--comparison-directory", type=str, default=None)
opt = parser.parse_args()

matplotlib.rcParams.update({'font.size': 18})

H1_snr_ref, L1_snr_ref, H1_snr_comp, L1_snr_comp = [], [], [], []
chsq_H1_ref, chsq_L1_ref, chsq_H1_comp, chsq_L1_comp = [], [], [], []
m1_comp, m2_comp, chirp_comp = [], [], []
time_comp = []
ifs_ref, ifs_comp = [], []
stat_comp, stat_ref = [], []

fols = glob.glob(str(opt.comparison_directory)+'/*INJ_coinc')
 
for fol in fols:
    #Obtaining relevant files
    comp = h5py.File(glob.glob(fol + '/*INJFIND*.hdf')[0], 'r')
    ref = h5py.File(glob.glob(str(opt.reference_directory) + '/' + ntpath.basename(fol) + '/*INJFIND*.hdf')[0], 'r')
    t_comp_H1 = h5py.File(glob.glob(fol+ '/H1*MERGE*.hdf')[0],'r')
    t_comp_L1 = h5py.File(glob.glob(fol+ '/L1*MERGE*.hdf')[0],'r')
    t_ref_H1 = h5py.File(glob.glob(str(opt.reference_directory) + '/' + ntpath.basename(fol) +'/H1*MERGE*.hdf')[0], 'r')
    t_ref_L1 = h5py.File(glob.glob(str(opt.reference_directory) + '/' + ntpath.basename(fol) +'/L1*MERGE*.hdf')[0], 'r')

    injid_ref = ref['found_after_vetoes/injection_index'][:]
    injid_comp = comp['found_after_vetoes/injection_index'][:]

    loc_comp = numpy.where(numpy.in1d(injid_comp, injid_ref))[0]
    loc_ref = numpy.where(numpy.in1d(injid_ref, injid_comp))[0]

    tid2_comp = comp['found_after_vetoes/trigger_id1'][:][loc_comp]
    tid1_comp = comp['found_after_vetoes/trigger_id2'][:][loc_comp]
    tid2_ref = ref['found_after_vetoes/trigger_id1'][:][loc_ref]
    tid1_ref = ref['found_after_vetoes/trigger_id2'][:][loc_ref]
    ifar_comp = comp['found_after_vetoes/ifar'][:][loc_comp]
    ifar_ref = ref['found_after_vetoes/ifar'][:][loc_ref]
    stat_comp_calc = comp['found_after_vetoes/stat'][:][loc_comp]
    stat_ref_calc = ref['found_after_vetoes/stat'][:][loc_ref]
    injid_comp = injid_comp[loc_comp]
    injid_ref = injid_ref[loc_ref]

    #Cutting the data based on IFAR
    v_ref = numpy.where(ifar_ref > 1.0) 
    v_comp = numpy.where(ifar_comp > 1.0)
    cut = numpy.where(numpy.in1d(injid_ref[v_ref], injid_comp[v_comp])) #fix
    tid2_comp = tid2_comp[cut]
    tid1_comp = tid1_comp[cut]
    tid2_ref = tid2_ref[cut]
    tid1_ref = tid1_ref[cut]
    loc_comp = loc_comp[cut]
    loc_ref = loc_ref[cut]
    ifar_comp = ifar_comp[cut]
    ifar_ref = ifar_ref[cut]
    stat_ref_calc = stat_ref_calc[cut]
    stat_comp_calc = stat_comp_calc[cut]
    injid_ref = injid_ref[cut]
    injid_comp = injid_comp[cut]

    ifs_ref = numpy.append(ifs_ref, ifar_ref)
    ifs_comp = numpy.append(ifs_comp, ifar_comp)
    stat_comp = numpy.append(stat_comp, stat_comp_calc)
    stat_ref = numpy.append(stat_ref, stat_ref_calc)

    #Chirp Mass Calculation
    m1_comp = numpy.append(m1_comp, comp['injections/mass1'][:][injid_comp])
    m2_comp = numpy.append(m2_comp, comp['injections/mass2'][:][injid_comp])
    m_total_comp = m1_comp + m2_comp
    eta_comp = (m1_comp * m2_comp) / (m_total_comp * m_total_comp)
    chirp_comp = numpy.append(chirp_comp, m_total_comp * (eta_comp)**(3.0/5.0))
 
    time_comp = numpy.append(time_comp, comp['found_after_vetoes/time1'][:][loc_comp])
    
    #Reduced Chi Squared Calculations
    chsq_H1_comp_calc = t_comp_H1['H1/chisq'][:][tid1_comp] / (2 * t_comp_H1['H1/chisq_dof'][:][tid1_comp] - 2)
    chsq_L1_comp_calc = t_comp_L1['L1/chisq'][:][tid2_comp] / (2 * t_comp_L1['L1/chisq_dof'][:][tid2_comp] - 2)
    chsq_H1_ref_calc = t_ref_H1['H1/chisq'][:][tid1_ref] / (2 * t_ref_H1['H1/chisq_dof'][:][tid1_ref] - 2)
    chsq_L1_ref_calc = t_ref_L1['L1/chisq'][:][tid2_ref] / (2 * t_ref_L1['L1/chisq_dof'][:][tid2_ref] - 2)
    chsq_H1_comp = numpy.append(chsq_H1_comp, chsq_H1_comp_calc)
    chsq_L1_comp= numpy.append(chsq_L1_comp, chsq_L1_comp_calc)
    chsq_H1_ref = numpy.append(chsq_H1_ref, chsq_H1_ref_calc)
    chsq_L1_ref = numpy.append(chsq_L1_ref, chsq_L1_ref_calc) 
  
    #Network SNR Calculations
    H1_snr_ref = numpy.append(H1_snr_ref, t_ref_H1['H1/snr'][:][tid1_ref])
    L1_snr_ref = numpy.append(L1_snr_ref, t_ref_L1['L1/snr'][:][tid2_ref])
    H1_snr_comp = numpy.append(H1_snr_comp, t_comp_H1['H1/snr'][:][tid1_comp])
    L1_snr_comp = numpy.append(L1_snr_comp, t_comp_L1['L1/snr'][:][tid2_comp])

df = ifs_comp / ifs_ref
s = df.argsort()[::-1]
stat_ref = stat_ref[s]
stat_comp = stat_comp[s]
chirp_comp = chirp_comp[s]
time_comp = time_comp[s]
chsq_H1_ref = chsq_H1_ref[s]
chsq_L1_ref = chsq_L1_ref[s]
chsq_H1_comp = chsq_H1_comp[s]
chsq_L1_comp = chsq_L1_comp[s]
H1_snr_ref = H1_snr_ref[s]
L1_snr_ref = L1_snr_ref[s]
H1_snr_comp = H1_snr_comp[s]
L1_snr_comp = L1_snr_comp[s]

combined_snr_ratio = (((H1_snr_comp)**2.0 + (L1_snr_comp)**2.0)**0.5) / (((H1_snr_ref)**2.0 + (L1_snr_ref)**2.0)**0.5)
stat_ratio = stat_comp / stat_ref
chsq_H1_ratio = chsq_H1_comp / chsq_H1_ref
chsq_L1_ratio = chsq_L1_comp / chsq_L1_ref
snr_H1_ratio = H1_snr_comp / H1_snr_ref
snr_L1_ratio = L1_snr_comp / L1_snr_ref


############### PLOTTING TIME ###############
params = [stat_ratio, chsq_H1_ratio, chsq_L1_ratio, snr_H1_ratio, snr_L1_ratio, combined_snr_ratio]
param_name = ['New SNR Ratio', 'H1 Red. Chisq Ratio', 'L1 Red. Chisq Ratio', 'H1 SNR Ratio', 'L1 SNR Ratio', 'Network SNR Ratio']

for param, name in zip(params, param_name):

    #linear plots vs time
    fig, ax = pylab.subplots(1, 1, figsize=[15,10])
    ax.hexbin(time_comp, param)
    ax.set_xlabel('Time')
    ax.set_ylabel(name)
    ax.set_title(name + ' vs. Time')
    fig.colorbar(ax.hexbin(time_comp, param), ax=ax)
    fig.savefig(name + ' vs. Time, lin.png')
    plt.close()

    #log plots vs time
    fig, ax = pylab.subplots(1, 1, figsize=[15,10])
    ax.hexbin(time_comp, param, yscale = 'log')
    ax.set_xlabel('Time')
    ax.set_ylabel(name)
    ax.set_title(name + ' vs. Time')
    fig.colorbar(ax.hexbin(time_comp, param, yscale = 'log'), ax=ax)
    fig.savefig(name + ' vs. Time, log.png')
    plt.close()
    
    #linear plots vs chirp mass
    fig, ax = pylab.subplots(1, 1, figsize=[15,10])
    ax.hexbin(chirp_comp, param)
    ax.set_xlabel('Chirp Mass')
    ax.set_ylabel(name)
    ax.set_xlim([numpy.min(chirp_comp), numpy.max(chirp_comp)])
    ax.set_title(name + ' vs. Chirp Mass')
    fig.colorbar(ax.hexbin(time_comp, param), ax=ax)
    fig.savefig(name + ' vs. Chirp Mass, lin.png')
    plt.close()
    
    #log plots vs chirp mass
    fig, ax = pylab.subplots(1, 1, figsize=[15,10])
    ax.hexbin(chirp_comp, param, yscale = 'log')
    ax.set_xlabel('Chirp Mass')
    ax.set_ylabel(name)
    ax.set_title(name + ' vs. Chirp Mass')
    ax.set_xlim([numpy.min(chirp_comp), numpy.max(chirp_comp)])
    fig.colorbar(ax.hexbin(time_comp, param, yscale = 'log'), ax=ax)
    fig.savefig(name + ' vs. Chirp Mass, log.png')
    plt.close()
