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
parser.add_argument("--output-directory", type = str, default = '.')
parser.add_argument("--cut-size", type = str, default=None)
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

    ifar_ref = ref['found_after_vetoes/ifar'][:][loc_ref]
    injid_ref = injid_ref[loc_ref]
    injid_comp = injid_comp[loc_comp]


    #Cutting the data based on IFAR
    if opt.cut_size == 'small':
        cut_ref = numpy.where(ifar_ref > 1.0)[0]
    else:
        cut_ref = numpy.logical_and(ifar_ref > 1.0, ifar_ref < 1000)
    cut_comp = numpy.where(numpy.in1d(injid_comp, injid_ref[cut_ref]))[0]
    loc_comp = loc_comp[cut_comp]
    loc_ref = loc_ref[cut_ref]
    ifar_ref = ifar_ref[cut_ref]
    injid_ref = injid_ref[cut_ref]
    injid_comp = injid_comp[cut_comp]

    tid2_comp = comp['found_after_vetoes/trigger_id1'][:][loc_comp]
    tid1_comp = comp['found_after_vetoes/trigger_id2'][:][loc_comp]
    tid2_ref = ref['found_after_vetoes/trigger_id1'][:][loc_ref]
    tid1_ref = ref['found_after_vetoes/trigger_id2'][:][loc_ref]
    stat_comp_calc = comp['found_after_vetoes/stat'][:][loc_comp]
    stat_ref_calc = ref['found_after_vetoes/stat'][:][loc_ref]
    ifar_comp = comp['found_after_vetoes/ifar'][:][loc_comp]

    ifs_ref = numpy.append(ifs_ref, ifar_ref)
    ifs_comp = numpy.append(ifs_comp, ifar_comp)
    stat_comp = numpy.append(stat_comp, stat_comp_calc)
    stat_ref = numpy.append(stat_ref, stat_ref_calc)

    #Chirp Mass Calculation
    m1_comp = numpy.append(m1_comp, comp['injections/mass1'][:][injid_comp])
    m2_comp = numpy.append(m2_comp, comp['injections/mass2'][:][injid_comp])
    m_total_comp = m1_comp + m2_comp
    eta_comp = (m1_comp * m2_comp) / (m_total_comp * m_total_comp)
    chirp_comp = m_total_comp * (eta_comp)**(3.0/5.0)
 
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

combined_snr_ratio = (((H1_snr_comp)**2.0 + (L1_snr_comp)**2.0)**0.5) / (((H1_snr_ref)**2.0 + (L1_snr_ref)**2.0)**0.5)
stat_ratio = stat_comp / stat_ref
chsq_H1_ratio = chsq_H1_comp / chsq_H1_ref
chsq_L1_ratio = chsq_L1_comp / chsq_L1_ref
snr_H1_ratio = H1_snr_comp / H1_snr_ref
snr_L1_ratio = L1_snr_comp / L1_snr_ref
ifar_ratio = ifar_comp / ifar_ref

############### PLOTTING TIME ###############
params = [chsq_H1_ratio, chsq_L1_ratio, combined_snr_ratio]
param_name = ['H1 Red. Chisq Ratio', 'L1 Red. Chisq Ratio', 'Network SNR Ratio']

for param, name in zip(params, param_name):

    fig, ax = pylab.subplots(1, 1, figsize=[15,10])
    ax.hexbin(param, stat_ratio)
    ax.set_ylabel('New SNR Ratio')
    ax.set_xlabel(name)
    ax.set_title('New SNR Ratio ' + name)
    fig.colorbar(ax.hexbin(param, stat_ratio), ax=ax)
#    fig.savefig(str(opt.output_directory)  + '/' + 'New SNR Ratio vs. ' + name + ' lin, hex.png')
    plt.close()

    fig, ax = pylab.subplots(1, 1, figsize=[15,10])
    ax.hexbin(param, stat_ratio, yscale='log')
    ax.set_ylabel('New SNR Ratio')
    ax.set_xlabel(name)
    ax.set_title('New SNR Ratio ' + name)
    fig.colorbar(ax.hexbin(param, stat_ratio, yscale='log'), ax=ax)
#    fig.savefig(str(opt.output_directory)  + '/' + 'New SNR Ratio vs. ' + name + ' log, hex.png')
    plt.close()

fig, ax = pylab.subplots(1, 1, figsize=[15,10])
ax.hist(combined_snr_ratio, bins = 20)
ax.set_xlabel('Network SNR Ratio (512s / 256s)')
fig.savefig(str(opt.output_directory)  + '/' + 'Network SNR Ratio, hist.png')
plt.close()

fig, ax = pylab.subplots(1, 1, figsize=[15,10])
ax.hist(stat_ratio, bins = 20)
ax.set_xlabel('New SNR Ratio (512s / 256s)')
fig.savefig(str(opt.output_directory)  + '/' + 'New SNR Ratio, hist.png')
plt.close()
print stat_ratio[numpy.where(stat_ratio > 1.05)[0]]

fig, ax = pylab.subplots(1, 1, figsize=[15,10])
ax.hist(ifar_ratio, bins = 20)
ax.set_xlabel('IFAR Ratio (512s / 256s)')
fig.savefig(str(opt.output_directory)  + '/' + 'IFAR Ratio, hist.png')
plt.close()

fig, ax = pylab.subplots(1, 1, figsize=[15,10])
ax.scatter(time_comp, combined_snr_ratio, c=ifs_ref, norm=matplotlib.colors.LogNorm())
ax.set_xlabel('Time')
ax.set_ylabel('Network SNR Ratio')
ax.set_yscale('log')
#fig.savefig(str(opt.output_directory)  + '/' + 'Network SNR vs. Time, scatter, log.png')
plt.close()

fig, ax = pylab.subplots(1, 1, figsize=[15,10])
ax.scatter(time_comp, combined_snr_ratio, c=ifs_ref, norm=matplotlib.colors.LogNorm())
ax.set_xlabel('Time')
ax.set_ylabel('Network SNR Ratio')
#fig.savefig(str(opt.output_directory)  + '/' + 'Network SNR vs. Time, scatter, lin.png')
plt.close()
