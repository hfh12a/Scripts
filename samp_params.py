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
parser.add_argument("--output-kind", type=str, default=None)
opt = parser.parse_args()

matplotlib.rcParams.update({'font.size': 18})


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
    missed_comp = numpy.where(numpy.in1d(injid_ref, injid_comp, invert = True))[0]

    tid2_ref = ref['found_after_vetoes/trigger_id1'][:][missed_comp]
    tid1_ref = ref['found_after_vetoes/trigger_id2'][:][missed_comp]
    ifar_ref = ref['found_after_vetoes/ifar'][:][missed_comp]
    injid_ref = injid_ref[missed_comp]

    #Cutting the data based on IFAR
    cut = numpy.where(ifar_ref > 1.0) 
    tid2_ref = tid2_ref[cut]
    tid1_ref = tid1_ref[cut]
    injid_ref = injid_ref[cut]
    ifar_ref = ifar_ref[cut]
    missed_comp = missed_comp[cut]
   
    ifar_ref = ref['found_after_vetoes/ifar'][:][missed_comp]
    stat_ref = ref['found_after_vetoes/stat'][:][missed_comp]
    fap_ref = ref['found_after_vetoes/fap'][:][missed_comp]
    exc_fap_ref = ref['found_after_vetoes/fap_exc'][:][missed_comp]
    exc_ifar_ref = ref['found_after_vetoes/ifar_exc'][:][missed_comp]
       
    #Reduced Chi Squared Calculations
    chsq_H1_ref= t_ref_H1['H1/chisq'][:][tid1_ref] / (2 * t_ref_H1['H1/chisq_dof'][:][tid1_ref] - 2)
    chsq_L1_ref= t_ref_L1['L1/chisq'][:][tid2_ref] / (2 * t_ref_L1['L1/chisq_dof'][:][tid2_ref] - 2)
   
    #Detector SNRs
    H1_snr_ref = t_ref_H1['H1/snr'][:][tid1_ref]
    L1_snr_ref = t_ref_L1['L1/snr'][:][tid2_ref]
    
    m1_ref = ref['injections/mass1'][:][injid_ref]
    m2_ref = ref['injections/mass2'][:][injid_ref]
    spin1z_ref = ref['injections/spin1z'][:][injid_ref]
    spin2z_ref = ref['injections/spin2z'][:][injid_ref]

    if opt.output_kind == 'table':
        for h in range(0, len(missed_comp)):
            if missed_comp[h]:
                print '|| {0} || {1} || {2} || {3} || {4} || {5} || {6} || {7} || {8} || {9} || {10} || {11} || {12} || {13} ||'.format(ntpath.basename(fol), ifar_ref[h], fap_ref[h], H1_snr_ref[h], L1_snr_ref[h], stat_ref[h], m1_ref[h], m2_ref[h], spin1z_ref[h], spin2z_ref[h], chsq_H1_ref[h], chsq_L1_ref[h], exc_ifar_ref[h], exc_fap_ref[h])
    else:
	for h in range(0, len(missed_comp)):
	    if missed_comp[h]:
	        print 'Injections missed in 2048 Hz run but not 4096 Hz run (With IFAR > 1):'
	        print "Injection: ", ntpath.basename(fol)
	        print "IFAR: ", ifar_ref[h]
	        print "FAP: ", fap_ref[h]
	        print "H1 SNR: ", H1_snr_ref[h]
	        print "L1 SNR: ", L1_snr_ref[h]
	        print "New SNR: ", stat_ref[h]
	        print "Mass 1: ", m1_ref[h]
	        print "Mass 2: ", m2_ref[h]
	        print "Spin 1z: ", spin1z_ref[h]
	        print "Spin 2z: ", spin2z_ref[h]
	        print "H1 Red. Chisq: ", chsq_H1_ref[h]
	        print "L1 Red. Chisq: ", chsq_L1_ref[h]
	        print "Exc. IFAR: ", exc_ifar_ref[h]
	        print "Exc. FAP: ", exc_fap_ref[h]
	        print "__________________________________________________________________________"
	
    H1_snr_ref, L1_snr_ref, H1_snr_comp, L1_snr_comp = [], [], [], []
    chsq_H1_ref, chsq_L1_ref, chsq_H1_comp, chsq_L1_comp = [], [], [], []
    m1_comp, m2_comp, chirp_comp, m1_ref, m2_ref = [], [], [], [], []
    fap_ref, exc_fap_ref, exc_ifar_ref = [], [], []
    ifar_ref, ifar_comp = [], []
    stat_comp, stat_ref = [], []
    spin1z_ref, spin2z_ref = [], []

