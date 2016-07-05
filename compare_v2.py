#!/usr/bin/env python
import h5py, numpy
import argparse
import matplotlib; matplotlib.use('Agg')
import matplotlib.colors
import pylab
import glob
import ntpath
import sys
from pycbc.events import newsnr

parser = argparse.ArgumentParser(usage='', description="Comparing outputs from two different runs")
parser.add_argument("--reference-directory", type=str, default=None)
parser.add_argument("--comparison-directory", type=str, default=None)
parser.add_argument("--html-output", type=str, default='test.html')
parser.add_argument("--text-output", type=str, default='output.txt')
opt = parser.parse_args()

f = open(str(opt.text_output), 'w')
sys.stdout = f

fols = glob.glob(str(opt.comparison_directory)+'/*INJ_coinc')
 
rd, md, ifs_comp, ifs_ref, pv = [], [], [], [], []
time = []
stat_diff = []
inc = []

o1, o2 = [], []
m1, m2 = [], []
chisq = []

g1, g2 = [], []
v1, v2 = [], []

total_found = 0
total_missed = 0
for fol in fols:
    print 'Injection: ', ntpath.basename(fol)
    comp = h5py.File(glob.glob(fol + '/*INJFIND*.hdf')[0], 'r')
    ref = h5py.File(glob.glob(str(opt.reference_directory) + '/' + ntpath.basename(fol) + '/*INJFIND*.hdf')[0], 'r')
    print 'Reference File: ' ,glob.glob(str(opt.reference_directory) + '/' + ntpath.basename(fol) + '/*INJFIND*.hdf')[0]
    print 'Comparison File: ', glob.glob(fol + '/*INJFIND*.hdf')[0]
    print 
    
    t1 = h5py.File(glob.glob(fol + '/H1*MERGE*.hdf')[0], 'r')
    t2 = h5py.File(glob.glob(fol + '/L1*MERGE*.hdf')[0], 'r')

    injid_comp = comp['found_after_vetoes/injection_index'][:]
    injid_ref = ref['found_after_vetoes/injection_index'][:]
    
    loc = numpy.where(numpy.in1d(injid_comp, injid_ref))[0]
    loc2 = numpy.where(numpy.in1d(injid_ref, injid_comp))[0]
    
    missed_comp = len(injid_comp) - numpy.sum( numpy.in1d(injid_comp, injid_ref) )
    total_missed += missed_comp

    tid2 = comp['found_after_vetoes/trigger_id1'][:][loc]
    tid1 = comp['found_after_vetoes/trigger_id2'][:][loc]
    ifar_comp = comp['found_after_vetoes/ifar'][:][loc]
    ifar_ref = ref['found_after_vetoes/ifar'][:][loc2]
    stat_comp = comp['found_after_vetoes/stat'][:][loc]
    stat_ref = ref['found_after_vetoes/stat'][:][loc2]
    injid_comp = injid_comp[loc]
    injid_ref = injid_ref[loc2]

    print 'Number of instances where IFAR in reference is greater than that in comparison: ' , (ifar_ref > ifar_comp).sum()
    print 'Number of instances where IFAR in comparison is greater than that in reference: ' , (ifar_comp > ifar_ref).sum()
    print 'Number of instances where IFAR is equal in both: ' , (ifar_comp == ifar_ref).sum()
    print 'Number of injections found after vetoes in comparison: ' , len(ifar_comp)
    print 'Number of injections found in reference, but missed in comparison: ', missed_comp
    
    v = numpy.where(ifar_ref > 1.0)[0]
    tid2 = tid2[v]
    tid1 = tid1[v]
    loc = loc[v]
    loc2 = loc2[v]
    ifar_comp = ifar_comp[v]
    ifar_ref = ifar_ref[v]
    stat_comp = stat_comp[v]
    stat_ref = stat_ref[v]
    injid_comp = injid_comp[v]
    ifs_ref = numpy.append(ifs_ref, ifar_ref) 
    ifs_comp = numpy.append(ifs_comp, ifar_comp)

    print "Number found with ifar_comp > 1.0: {0}".format( len(ifar_comp) )
    total_found += len(ifar_comp)

    #print 'len ifs: ', len(ifs_ref),'len ifs: ',  len(ifs)    
    print '__________________________________________________________________________________________________________________________'
    time = numpy.append(time, ref['found_after_vetoes/time1'][:][loc])
    c1 = t1['H1/chisq'][:][tid1] / (2 * t1['H1/chisq_dof'][:][tid1] - 2)
    c2 = t2['L1/chisq'][:][tid2] / (2 * t2['L1/chisq_dof'][:][tid2] - 2)
    chisq = numpy.append(chisq, numpy.maximum(c1, c2))
  
    o2 = numpy.append(o2, comp['injections/optimal_snr_1'][:][injid_comp])
    o1 = numpy.append(o1, comp['injections/optimal_snr_2'][:][injid_comp])

    gs1 = t1['H1/sigmasq'][:][tid1] ** 0.5
    gs2 = t2['L1/sigmasq'][:][tid2] ** 0.5

    vs1 = t1['H1/snr'][:][tid1]
    vs2 = t2['L1/snr'][:][tid2]
     
    rd = numpy.append(rd, ifar_comp/ifar_ref) 
    md = numpy.append(md, (vs1**2.0 + vs2**2.0) ** 0.5)

    g1 = numpy.append(g1, gs1 / vs1)
    g2 = numpy.append(g2, gs2 / vs2)
    v1 = numpy.append(v1, vs1)
    v2 = numpy.append(v2, vs2)
    pv = numpy.append(pv, stat_comp**2.0 - newsnr(vs1, c1)**2.0 - newsnr(vs2, c2)**2.0) 

    stat_diff = numpy.append(stat_diff, abs(stat_ref - stat_comp))


#print "SIZE", len(o1)
print "Total number found with reference IFAR >1 and comparison IFAR <1000: {0}".format(total_found)
print "Total number found (at all) in reference but missed in comparison: {0}".format(total_missed)

df = ifs_comp / ifs_ref
s = df.argsort()[::-1]
df = df[s]
rd = rd[s]
time = time[s]
ifs_comp = ifs_comp[s]
ifs_ref = ifs_ref[s]
g1 = g1[s]
g2 = g2[s]
v1 = v1[s]
v2 = v2[s]

#print "SIZE, sort", len(o1)

#Creating Plots
fig, axs = pylab.subplots(2, 2,figsize=[15,15])

ax = axs[0,0]
ax.scatter(ifs_ref, stat_diff,c=ifs_ref, norm=matplotlib.colors.LogNorm())
ax.set_ylabel('abs(difference)')
ax.set_xlabel('IFAR, Reference')
ax.set_xlim([0,6000])
ax.set_title('Statistic Difference vs. IFAR, Reference')

ax = axs[0,1]
ax.scatter( 4.0 * 256.0 * numpy.minimum((v1/g1),(v2/g2)), stat_diff, c=ifs_ref, norm=matplotlib.colors.LogNorm())
ax.set_xlabel('Minimum New SNR')
ax.set_ylabel('abs(difference)')
ax.set_title('Statistic Difference vs. Min of New SNR')

ax = axs[1,0]
ax.scatter(ifs_ref, rd, c=ifs_ref, norm=matplotlib.colors.LogNorm())
ax.set_xlabel('IFAR, Reference')
ax.set_ylabel('IFAR Comparison / IFAR Reference')
ax.set_title('IFAR Ratio  vs. IFAR, reference')

ax = axs[1,1]
points = ax.scatter( 4.0 * 256.0 * numpy.minimum((v1/g1),(v2/g2)), rd, c=ifs_ref,  norm=matplotlib.colors.LogNorm())
ax.set_xlabel('Minimum New SNR')
ax.set_ylabel('IFAR Comparison / IFAR Reference')
ax.set_title('IFAR Ratio vs. Minimum New SNR')

#Links Plots
import mpld3, mpld3.plugins
mpld3.plugins.connect(fig, mpld3.plugins.LinkedBrush(points))
#pylab.savefig("test.png")
mpld3.save_html(fig, str(opt.html_output)) #Make an argument for where to save it.

f.close()
