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
opt = parser.parse_args()

f = open('output.txt', 'w')
sys.stdout = f

fols = glob.glob(str(opt.comparison_directory)+'/*INJ_coinc')
#print 'fols: ', fols
#print 
pd, td, rd, md, ifs, ifs2 = [], [], [], [], [], []
time = []
stat_diff = []
inc = []

o1, o2 = [], []
m1, m2 = [], []
chisq = []

g1, g2 = [], []
v1, v2 = [], []

pv = []
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

    injid = comp['found_after_vetoes/injection_index'][:]
    injid2 = ref['found_after_vetoes/injection_index'][:]
    
    loc = numpy.where(numpy.in1d(injid, injid2))[0]
    loc2 = numpy.where(numpy.in1d(injid2, injid))[0]
    
    missed_ref = len(injid) - numpy.sum( numpy.in1d(injid, injid2) )
    total_missed += missed_ref

    tid2 = comp['found_after_vetoes/trigger_id1'][:][loc]
    tid1 = comp['found_after_vetoes/trigger_id2'][:][loc]
    ifar = comp['found_after_vetoes/ifar'][:][loc]
    ifar2 = ref['found_after_vetoes/ifar'][:][loc2]
    stat = comp['found_after_vetoes/stat'][:][loc]
    stat2 = ref['found_after_vetoes/stat'][:][loc2]
    injid = injid[loc]
    injid2 = injid2[loc2]

    print 'Number of instances where IFAR in reference is greater than that in comparison: ' , (ifar2 > ifar).sum()
    print 'Number of instances where IFAR in comparison is greater than that in reference: ' , (ifar > ifar2).sum()
    print 'Number of instances where IFAR is equal in both: ' , (ifar == ifar2).sum()
    print 'Number of injections found after vetoes in comparison: ' , len(ifar)
    print 'Number of injections found in reference, but missed in comparison: ', missed_ref
    
    v = numpy.where(ifar2 > 1.0)[0]
    tid2 = tid2[v]
    tid1 = tid1[v]
    loc = loc[v]
    loc2 = loc2[v]
    ifar = ifar[v]
    stat = stat[v]
    stat2 = stat2[v]
    ifar2 = ifar2[v]
    injid = injid[v]
    ifs2 = numpy.append(ifs2, ifar2) 
    ifs = numpy.append(ifs, ifar)

    print "Number found with ifar2 > 1.0: {0}".format( len(ifar) )
    total_found += len(ifar)

    
   # inc = numpy.append(inc, s['injections/inclination'][:][injid])
   # m1 = numpy.append(m1, s['injections/mass1'][:][injid])
   # m2 = numpy.append(m2, s['injections/mass2'][:][injid])

    print 'len ifs2: ', len(ifs2),'len ifs: ',  len(ifs)    
    print '__________________________________________________________________________________________________________________________'
    time = numpy.append(time, ref['found_after_vetoes/time1'][:][loc])
    c1 = t1['H1/chisq'][:][tid1] / (2 * t1['H1/chisq_dof'][:][tid1] - 2)
    c2 = t2['L1/chisq'][:][tid2] / (2 * t2['L1/chisq_dof'][:][tid2] - 2)
    chisq = numpy.append(chisq, numpy.maximum(c1, c2))
  
    o2 = numpy.append(o2, comp['injections/optimal_snr_1'][:][injid])
    o1 = numpy.append(o1, comp['injections/optimal_snr_2'][:][injid])

    gs1 = t1['H1/sigmasq'][:][tid1] ** 0.5
    gs2 = t2['L1/sigmasq'][:][tid2] ** 0.5

    vs1 = t1['H1/snr'][:][tid1]
    vs2 = t2['L1/snr'][:][tid2]

    #pd = numpy.append(pd, (t1['H1/coa_phase'][:][tid1] - t2['L1/coa_phase'][:][tid2]) % (numpy.pi * 2.0))
    #td = numpy.append(td, (t1['H1/end_time'][:][tid1] - t2['L1/end_time'][:][tid2]))
     
    rd = numpy.append(rd, ifar/ifar2) 
    md = numpy.append(md, (vs1**2.0 + vs2**2.0) ** 0.5)

    g1 = numpy.append(g1, gs1 / vs1)
    g2 = numpy.append(g2, gs2 / vs2)
    v1 = numpy.append(v1, vs1)
    v2 = numpy.append(v2, vs2)
    pv = numpy.append(pv, stat**2.0 - newsnr(vs1, c1)**2.0 - newsnr(vs2, c2)**2.0) 

    stat_diff = numpy.append(stat_diff, abs(stat2 - stat))
    #print pv
    #print v1[0:3], v2[0:3], o1[0:3], o2[0:3]


print "SIZE", len(o1)
print "Total number found with reference IFAR >1 and comparison IFAR <1000: {0}".format(total_found)
print "Total number found (at all) in reference but missed in comparison: {0}".format(total_missed)

df = ifs / ifs2
s = df.argsort()[::-1]
#o1 = o1[s]
#o2 = o2[s]
df = df[s]
rd = rd[s]
#md = md[s]
time = time[s]
#chisq = chisq[s]
ifs = ifs[s]
ifs2 = ifs2[s]
g1 = g1[s]
g2 = g2[s]
v1 = v1[s]
v2 = v2[s]
#pv = pv[s]
print "SIZE, sort", len(o1)

#Creating Plots
fig, axs = pylab.subplots(2, 2,figsize=[15,15])

ax = axs[0,0]
ax.scatter(ifs2, stat_diff,c=ifs2, norm=matplotlib.colors.LogNorm())
ax.set_ylabel('abs(difference)')
ax.set_xlabel('IFAR, Reference')
ax.set_xlim([0,6000])
ax.set_title('Statistic Difference vs. IFAR, Reference')

ax = axs[0,1]
ax.scatter( 4.0 * 256.0 * numpy.minimum((v1/g1),(v2/g2)), stat_diff, c=ifs2, norm=matplotlib.colors.LogNorm())
ax.set_xlabel('Minimum New SNR')
ax.set_ylabel('abs(difference)')
ax.set_title('Statistic Difference vs. Min of New SNR')

ax = axs[1,0]
ax.scatter(ifs2, rd, c=ifs2, norm=matplotlib.colors.LogNorm())
ax.set_xlabel('IFAR, Reference')
ax.set_ylabel('IFAR Comparison / IFAR Reference')
ax.set_title('IFAR Ratio  vs. IFAR, reference')

ax = axs[1,1]
points = ax.scatter( 4.0 * 256.0 * numpy.minimum((v1/g1),(v2/g2)), rd, c=ifs2,  norm=matplotlib.colors.LogNorm())
ax.set_xlabel('Minimum New SNR')
ax.set_ylabel('IFAR Comparison / IFAR Reference')
ax.set_title('IFAR Ratio vs. Minimum New SNR')

#Links Plots
import mpld3, mpld3.plugins
mpld3.plugins.connect(fig, mpld3.plugins.LinkedBrush(points))
#pylab.savefig("test.png")
mpld3.save_html(fig, str(opt.html_output)) #Make an argument for where to save it.

f.close()
