import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

### Data
paths = ["" for x in range(41)]
for i in range(41):
	paths[i] = "SiPM_refl_%s_z21"%(900-i)

data  = np.zeros(( len(paths),37 ))
dataS = np.zeros(( len(paths),37 ))
sort_temp = np.zeros(( 37 ))
for i in range( len(paths) ):
	load = np.loadtxt(paths[i])
	for j in range( 37 ):
		data[i][j]   = load[j][5]/5.
		sort_temp[j] = data[i][j]
	sort_temp = sorted(sort_temp)
	for j in range( 37 ):
		dataS[i][j] = sort_temp[j]

Angle = np.zeros(( len(paths) ))
for i in range( len(paths) ):
	Angle[i] = 0.45 * i

### Calculate
Max = np.zeros(( len(paths),37 ))
for i in range( len(paths) ):
	for j in range(37):
		for k in range(j+1):
			Max[i][j] = Max[i][j] + dataS[i][36-k]
		Max[i][j] = Max[i][j]/(j+1)

### Write
target = open("SiPM_refl_zero", 'w')
target.truncate()

for i in range( len(paths) ):
	target.write( "%f	%f \n" % (Angle[i], Max[i][2]) )

target.close()

### Plot peaks
ZeroPlot = pl.figure(num=None, figsize=(14, 7), dpi=80, facecolor='w', edgecolor='k')

PlotColor	= [ '#2f4343', '#e30700', '#4387e3', '#47af27', '#ff6f17', '#db3b8f', '#00dddd', '#d7a713' ]

for j in range(30):
	PMTrate	= np.array([ Max[i][j] for i in range( len(paths) ) ])
	pl.plot( Angle, PMTrate, linestyle='None', marker='x', color=PlotColor[j%8], label='No. of Max: %s'%(j+1) )

### Format plot
pl.title("SiPM")
pl.xlabel("PMT-angle [deg]")
pl.ylabel("PMT rate [Hz]")
pl.xlim(3,20)
pl.ylim(0,16000)
pl.xticks(np.arange(5, 20, 5))

pl.legend(loc='upper right', fontsize = 'small', ncol=2, handleheight=2., labelspacing=0.01)

pl.grid()

### Save & Show Plot
pl.savefig("Zero_peaks.pdf")
pl.show()





