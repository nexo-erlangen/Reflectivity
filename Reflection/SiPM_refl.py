import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import math
from scipy.optimize import curve_fit

################
###   Data   ###
################

# Define Data paths
paths  = ["../SiPM_refl_z2/SiPM_refl_zero",					"SiPM_refl_Xem40_s11",	"SiPM_refl_Xem35_s12",	"SiPM_refl_Xem30_s11",
		  "SiPM_refl_Xem25_s11",	"SiPM_refl_Xem20_s11",	"SiPM_refl_Xem15_s11",	"SiPM_refl_Xem10_s11",  "SiPM_refl_Xem05_s11",
		  "SiPM_refl_Xen00_s12",	"SiPM_refl_Xep05_s11",	"SiPM_refl_Xep10_s11",  "SiPM_refl_Xep15_s11",	"SiPM_refl_Xep20_s11" ]

# Define Xenon angles
XeAngle = [00, -40, -35, -30, -25, -20, -15, -10, -05, +00, +05, +10, +15, +20]

# Initiate arrays
data		= np.zeros(( len(paths),200,2 ))			# 3D data array
DataLowVal	= np.zeros(( len(paths) ))					# Angle of first data line
DataUpBin	= np.zeros(( len(paths) ), dtype=np.int)	# Number of last line

# Fill data arrays
for i in range( len(paths) ):
	load = np.loadtxt(paths[i])											# Load path i
	if i==0:
		for j in range( len(load) ):									# Load zero data (already normalized to Hz)
			data[i][j][0] = load[j][0]									#	Angle 
			data[i][j][1] = load[j][1]									#	Counts
	else:
		for j in range( len(load) ):									# Load reflection data (normalize to Hz!)
			data[i][j][0] = load[j][2]
			data[i][j][1] = load[j][5]/10.
	ytemp = np.array( [ data[i][k][1] for k in range(len(data[i])) ] )	# Fill temporary array
	DataLowVal[i] = data[i][0][0]										# Write lowest angle
	DataUpBin[i]  = np.argmax(ytemp==0)									# Get highest line with count != 0

###############
###   Fit   ###
###############

# Initiate fit parameter array
peak_para = np.zeros(( len(paths),6 ))

# Define fit function
def sinc(x, a, b, c):
	return a * np.sinc((x-b)*c) + 18.

# Fit procedure
for i in range( len(paths) ):
	PMTangle	= np.array([		 data[i][j][0]  for j in range(0,DataUpBin[i]) ])	# New array for PMT angle
	PMTrate		= np.array([		 data[i][j][1]  for j in range(0,DataUpBin[i]) ])	# New array for PMT rate
	PMTrate_err	= np.array([ np.sqrt(data[i][j][1]) for j in range(0,DataUpBin[i]) ])	# New array for error on PMT rate (Poisson)
	init_a = 4000.					# Set fit amplitude start value 
	if i==0:
		init_b = 9.3				# Set fit position start value for zero peak
	else:							
		init_b = 15. + 10.2 * i		# Guess fit postitions for reflection peaks
	init_c = 1.						# Set fit width start value
	for k in range(3):
		if k==0:									# First fit routine
			fitrange = 3.							# Set Fit range around estimated position start values
			init_vals = [init_a, init_b, init_c]	# Set estimated start values as fit initials
			popt, pcov = curve_fit(sinc, PMTangle[(PMTangle>(init_b-fitrange)) & (PMTangle<(init_b+fitrange))],
										 PMTrate[ (PMTangle>(init_b-fitrange)) & (PMTangle<(init_b+fitrange))],
										 p0 = init_vals)
		else:										# Following fit routines
			fitrange = 2.							# Set new fit range
			init_vals = popt						# Set fit parameters of first fit routine as fit initials
			init_b = popt[1]						# Set position start value based on first fit routine
			popt, pcov = curve_fit(sinc, PMTangle[(PMTangle>(init_b-fitrange)) & (PMTangle<(init_b+fitrange))],
										 PMTrate[ (PMTangle>(init_b-fitrange)) & (PMTangle<(init_b+fitrange))],
										 p0 = init_vals,
										 sigma=PMTrate_err[(PMTangle>(init_b-fitrange)) & (PMTangle<(init_b+fitrange))] )
	peak_para[i][0] = popt[0]						# Write fitted amplitudes
	peak_para[i][1] = popt[1]						# Write fitted positions
	peak_para[i][2] = popt[2]						# Write fitted widths
	peak_para[i][3] = np.sqrt( pcov[0][0] )			# Write fitted amplitude errors
	peak_para[i][4] = np.sqrt( pcov[1][1] )			# Write fitted position errors
	peak_para[i][5] = np.sqrt( pcov[2][2] )			# Write fitted width errors

#####################
###   Integrate   ###
#####################

# Initiate arrays
intrange = 9.											# Integral range around fitted peak position
PeakLowLim	= np.zeros(( len(paths) ), dtype=np.int)	# Lowest peak bin
PeakUpLim	= np.zeros(( len(paths) ), dtype=np.int)	# Highest peak bin
SumPeak		= np.zeros(( len(paths) ), dtype=np.int)	# Overall counts in peak range
BG			= np.zeros(( len(paths) ), dtype=np.int)	# Overall background counts
PeakBG		= np.zeros(( len(paths) ), dtype=np.int)	# Overall background counts in peak
Peak		= np.zeros(( len(paths) ), dtype=np.int)	# Overall counts in peak w/o background

# Fill arrays
for i in range( len(paths) ):
	PeakLowLim[i] = int ( math.ceil(  (peak_para[i][1] - intrange*abs(peak_para[i][2]) - DataLowVal[i]) / 0.45 ) )
	PeakUpLim[i]  = int ( math.floor( (peak_para[i][1] + intrange*abs(peak_para[i][2]) - DataLowVal[i]) / 0.45 ) )
	for j in range( PeakLowLim[i], (PeakUpLim[i]+1) ):
		SumPeak[i] = SumPeak[i] + data[i][j][1]							# Get overall peak counts by integrating over entire peak range
	for j in range( 0, PeakLowLim[i] ):
		BG[i] = BG[i] + abs(data[i][j][1])								# Get overall background by integration over all bins before...
	for j in range( (PeakUpLim[i]+1), DataUpBin[i] ):
		BG[i] = BG[i] + abs(data[i][j][1])								# ...and after peak range
	BG[i] = BG[i]/( DataUpBin[i] - (PeakUpLim[i]+1) + PeakLowLim[i] )	# Get background per bin
	PeakBG[i] = BG[i] * ( PeakUpLim[i]+1 - PeakLowLim[i] )				# Get background within peak range
	Peak[i] = SumPeak[i] - PeakBG[i]									# Get background-cleaned peak counts

#####################
###   Calculate   ###
#####################

# Initiate arrays
Angle		= np.zeros(( len(paths)-1 ))
Angle_err	= np.zeros(( len(paths)-1 ))
Refl		= np.zeros(( len(paths)-1 ))
Refl_err	= np.zeros(( len(paths)-1 ))

# Fill arrays
for i in range( len(paths)-1 ):
	Angle[i] = ( 180 - ( peak_para[i+1][1] - peak_para[0][1] ) )/2					# Formula for incident angle (see also PhD thesis Cecilia)
	Refl[i] = float(Peak[i+1]) / Peak[0]											# Get reflection relative to zero peak reference
	Angle_err[i] = 0.5 * np.sqrt( (peak_para[i+1][4])**2 + (peak_para[0][4])**2 )	# Calculate errors based on Gaussian error propagation
	Refl_err[i] = 1/float(Peak[0]) * np.sqrt( SumPeak[i+1] + PeakBG[i+1] + ( Peak[i+1]/Peak[0] )**2 * (SumPeak[0] + PeakBG[0]) )

######################
###   Plot peaks   ###
######################

# Initiate plot
PeakPlot = pl.figure(num=None, figsize=(14, 7), dpi=80, facecolor='w', edgecolor='k')

# Define RCT colour palette
PlotColor	= [ '#2f4343', '#e30700', '#4387e3', '#47af27', '#ff6f17', '#db3b8f', '#00dddd', '#d7a713' ]
FitColor	= [ '#000000', '#8f0000', '#3f4377', '#1f7b00', '#b74700', '#93074b', '#00a0a0', '#7f571f' ]

# Plot all data and fit functions
for i in range( len(paths) ):
	PMTangle = np.array([ data[i][j][0] for j in range(0,DataUpBin[i]) ])		# Temporary arrays for plot procedure
	PMTrate  = np.array([ data[i][j][1] for j in range(0,DataUpBin[i]) ])
	if i==0:
		pl.plot( PMTangle, PMTrate, linestyle='None', marker='x', color=PlotColor[i%8], label='Zero-peak' )
	else:
		pl.plot( PMTangle, PMTrate, linestyle='None', marker='x', color=PlotColor[i%8], label='Xe-anlge: %s'%XeAngle[i] )
	xfine = np.linspace(peak_para[i][1]-1.5*fitrange, peak_para[i][1]+1.5*fitrange, 10000)					# Temporary array as fit function sample
	pl.plot( xfine, sinc(xfine, peak_para[i][0], peak_para[i][1], peak_para[i][2]), color=FitColor[i%8] )	# Plot fit functions based on sample

#######################
###   Format plot   ###
#######################

# Basic formating
pl.title("SiPM")
pl.xlabel("PMT-angle [deg]")
pl.ylabel("PMT rate [Hz]")
pl.xlim(-10,160)
pl.ylim(0,14000)
pl.xticks(np.arange(0, 170, 10))

# Legend and grid
pl.legend(loc='upper right')
pl.legend(ncol=2, handleheight=2., labelspacing=0.02)
pl.grid()

# Save and show plot
pl.savefig("SiPM_peaks.pdf")
pl.show()

###########################
###   Plot reflection   ###
###########################

# Initiate reflection plot
ReflPlot = pl.errorbar(Angle, Refl, xerr=[Angle_err,Angle_err], yerr=[Refl_err, Refl_err], linestyle='None', marker='x', color='b', label="Reflectivity")

# Basic formating
pl.title("SiPM")
pl.xlabel("Incident angle [deg]")
pl.ylabel("Reflectivity")
pl.xlim(0,90)
pl.ylim(0.1,0.6)
pl.xticks(np.arange(0, 90, 10))

# Legend and grid
pl.legend(loc='upper left')
pl.grid()

# Optional formating
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((60,.1), 30, .6, edgecolor="none", facecolor="blue", alpha=0.2))
currentAxis.text(62.5,0.55,'- Beam shadowed -', fontsize=14)
currentAxis.text(15,0.45,'PRELIMINARY', fontsize=50, color="r", rotation=-45, alpha=0.15)

# Save & Show Plot
pl.savefig("SiPM_refl.pdf")
pl.show()











