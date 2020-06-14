# A simple program to read and plot pulsar data files
# The program uses a manually entered period and phase offset, rather than tempo.
# It does not accommodate the possibilty of of reading a file while it is
# being written.

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import ClassReadDynSpecFileRDMA as crdsf
import argparse
from collections import OrderedDict

# Define argparse arguments
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file", type=str, required=True)
parser.add_argument("-p", "--period", type=float, required=True)
parser.add_argument("-d", "--dm", '--DM', type=float, default=0)
parser.add_argument("-n", "--nbins", '--nBins', type=int, default=200)

args = parser.parse_args()

file=args.file
period=args.period
dm=args.dm
nbins=args.nbins

plt.rcParams.update({'font.size': 18})
#plt.rc('text', usetex=True)
#plt.rcParams['text.latex.preamble'] = [r'\boldmath']

# Change argument's data type
period = float(period)

# Define the number of phase bins and frequencies
DM, nBins, nFreqs = float(dm), int(nbins), 512

def shift_col_up(phase_num,n_rows):
    sumArray[:,phase_num] = np.roll(sumArray[:,phase_num], -n_rows)
def estimate_S_N(averagePhase,mask=np.array([])):
    """Simple estimation of signal_to_noise, with optional masking.
    If mask not given, then all channels will be used in estimating noise
    (will drastically underestimate S:N! Not robust to outliers!)
    Inputs: averagePhase (1D array)
    mask: optional mask of ignored channels (non-zero values are not used by estimator)"""
    if mask.size == 0:
        mask = np.zeros_like(averagePhase)

    noise = np.std((averagePhase[2:]-averagePhase[:-2])[mask[1:-1] == 0])/np.sqrt(2)
    background = np.nanmean(averagePhase[mask == 0])
    return (averagePhase-background)/noise

def sizeof_fmt(num, suffix='B'):
    for unit in ['',' Ki',' Mi',' Gi',' Ti',' Pi',' Ei',' Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

# Create a new instance of the file reader
f = crdsf.RDMASpecFile()

# Specify the path and the file
path = './'

# We take the full frequency range
fmin, fmax = 0., 9999.
frange = (fmin * 1.0E6, fmax * 1.0E6)

# Open the specified file
fi = f.gopenrt(file, path=path, frange=frange)

with open(path+file[:14]+'.inf', 'r') as metadata:
    for line in metadata:
        if 'dt' in line:
            dt = float(list(OrderedDict.fromkeys(line.split(' ')))[2])
        elif 'start_utc' in line:
            start_date = str(list(OrderedDict.fromkeys(line.split(' ')))[2])
            start_time = str(list(OrderedDict.fromkeys(line.split(' ')))[3])
        elif 'duration' in line:
            nRows = int(list(OrderedDict.fromkeys(line.split(' ')))[2])
        elif 'ChanSel' in line:
            polChan = str(list(OrderedDict.fromkeys(line.split(' ')))[2])
            if polChan == '1':
                polChan = 'A'
            elif polChan == '2':
                polChan = 'B'
            elif polChan == '3':
                polChan = 'AB'
            else:
                polChan = 'Unknown ("'+polChan+'")'

obs_duration = str(nRows*dt)

# prepare arrays of phase vs. frequency
sumArray    = np.zeros((nBins,nFreqs))
countsArray = np.zeros(nBins)

#Fill the Phase vs. Frequency array
i = 0

# The value for the period below is hardcoded for J0332+5434
# for the particular file that has been chosen.  The precise value
# will vary depending on the time etc.

while True:
    data = f.gnextrt()		# get new data
    if data[0] == 0: break 	# if no new data, then exit
    ts = data[1]                # the array of times.  NOTE: times are in Modified Julian Day (MJD)
    if i == 0 : ts0 = ts[0]     # save the first time reading
    signal = data[2]            # the FFT

    times = 86400.*(ts - ts0)
    phases = np.fmod(times,period)/period
    iPhases = np.int64(nBins*phases)
    for j, iPhase in enumerate(iPhases):
        sumArray[iPhase] += np.flipud(signal[j][1:])
        countsArray[iPhase] += 1
    i += 1

# Compute average FFT for each phase
for i in range(nBins):
    sumArray[i] /= countsArray[i]

# Compute the average signal for each frequency bin (zero DM)
averageFFT = np.average(sumArray,axis=0)

# Fold data (zero DM)
no_dm_averagePhase = np.nanmean(sumArray,axis=1)

# Plot averageFFT
freqs = np.linspace(1275.,1525.,len(averageFFT))

# Mask channels wiht RFI
#for i in range(sumArray.shape[1]):
#    if np.mean(sumArray[:, i]) < 0.3:
#        sumArray[:, i] = np.nan


# Perform de-dispersion
deltaF = float(np.max(freqs)-np.min(freqs))/sumArray.shape[1]
f_start = np.min(freqs)
for phase_num in range(sumArray.shape[1]):
    f_chan = f_start+(phase_num)*deltaF
    deltaT = 4149*dm*((1/(f_chan**2))-(1/(np.max(freqs)**2)))
    n = int((float(deltaT)/(float(1)/sumArray.shape[0])))
    shift_col_up(phase_num,n)

# Compute the average signal for each frequency bin
averageFFT = np.average(sumArray,axis=0)

# Subtract the averaged spectrum from all rows
for i in range(nBins):
    sumArray[i] -= averageFFT

# Fold data (de-dispersed)
averagePhase = np.nanmean(sumArray,axis=1)

# Get the index of the pulse (i.e. global maximum)
max_index = np.argmax(averagePhase)

# Same for zero DM
no_dm_max_index = np.argmax(no_dm_averagePhase)

fig = plt.figure(figsize=(20,15))
gs = gridspec.GridSpec(nrows=2, ncols=2, width_ratios=[1, 1.8], height_ratios=[1.5, 3])
#gs.update(wspace=5, hspace=1)

#ax2 = fig.add_subplot(2,1,2)
#xaxis = np.arange(np.min(freqs),np.max(freqs),0.1)

ax3 = fig.add_subplot(gs[1, 1])
ax3.plot(freqs, averageFFT, '#3182bd', linewidth=2.5)
ax3.fill_between(freqs, averageFFT, color='#deebf7')
#ax3.plot(freqs, averageFFT, '#bd3131', linewidth=2.5)
#ax3.fill_between(freqs, averageFFT, color='#f7dede')
ax3.set_xlim(np.min(freqs),np.max(freqs))
ax3.set_title('Averaged Spectrum')
ax3.set_xlabel('Frequency (MHz)')
ax3.set_ylabel('Power (arb. units)')
ax3.grid()

ax4 = fig.add_subplot(gs[0, 1])
ax4.set_title('Observation Parameters')
ax4.annotate(r'Observation name: '+file+' ('+str(sizeof_fmt(os.path.getsize(path+file)))+')\nTelescope: TLM-18\nDatetime of first sample (UTC): '+start_date+' '+start_time+'\nObserving duration: '+obs_duration+' sec\n$\mathrm{T}_{sample}$: '+str(dt)+' sec\nFFT samples folded: '+str(nRows)+'\nObserving band: L ('+str(np.min(freqs))+'-'+str(np.max(freqs))+' MHz)\nObserving bandwidth: '+str(np.max(freqs)-np.min(freqs))+' MHz\nFrequency channels: '+str(nFreqs)+'\nPolarization channel: '+polChan+'\n\nFolding period $(P_{0})$: '+str(period)+' sec $(f_{0} = $'+str(1/period)+' Hz$)$\nApplied DM: '+str(DM)+' $pc/cm^{3}$\nDetected pulse phase: '+str(float(max_index)/nBins), xy=(0, 312), xycoords='axes points', size=18, ha='left', va='top', color='black')
ax4.axis('off')

# Apply mask for S/N estimator
mask = np.zeros_like(averagePhase)
for i in range(max_index-int(nBins/10),max_index+int(nBins/10)):
    mask[i] = 1

no_dm_mask = np.zeros_like(no_dm_averagePhase)
for i in range(no_dm_max_index-int(nBins/10),no_dm_max_index+int(nBins/10)):
    no_dm_mask[i] = 1

# Get the value of the S/N of the pulse and append it to file
max_snr = np.amax(estimate_S_N(averagePhase, mask))

with open('results.txt', 'a') as myfile:
    myfile.write(str(nBins)+','+str(dm)+','+str(period)+','+str(max_snr)+'\n')

no_dm_max_snr = np.amax(estimate_S_N(no_dm_averagePhase, no_dm_mask))

# Rotate sumArray by 90 deg to match Phase axis values
sumArray = np.rot90(sumArray)

# Flip sumArray to match Frequency axis values
sumArray = np.flipud(sumArray)

# Define frequency step for axis
freq_step = float(np.max(freqs)-np.min(freqs))/3

#ax3 = fig.add_subplot(2,2,4)
ax2 = fig.add_subplot(gs[1, 0])
ax2.imshow(sumArray, origin="lower", interpolation="None", aspect="auto", extent=[0, 1, np.min(freqs), np.max(freqs)])
ax2.set_xlabel('Phase')
ax2.set_ylabel('Frequency (MHz)')

#cbar = fig.colorbar(cax)
#cbar.ax.tick_params(size=0)
#cbar.set_ticks([])

# Plot Pulse Profile
ax1 = fig.add_subplot(gs[0, 0])
ax1.axvline(x=max_index, color='brown', linestyle='--', linewidth=2)
ax1.plot(estimate_S_N(averagePhase, mask),'ro-')
ax1.plot(estimate_S_N(no_dm_averagePhase, mask),'go-')
ax1.annotate('S/N (de-dispersed) = '+str('{:.2f}'.format(max_snr)), xy=(8, 316), xycoords='axes points', size=18, ha='left', va='top', color='brown')
ax1.annotate('S/N (zero DM) = '+str('{:.2f}'.format(no_dm_max_snr)), xy=(8, 295), xycoords='axes points', size=18, ha='left', va='top', color='darkgreen')
ax1.set_title('Pulse Profile')
ax1.set_ylabel('Signal-to-Noise Ratio (S/N)')
ax1.grid()

ax1.get_shared_x_axes().join(ax1, ax3)
ax1.set_xticklabels([])

# Plot dispersion limits
left_index = np.where(estimate_S_N(averagePhase, mask) >= 5)[0][0]
right_index = np.where(estimate_S_N(averagePhase, mask) >= 5)[0][-1]

ax2.axvline(x=float(left_index)/nBins, color='lime', linewidth=1.5)
ax2.axvline(x=float(right_index)/nBins, color='lime', linewidth=1.5)

# Save to file
plt.tight_layout()
plt.subplots_adjust(wspace=0.15, hspace=0.01)
plt.savefig("plot.png")
