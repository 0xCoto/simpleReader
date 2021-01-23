### NOTE: This work is currently under development and requires further documentation.

# simpleReader
**simpleReader** is the software used for the analysis of pulsar data obtained with the ISEC TLM-18 Telescope. Unlike other tools like PRESTO (which definitely offers a richer variety of plots), it calculates and displays the average spectrum of the observation, which makes RFI detection easier. Most of the data processing is carried out using [`NumPy`](https://numpy.org/) arrays, which makes the algorithms run faster, with a relatively low computational expense. The two main functionalities are **incoherent de-dispersion**, **epoch folding**, **S/N detection**, **visual indication of the pulse width/DM**, **Power vs Frequency vs Phase**, and **automatic pulse phase detection**, which will be useful for glitch detection in the near future (although a method for pulse phase calibration shall be investigated in order to account for imperfections in the timing accuracy of the instrumentation of the telescope (particularly the superheterodyne receiver)). **Power vs Frequency vs Time** (subtracts data unassociated with the pulsar pulses) is already built, and will soon be implemented. PoC example:
<p align="center">
  <img src="https://i.imgur.com/xE3SMTd.png" alt="Power vs Frequency vs Time"/>
</p>

Upcoming features (not yet implemented) include:
- **Period determination** (e.g. by taking the Fourier transformation of the time series of the observation in the post-dedispersion state and/or by brute-forcing values around the topocentric period (obtained using e.g. TEMPO) or the barycentric period (taken by e.g. the ATNF Pulsar Database). This can also be done by measuring the distance between individual pulses in the time series, although this is only applicable to observations where the invidivdual-pulse-to-noise ratio is high (i.e. pulsars with a high flux density at 1400 MHz like PSR B0329+54 (J0332+5434)). The code for this has already been written for the Green Bank 20m Telescope and tested successfully:
<p align="center">
<img src="https://i.imgur.com/e2FTa5Q.png" alt="Proof of concept"/>
</p>

- - **P-search update:** An initial period determination algorithm has been composed for the pulsar monitoring project conducted by TLM-18. It works by brute-forcing a list of periods (e.g. barycentric period ± 0.1 sec with a user-defined step) and seeing where the S/N maximizes:
<p align="center">
<img src="https://i.imgur.com/OLRfRDL.png" alt="Proof of concept"/>
</p>

The period-search algorithm and the plotting script can be found in `psearch.py` and `period_plot.py` respectively.

- **Dispersion measure determination** (e.g. by brute-forcing DM values after the appropriate topocentric period has been applied (for appropriate pulse folding), and/or with another less computationally-expensive algorithm (TBD))
- **Incoherent de-dispersion** (by appropriately delaying each frequency channel on the dynamic spectrum (waterfall)):
<p align="center">
  <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/361615c3d3a7874fd98554fb9dbe22cb8267ea36" alt="DM Equation"/>
</p>

(code will borrow [lines 7-18](https://github.com/0xCoto/PSR-Toolkit/blob/master/psr_toolkit.py#L7-L18) in [`psr_toolkit.py`](https://github.com/0xCoto/PSR-Toolkit/blob/master/psr_toolkit.py) from [`PSR-Toolkit`](https://github.com/0xCoto/PSR-Toolkit))

- - **Incoherent de-dispersion update:** An incoherent-dedispersion algorithm has been applied to `simpleReader.py`. Both the un-de-dispersed and de-dispersed pulse profiles are shown plotted. The data are de-dispersed with a user-defined DM.

- - **S/N vs Pulse profile bins**
<p align="center">
<img src="https://i.imgur.com/cxqP5jB.png" alt="Proof of concept"/>
</p>

- (More functionalities to be added soon)

### Example observation of PSR B0329+54 (J0332+5434) using simpleReader
<p align="center">
  <img src="https://i.imgur.com/3zknvDu.png" alt="Example Observation"/>
</p>

### Usage
```
$ python simpleReader.py -h
usage: simpleReader.py [-h] -f FILE -p PERIOD [-d DM] [-n NBINS]

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE
  -p PERIOD, --period PERIOD
  -d DM, --dm DM, --DM DM
  -n NBINS, --nbins NBINS, --nBins NBINS
```

To analyze an observation, run:
```
python simpleReader.py -f OBS_FILENAME --period FOLDING_PERIOD --nBins NUMBER_OF_BINS
```
Example:
```
python simpleReader.py -f J0332+5434_194_chB.sdf --period 0.71459 --nBins 200
```

### Credits
This tool is developed by [Apostolos Spanakis-Misirlis](mailto:0xcoto@protonmail.com) and [Prof. Daniel R. Marlow](mailto:marlow@princeton.edu). We would like to acknowledge the help provided by [Norman C. Jarosik](mailto:jarosik@Princeton.EDU) regarding the software used for the acquisition of the data.
