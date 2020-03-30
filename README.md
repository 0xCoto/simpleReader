### NOTE: This work is currently under development and requires further documentation.

# simpleReader
**simpleReader** is the software used for the analysis of pulsar data obtained with the ISEC TLM-18 Telescope. Unlike other tools like PRESTO (which definitely offers a richer variety of plots), it calculates and displays the average spectrum of the observation, which makes RFI detection easier. Most of the data processing is carried out using [`NumPy`](https://numpy.org/) arrays, which makes the algorithms run faster, with a relatively low computational expense. The two main functionalities are **epoch folding**, **S/N detection**, **visual indication of the pulse width/DM**, **Power vs Frequency vs Phase**, and **automatic pulse phase detection**, which will be useful for glitch detection in the near future (although a method for pulse phase calibration shall be investigated in order to account for imperfections in the timing accuracy of the instrumentation of the telescope (particularly the superheterodyne receiver)). **Power vs Frequency vs Time** (subtracts data unassociated with the pulsar pulses) is already built, and will soon be implemented. PoC example:
<p align="center">
  <img src="https://i.imgur.com/xE3SMTd.png" alt="Power vs Frequency vs Time"/>
</p>

Upcoming features (not yet implemented) include:
- **Period determination** (e.g. by taking the Fourier transformation of the time series of the observation in the post-dedispersion state and/or by brute-forcing values around the topocentric period (obtained using e.g. TEMPO) or the barycentric period (taken by e.g. the ATNF Pulsar Database). This can also be done by measuring the distance between individual pulses in the time series, although this is only applicable to observations where the invidivdual-pulse-to-noise ratio is high (i.e. pulsars with a high flux density at 1400 MHz like PSR B0329+54 (J0332+5434)). The code for this has already been written for the Green Bank 20m Telescope and tested successfully:
<p align="center">
  <img src="https://i.imgur.com/e2FTa5Q.png" alt="Proof of concept"/>
</p>

- **Dispersion measure determination** (e.g. by brute-forcing DM values after the appropriate topocentric period has been applied (for appropriate pulse folding), and/or with another less computationally-expensive algorithm (TBD))
- **Incoherent de-dispersion** (by appropriately delaying each frequency channel on the dynamic spectrum (waterfall):
<p align="center">
  <img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/361615c3d3a7874fd98554fb9dbe22cb8267ea36" alt="DM Equation"/>
</p>

(code will borrow [lines 7-18](https://github.com/0xCoto/PSR-Toolkit/blob/master/psr_toolkit.py#L7-L18) in [`psr_toolkit.py`](https://github.com/0xCoto/PSR-Toolkit/blob/master/psr_toolkit.py) from [`PSR-Toolkit`](https://github.com/0xCoto/PSR-Toolkit))

- (More functionalities to be added soon)

### Example observation of PSR B0329+54 (J0332+5434) using simpleReader
<p align="center">
  <img src="https://i.imgur.com/GZXN2tW.png" alt="Example Observation"/>
</p>

### Credits
This tool is built by [Apostolos Spanakis-Misirlis](mailto:0xcoto@protonmail.com) and [Prof. Daniel R. Marlow](mailto:marlow@princeton.edu).
