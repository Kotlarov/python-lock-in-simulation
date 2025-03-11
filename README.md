# python-lock-in-simulation
This project is a simulation of the measurement of CO2 with wavelength-modulated spectroscopy. I worked on this experiment at the Technical University of Munich and wanted a simulation to answer my questions about the results I got. The experiment consisted of a tube filled with a set concentration of CO2 and N2 with an Infrared  (IR) Laser to measure the concentration. Lock-in amplification was used to increase the sensitivity: the IR laser was modulated from 2003.25 to 2003.75 nm right over an absorption peak, see Figures 1 and 2. This induces a signal at the detector with double the frequency (Figure 3). The Lock-in amplifier acts as a strong frequency filter, thus increasing the signal-to-noise ratio. The Script can be separated into 5 steps:

1. Use HITRAN to get CO2 absorption data and do a Voigt fit around the wavelength of interest (Figure 1). 
2. Modulate the laser wavelength with a sin wave. This signal also acts as a reference signal for the lock-in amplifier (Figure 2).
3. Simulate a signal that is proportional to the actual detector signal. This is done by looking at the composition of the Voigt fit and the modulation: Voigt(modulation(time)) (Figure 3).
4. Simulate the Lock-In-amplifier:
   - Correct the offset of the detector signal and multiply it by the reference signal. This is called the mixed 2f signal "2f" for double frequency.
   - Apply a low pass filter to the mixed 2f signal. This can be seen as a moving average over the mixed signal. The filtered 2f signal takes a short time to settle and results in a constant signal that is  proportional to the CO2 concentration (Figure 4).
