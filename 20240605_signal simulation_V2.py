import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.special import wofz
from scipy.signal import butter, filtfilt

# settings
fs = 51200  # sampling frequency 51200 in LabView
f_triangle = 0  # in Hz, not working so set to 0 for now
f_sin = 200 # in Hz, set to 1000 at the moment
ampl_triangle = 0  # in nm, 0.25, note that inertia of laser is not taken into account
ampl_sin = 0.25  # in nm, 0.1 realistic with triangle, 0.25 without
modulation_offset = 2003.50505  # center wavelength
loss = 0  # 0 no loss, 1 no signal
time_constant = 0.001  # in seconds, cutoff_freq = 1 / (2 * np.pi * time_constant), currently 100ms in LabView
phase_shift = 1 / (3 *  np.pi * f_sin) # between reference and detector signal, should be set to maximize filtered signal
filter_order = 5 #  filter_oder * 20 dB/decade, 10 für 200dB, 5 für 100dB
t = np.arange(0, 1, 1/fs)  # step size: 1/sampling frequency

# modulation signal
def mod_triangle(t):
    return ampl_triangle * signal.sawtooth(2 * np.pi * f_triangle * t,width=0.5)
def mod_sin(t):
    return ampl_sin * np.sin(2 * np.pi * f_sin * t)
def modulation(t):
    return mod_triangle(t) + mod_sin(t)

def voigt_profile(x, amplitude, center, sigma, gamma):
    z = ((x - center) + 1j * gamma) / (sigma * np.sqrt(2))
    v = amplitude * wofz(z).real / (sigma * np.sqrt(2 * np.pi))
    return v

# HITRAN data
with open("C:/Users/alexi/Documents/Python/20240605_Signal simulation triangle sinus/1p CO2 transmission.txt") as f:
    lines = f.readlines()
    CO2trans, CO2abs, wavenumber, wavelength = [], [], [], []
    nblines = len(lines)
    print("number lines:", nblines)
    for line in lines[7:(nblines-0)]:   #Bereich anpassen  167-160 für nur mittleren peak
        wavenumber.append(float(line.split()[0]))
        CO2abs.append(1-float(line.split()[1]))
        CO2trans.append(float(line.split()[1]))
        wavelength.append(10000000/float(line.split()[0]))  #in nm

# voigt fit 2003 nm peak
initial_guess = (0.023975662, 2003.00076, 0.0000028, 0.032)  # Initial parameter guess amplitude center sigma gamma
params, covariance = curve_fit(voigt_profile, wavelength, CO2abs, p0=initial_guess)
# Extract the fitted parameters
amplitude_fit, center_fit, sigma_fit, gamma_fit = params
print("Amplitude:", amplitude_fit, "Center:", center_fit, "Sigma:", sigma_fit, "Gamma:", gamma_fit)

# voigt fit 2003.5 nm peak
initial_guess2 = (0.024621842, 2003.50505, 0.0000028, 0.032)  # Initial parameter guess amplitude center sigma gamma
params2, covariance2 = curve_fit(voigt_profile, wavelength, CO2abs, p0=initial_guess2)
# Extract the fitted parameters
amplitude_fit2, center_fit2, sigma_fit2, gamma_fit2 = params2
print("Amplitude2:", amplitude_fit2, "Center2:", center_fit2, "Sigma2:", sigma_fit2, "Gamma2:", gamma_fit2)

# voigt fit 2004 nm peak
initial_guess3 = (0.024902473, 2004.02107, 0.0000028, 0.032)  # Initial parameter guess amplitude center sigma gamma
params3, covariance3 = curve_fit(voigt_profile, wavelength, CO2abs, p0=initial_guess3)
# Extract the fitted parameters
amplitude_fit3, center_fit3, sigma_fit3, gamma_fit3 = params3
print("Amplitude3:", amplitude_fit3, "Center3:", center_fit3, "Sigma3:", sigma_fit3, "Gamma3:", gamma_fit3)

# summ of the 3 voigt fits
def voigt(x):
    y = voigt_profile(x, *params) + voigt_profile(x, *params2) + voigt_profile(x, *params3)
    return y

# find absorption peaks, use for voigt fit
peaks = find_peaks(CO2abs, height=0.2)
peakheight = peaks[1]['peak_heights']  # list containing the height of the peaks
print("y_peaks:", peakheight)
peakpos, peakposlabel = list(), list()
for j in range(len(peaks[0])):
    peakpos.append(wavelength[peaks[0][j]])
    peakposlabel.append(round(peakpos[j], 5))
print("x_peaks:", peakposlabel)

# detector signal, 100- because voigt is over absorption, + modulation_offset because modulation is around 0
def detectorCO2signal(t):
    y = 100 - (1 - loss) * (voigt(modulation(t) + modulation_offset))
    return y

# Remove the offset
def detectorCO2signal_corr(t):
    return detectorCO2signal(t + phase_shift) - np.mean(detectorCO2signal(t))  # is this nessesary?

# Multiply the noisy signal by the reference signal
mixed_signal1f = detectorCO2signal_corr(t) * mod_sin(t)
mixed_signal2f = detectorCO2signal_corr(t) * mod_sin(2 * t)
# mixed_signal2f = mod_sin(t) * mod_sin(t)  # to test filter

# Low-pass filter, calculate the cutoff frequency from the time constant
cutoff_freq = 1 / (2 * np.pi * time_constant)  # für IRR Filter
print('cutoff freq:', cutoff_freq)
b, a = butter(filter_order, cutoff_freq / (0.5 * fs), btype='low')

# Apply the low-pass filter to the mixed signal
filtered_signal1f = filtfilt(b, a, mixed_signal1f) # is 0 because only CO2 effect is simulated
filtered_signal2f = filtfilt(b, a, mixed_signal2f)

# plot CO2 absorption and voigt fit
# """
plt.figure(figsize=(6, 3.5))
plt.plot(wavelength, CO2abs, label=r"1% CO$_2$ rest N$_2$")
#parameters = np.array((0.23975662, 2003, 0.5, 0.2)) # zum manuell testen
#plt.plot(wavelength, voigt_profile(wavelength, *parameters), 'r-', label='Voigt fit')
plt.plot(wavelength, voigt(wavelength), 'k--', label='Voigt fit')
#plt.plot(wavelength, voigt_profile(wavelength, *params2), 'r-', label='Voigt fit')
plt.legend(loc="best")
plt.xticks(np.arange(2003, 2004.5, 0.5))
plt.xlabel("Wavelength in nm")
plt.ylabel("Absorption in % (2 m)")
plt.legend(loc="upper right")
plt.tight_layout()
plt.show()
# """

# plot modulation around 2003.5 nm
# """
fig, ax = plt.subplots(figsize=(6,3))
plt.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
ax.plot(t, modulation(t))
# ax.plot(t, detectorCO2signal_corr(t))  # to visualize phase shift
# plt.yticks([2003,2003.5, 2004], ['2003.0','2003.5','2004.0'], rotation=0)
plt.ylabel('Wavelength in nm +2003.5')
plt.xlabel('Time in s')
fig.tight_layout()
plt.show()
# """

# plot detector signal
# """
fig, ax = plt.subplots(figsize=(5,3))
ax.plot(t, detectorCO2signal(t))
plt.title('CO2 detector signal')
plt.ylabel(' ')
plt.xlabel('Time in s')
fig.tight_layout()
plt.show()
# """

# plot 1f 2f signal
# """
plt.figure(figsize=(8, 4))

plt.subplot(2, 1, 1)
plt.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
plt.plot(t, mixed_signal2f)
plt.ylabel('Mixed 2f signal')
plt.xlabel('Time in s')

plt.subplot(2, 1, 2)
plt.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
plt.plot(t, filtered_signal2f)
plt.ylabel('Filtered 2f signal')
plt.xlabel('Time in s')
plt.tight_layout()
plt.show()
# """