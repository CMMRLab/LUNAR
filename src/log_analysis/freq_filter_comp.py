# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
April 3rd, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This script simply tests different frequency based filtering methods for stress-strain
data. If using from the Spyder IDE these commands can be useful for figure pop-ups:

    %matplotlib qt
    %matplotlib inline


Resources for generating stress-strain data:
    https://en.wikipedia.org/wiki/Ramberg%E2%80%93Osgood_relationship
    https://www.quadco.engineering/en/know-how/fea-ramberg-osgood-equation.htm
    https://learnfea.com/stress-strain-curve-approximation-2/
"""

##############################
# Import Necessary Libraries #
##############################
import src.log_analysis.signal_processing as signals
import matplotlib.pyplot as plt
import numpy as np


# Generate stress and strain 
max_strain = 0.5
max_stress = 150
strain = np.linspace(0, 0.5, 2500)
stress  = max_stress*(1 - np.exp(-strain*max_strain/(max_strain**4)))


# Add noise to signal
np.random.seed(12345)
dx = np.max(strain) - np.min(strain)
dy = np.max(stress) - np.min(stress)
x_mag = dx/1000
y_mag = dy/6
y_noisy = stress + y_mag*np.random.randn(strain.size) + y_mag*np.sin(2*np.pi*x_mag)


# Setup figure options for filtering
savefig = False
figname = ''
dpi = 300


# Filter with the inverse FFT method
iFFT_qm = 'msr'
iFFT_threshold = 'mean:1.0'
y_iFFT, qm_iFFT = signals.iFFT_filter(strain, y_noisy, iFFT_threshold, iFFT_qm, savefig, figname, dpi, plot_PSD=False)


# Filter with the Butterworth method
BW_wn = 'op'
BW_qm = 'msr'
BW_order = 2
write_data = False
if BW_wn.startswith('op'):
    BW_wn = signals.butter_optimize_wn_with_power_spectrum(strain, y_noisy, BW_wn, write_data, savefig, figname+'_FFT_PSD_wn', dpi)
y_BW, qm_BM = signals.butter_lowpass_filter(strain, y_noisy, BW_wn, BW_order, BW_qm, write_data, savefig, figname, dpi, pflag=True)



# Plot results
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 9))
ax1.plot(strain, y_noisy, '-', lw=2.0, color='tab:cyan', label='Noisy Signal')
ax1.plot(strain, stress, '-', lw=2.0, color='tab:blue', label='True Signal')
ax1.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
ax1.set_xlabel('Time')
ax1.set_ylabel('Signal')

ax2.plot(strain, y_iFFT, '-', lw=2.0, color='tab:red', label='iFFT filter qm={}={}'.format(iFFT_qm, qm_iFFT))
ax2.plot(strain, y_BW, '-', lw=2.0, color='tab:orange', label='Butterworth filter qm={}={}; order={}'.format(BW_qm, BW_order, qm_BM))
ax2.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
ax2.set_xlabel('Time')
ax2.set_ylabel('Cleaned Signal')

ax3.plot(strain, stress, '-', lw=2.0, color='tab:blue', label='True Signal')
ax3.plot(strain, y_iFFT, '-', lw=2.0, color='tab:red', label='iFFT filter qm={}={}'.format(iFFT_qm, qm_iFFT))
ax3.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
ax3.set_xlabel('Time')
ax3.set_ylabel('Cleaned Signal')

ax4.plot(strain, stress, '-', lw=2.0, color='tab:blue', label='True Signal')
ax4.plot(strain, y_BW, '-', lw=2.0, color='tab:orange', label='Butterworth filter qm={}={}; order={}'.format(BW_qm, BW_order, qm_BM))
ax4.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
ax4.set_xlabel('Time')
ax4.set_ylabel('Cleaned Signal')
fig.tight_layout()