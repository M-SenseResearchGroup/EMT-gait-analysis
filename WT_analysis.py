# Walk and Turn.
# Code initially by Chris Petrillo
# Python Version: Lukas Adamowicz
# 7/25/18

import sys
sys.path.append("C:\\Users\\Lukas Adamowicz\\Dropbox\\Masters\\MC10py")
import MC10py
# from scipy.signal import butter, freqz
# from numpy import loadtxt
# import matplotlib.pyplot as pl

# def load_mc10(study_dir, segment=True, sync=True, save=True, save_loc=None, save_subj=False, return_data=True):
paths = MC10py.LoadMC10('C:\\Users\\Lukas Adamowicz\\Documents\\Study Data\\EMT\\ASYM_OFFICIAL', segment=True,
                        sync=True, save=True, save_loc='import', save_subj=False, return_data=False)

# bl, al = butter(1, 8/62.5, 'lowpass')  # half order since will use filtfilt
