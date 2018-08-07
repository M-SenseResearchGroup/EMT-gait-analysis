# Walk and Turn.
# Code initially by Chris Petrillo
# Python Version: Lukas Adamowicz
# 7/25/18


import numpy as np
from scipy.signal import butter, filtfilt, argrelmax, argrelmin
import sys
sys.path.append("C:\\Users\\Lukas Adamowicz\\Dropbox\\Masters\\MC10py")
import MC10py
import matplotlib.pyplot as pl

# def load_mc10(study_dir, segment=True, sync=True, save=True, save_loc=None, save_subj=False, return_data=True):
# paths = MC10py.LoadMC10('C:\\Users\\Lukas Adamowicz\\Documents\\Study Data\\EMT\\ASYM_OFFICIAL', segment=True,
#                         sync=True, save=True, save_loc='import', save_subj=False, return_data=False)

data = MC10py.OpenMC10('C:\\Users\\Lukas Adamowicz\\Documents\\Study Data\\EMT\\ASYM_OFFICIAL\\data.pickle')

# plot walk and turn data for each subject

# Automatic split between before and after turn (ie turn detection)

def turn_detect(data):
    pl.close('all')
    wdb, wda = butter(1, 1.5/62.5, 'lowpass')  # filter gyroscopic data to see turn clearly

    ret_data = dict()
    for sub in data.keys():
        f, ax = pl.subplots(4, figsize=(14, 11))
        i = 0
        try:
            events = data[sub]['sacrum']['gyro'].keys()
            ret_data[sub] = dict()
            for ev in events:
                ret_data[sub][ev] = dict()
                if 'Walk and Turn' in ev:
                    iz = np.argmax(abs(np.mean(data[sub]['sacrum']['accel'][ev][:100, 1:], axis=0))) + 1  # z axis index
                    time = data[sub]['sacrum']['gyro'][ev][:, 0]
                    # ax[i].plot(time, data[sub]['sacrum']['gyro'][ev][:, iz])

                    fd = abs(filtfilt(wdb, wda, data[sub]['sacrum']['gyro'][ev][:, iz]))
                    mfd = np.mean(fd)
                    _ips = argrelmax(fd, order=63, mode='wrap')
                    ip = _ips[0][np.argwhere(fd[_ips] > 2*mfd)].flatten()
                    it = []
                    for j, k in zip(ip, range(len(ip))):
                        try:
                            it.append(np.argwhere(fd[:j] < 1.25 * mfd)[-1])
                        except:
                            ip = np.delete(ip, k)
                            continue
                        try:
                            it.append(np.argwhere(fd[j:] < 1.25 * mfd)[0] + j)
                        except:
                            pass
                    it = np.asarray(it).flatten()
                    ax[i].plot(time, fd)
                    ax[i].plot(time[ip], fd[ip], 'ro')
                    ax[i].plot(time[it], fd[it], 'ko')
                    ax[i].axhline(1.25*mfd, color='k', ls='--')
                    i += 1
        except KeyError:
            pass


turn_detect(data)
# bl, al = butter(1, 8/62.5, 'lowpass')  # half order since will use filtfilt
