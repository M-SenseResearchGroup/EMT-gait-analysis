# Walk and Turn.
# Code initially by Chris Petrillo
# Python Version: Lukas Adamowicz
# 7/25/18


import numpy as np
from scipy import signal
import pywt  # pywavelets
import sys
sys.path.append("C:\\Users\\Lukas Adamowicz\\Dropbox\\Masters\\MC10py")
import MC10py
import matplotlib.pyplot as pl
from matplotlib.patches import Patch

# def load_mc10(study_dir, segment=True, sync=True, save=True, save_loc=None, save_subj=False, return_data=True):
# paths = MC10py.LoadMC10('C:\\Users\\Lukas Adamowicz\\Documents\\Study Data\\EMT\\ASYM_OFFICIAL', segment=True,
#                         sync=True, save=True, save_loc='import', save_subj=False, return_data=False)

raw_data = MC10py.OpenMC10('C:\\Users\\Lukas Adamowicz\\Documents\\Study Data\\EMT\\ASYM_OFFICIAL\\data.pickle')

# plot walk and turn data for each subject

# Automatic split between before and after turn (ie turn detection)


def turn_detect(data, plot=False):
    pl.close('all')
    wdb, wda = signal.butter(1, 1.5/62.5, 'lowpass')  # filter gyroscopic data to see turn clearly

    ret_data = dict()
    for sub in data.keys():
        if plot:
            f, ax = pl.subplots(4, figsize=(14, 11))
        i = 0
        try:
            events = data[sub]['sacrum']['gyro'].keys()
            ret_data[sub] = dict()
            for ev in events:
                if 'Walk and Turn' in ev:
                    ret_data[sub][ev] = {'gyro': {'1': dict(), '2': dict()}, 'accel': {'1': dict(), '2': dict()}}
                    iz = np.argmax(abs(np.mean(data[sub]['sacrum']['accel'][ev][:100, 1:], axis=0))) + 1  # z axis index
                    time = data[sub]['sacrum']['gyro'][ev][:, 0]

                    fd = abs(signal.filtfilt(wdb, wda, data[sub]['sacrum']['gyro'][ev][:, iz]))
                    mfd = np.mean(fd)
                    _ips = signal.argrelmax(fd, order=63, mode='wrap')
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
                    if plot:
                        ax[i].plot(time, fd, label='Filtered')
                        ax[i].plot(time[ip], fd[ip], 'ro', label='max')
                        ax[i].plot(time[it], fd[it], 'ko', label='min')
                        ax[i].axhline(1.25*mfd, color='k', ls='--', label='1.25 * mean(filt)')

                    if len(it) == 2:
                        for loc in data[sub].keys():

                            ret_data[sub][ev]['accel']['1'][loc] = data[sub][loc]['accel'][ev][:it[0], :]
                            ret_data[sub][ev]['gyro']['1'][loc] = data[sub][loc]['gyro'][ev][:it[0], :]

                            ret_data[sub][ev]['accel']['2'][loc] = data[sub][loc]['accel'][ev][it[1]:, :]
                            ret_data[sub][ev]['gyro']['2'][loc] = data[sub][loc]['gyro'][ev][it[1]:, :]
                    else:
                        for loc in data[sub].keys():
                            ret_data[sub][ev]['accel']['1'][loc] = data[sub][loc]['accel'][ev][:it[0], :]
                            ret_data[sub][ev]['gyro']['1'][loc] = data[sub][loc]['gyro'][ev][:it[0], :]

                            ret_data[sub][ev]['accel']['2'][loc] = data[sub][loc]['accel'][ev][it[1]:it[2], :]
                            ret_data[sub][ev]['gyro']['2'][loc] = data[sub][loc]['gyro'][ev][it[1]:it[2], :]

                    if plot:
                        ax[i].plot(ret_data[sub][ev]['gyro']['1'][:, 0], ret_data[sub][ev]['gyro']['1'][:, iz],
                                   alpha=0.5, linewidth=5, color='r')
                        ax[i].plot(ret_data[sub][ev]['gyro']['2'][:, 0], ret_data[sub][ev]['gyro']['2'][:, iz],
                                   alpha=0.5, linewidth=5, color='r')
                    i += 1
        except KeyError:
            pass

    return ret_data


def extract_gait_params(data, plot=False):
    gait = {i: {j: dict() for j in data[i].keys()} for i in data.keys()}

    pl.close('all')
    for sub in data.keys():
        if plot:
            f, ax = pl.subplots(nrows=4, ncols=2, figsize=(16, 9))
            f.suptitle(f'Subject {sub}')
        n = 0  # figure tracking (for events)
        for ev in data[sub].keys():
            m = 0  # figure tracking (for left/right)
            gait[sub][ev] = dict()
            sens = 'proximal_lateral_shank'
            for side in ['left', 'right']:
                gait[sub][ev][side] = dict()
                loc = sens + '_' + side

                time = data[sub][ev]['gyro']['1'][loc][:, 0]

                wave = 'mexh'
                coefs, freq = pywt.cwt(data[sub][ev]['gyro']['1'][loc][:, 3], np.arange(1, 65), wave,
                                       sampling_period=1 / 62.5)

                m_ind = np.argmax(abs(coefs[5, :]))
                mc = np.mean(abs(coefs[5, :]))
                if coefs[5, m_ind] < 0:
                    pks, _ = signal.find_peaks(-coefs[5, :], height=mc, distance=25)
                else:
                    pks, _ = signal.find_peaks(coefs[5, :], height=mc, distance=25)

                # rough cadence frequency
                cf = 1 / (np.mean(np.diff(time[pks])) / 1000)
                ic = np.argmin(abs(freq - (2 * cf)))

                mc = np.mean(abs(coefs[ic, :]))
                if coefs[5, m_ind] < 0:
                    lmx, _ = signal.find_peaks(-coefs[ic, :], height=mc, distance=25)
                    lmn, _ = signal.find_peaks(coefs[ic, :], distance=25)
                else:
                    lmx, _ = signal.find_peaks(coefs[ic, :], height=mc, distance=25)
                    lmn, _ = signal.find_peaks(-coefs[ic, :], distance=25)

                gait[sub][ev][side]['step time'] = np.diff(time[lmx])/1000  # still have to convert to seconds
                _tr = []
                for pk in lmx:
                    clind = np.argmin(abs(lmn - pk))
                    if lmn[clind] < pk:
                        _tr.append(lmn[clind])
                        try:
                            _tr.append(lmn[clind + 1])
                        except:
                            pass
                    else:
                        _tr.append(lmn[clind - 1])
                        _tr.append(lmn[clind])
                tr = np.unique(_tr)

                _ext = np.append(tr, lmx)
                _extt = np.append(['min'] * len(tr), ['max'] * len(lmx))

                srt = np.argsort(_ext)
                ext = _ext[srt]
                extt = _extt[srt]

                stance = []
                swing = []
                gait[sub][ev][side]['stance time'] = []
                gait[sub][ev][side]['swing time'] = []
                for i in range(len(ext) - 2):
                    if extt[i] == 'min' and extt[i + 1] == 'max':
                        swing.append((ext[i], ext[i + 2]))
                        gait[sub][ev][side]['swing time'].append(time[ext[i+2]] - time[ext[i]])
                    elif extt[i] == 'min' and extt[i + 1] == 'min':
                        stance.append((ext[i], ext[i + 1]))
                        gait[sub][ev][side]['stance time'].append(time[ext[i+1]] - time[ext[i]])

                gait[sub][ev][side]['stance time'] = np.asarray(gait[sub][ev][side]['stance time'])
                gait[sub][ev][side]['swing time'] = np.asarray(gait[sub][ev][side]['swing time'])

                if plot:
                    for st in stance:
                        ax[n, m].plot(time[st[0]:st[1]], coefs[ic, st[0]:st[1]], linewidth=8, alpha=0.4, color='r')
                    for sw in swing:
                        ax[n, m].plot(time[sw[0]:sw[1]], coefs[ic, sw[0]:sw[1]], linewidth=8, alpha=0.4, color='g')

                    ax[n, m].plot(time, data[sub][ev]['gyro']['1'][loc][:, 3], alpha=0.6)
                    line1, = ax[n, m].plot(time, coefs[ic, :], label=f'{freq[ic]:.1f} Hz')
                    ax[n, m].set_title(f'{ev}\n{loc}')

                    red = Patch(color='r', alpha=0.4, label='Stance')
                    green = Patch(color='g', alpha=0.4, label='Swing')
                    ax[n, m].legend(handles=[line1, red, green])
                m += 1
            n += 1

        if plot:
            f.tight_layout()

    return gait


spld = turn_detect(raw_data, plot=False)  # split data about turn


gp = extract_gait_params(spld, plot=False)
