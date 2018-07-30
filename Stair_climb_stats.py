import numpy as np
import scipy.stats as stats
import pickle
import matplotlib.pyplot as pl

def import_stuff():
    path = 'C:\\Users\\Lukas Adamowicz\\Documents\\Study Data\\EMT\\Stair Climb Results\\Results.csv'

    rdata = np.loadtxt(path, delimiter=',', dtype=float, skiprows=3)
    # subject, foot side, bag type, direction, step length, lateral deviation, step height, max swing velocity,
    # foot attack angle, contact time, step time, cadence

    subjs = (np.unique(rdata[:, 0]).astype(int)).astype(str)
    sides = ['right', 'left']
    bagtypes = ['hb', 'lb', 'hs', 'ls']
    wdirs = ['up', 'down']

    data = {str(i): {j: {k: {l: [] for l in wdirs}for k in bagtypes} for j in sides} for i in subjs}

    for i in range(len(rdata[:, 0])):
        subj = str(int(rdata[i, 0]))
        side = sides[int(rdata[i, 1])-1]
        btype = bagtypes[int(rdata[i, 2])-1]
        wdir = wdirs[int(rdata[i, 3])-1]
        data[subj][side][btype][wdir].append(rdata[i, 4:])

    for subj in subjs:
        for side in sides:
            for btype in bagtypes:
                for wdir in wdirs:
                    data[subj][side][btype][wdir] = np.array(data[subj][side][btype][wdir], copy=False)

    return data


# data = import_stuff()
#
# fid = open('..\\Data\\stair_climb_results.pickle', 'wb')
# pickle.dump(data, fid)
# fid.close()

fid = open('..\\Data\\stair_climb_results.pickle', 'rb')
data = pickle.load(fid)
fid.close()

subjs = data.keys()
sides = ['right', 'left']
bagtypes = ['hb', 'lb', 'hs', 'ls']
wdirs = ['up', 'down']

# testing for normality and difference between feet, but two few samples for normality test
# dfs = dict()  # different foot stats
#
# for subj in subjs:
#     dfs[subj] = dict()
#     for bt in bagtypes:
#         dfs[subj][bt] = dict()
#         for wd in wdirs:
#             dfs[subj][bt][wd] = np.zeros((8,))
#             pr = np.zeros(8)
#             pl = np.zeros(8)
#             for i in range(8):
#                 try:
#                     _, pr[i] = stats.normaltest(data[subj]['right'][bt][wd][:, i])
#                     _, pl[i] = stats.normaltest(data[subj]['left'][bt][wd][:, i])
#                     _, dfs[subj][bt][wd][i] = stats.ranksums(data[subj]['right'][bt][wd][:, i],
#                                                              data[subj]['left'][bt][wd][:, i])
#                 except IndexError:
#                     dfs[subj][bt][wd][i] = 1.5
#
#             print(f'{subj}  {bt}  {wdir}\n pr \n pl\n {dfs[subj][bt][wd]}\n')

comb = dict()  # combined left and right feet
for subj in subjs:
    comb[subj] = dict()
    for bt in bagtypes:
        comb[subj][bt] = dict()
        for wd in wdirs:
            comb[subj][bt][wd] = np.concatenate((data[subj]['right'][bt][wd], data[subj]['left'][bt][wd]), axis=0)

# test for normality and difference up vs down
# not enough samples for many, and many not normally distributed
# for subj in subjs:
#     for bt in bagtypes:
#         try:
#             _, pu = stats.normaltest(comb[subj][bt]['up'], axis=0)
#             _, pd = stats.normaltest(comb[subj][bt]['down'], axis=0)
#             print(f'Subject: {subj}  Bagtype: {bt}\n {pu<0.05} \n {pd<0.05}\n')
#         except ValueError:
#             print(f'Subject: {subj}  Bagtype: {bt}\n too few samples \n')

ud_stats = dict()  # up down statistics (do the ascending and descending values come from the same distribution)?
# majority say different distribution.
"""
should use stats.wilcoxon ideally, but there are differences in size between up and down steps, as well as too few
data points for the normalcy tests
"""
for subj in subjs:
    ud_stats[subj] = dict()
    for bt in bagtypes:
        ud_stats[subj][bt] = np.zeros(8)
        for i in range(8):
            try:
                _, ud_stats[subj][bt][i] = stats.ranksums(comb[subj][bt]['up'][:, i], comb[subj][bt]['down'][:, i])
            except IndexError:
                ud_stats[subj][bt][i] = 1.5
        # print(f'{subj}, {bt}, same distribution: {ud_stats[subj][bt]>0.05}\n')

# intra-subject bag type differences, ascending
"""
From skimming, looks about even split of differences for different metrics for all the subjects
"""
bda_stats = dict()
compars = [('hblb', 'hb', 'lb'), ('hbhs', 'hb', 'hs'), ('hbls', 'hb', 'ls'), ('lbhs', 'lb', 'hs'), ('lbls', 'lb', 'ls'),
           ('hsls', 'hs', 'ls')]
for subj in subjs:
    bda_stats[subj] = dict()
    for c in compars:
        bda_stats[subj][c[0]] = np.zeros(8)
        for i in range(8):
            try:
                _, bda_stats[subj][c[0]][i] = stats.ranksums(comb[subj][c[1]]['up'][:, i], comb[subj][c[2]]['up'][:, i])
            except IndexError:
                bda_stats[subj][c[0]][i] = 999  # data is missing from something
        # print(f'{subj}  {c[0]}\n {bda_stats[subj][c[0]]}\n')

# intra-subject bag type differences descending
"""
For the most part no differences descending from skimming numbers
"""
bdd_stats = dict()
for subj in subjs:
    bdd_stats[subj] = dict()
    for c in compars:
        bdd_stats[subj][c[0]] = np.zeros(8)
        for i in range(8):
            try:
                _, bdd_stats[subj][c[0]][i] = stats.ranksums(comb[subj][c[1]]['down'][:, i],
                                                             comb[subj][c[2]]['down'][:, i])
            except IndexError:
                bdd_stats[subj][c[0]][i] = 999  # data is missing from something
        # print(f'{subj}  {c[0]}\n {bdd_stats[subj][c[0]]}\n')

# fid = open('..\\Data\\stair_climb_stats_prelim.csv', 'w')
#
# fid.write('Subject,Direction,Comparison,Step_Length,Lateral_Dev,Step_Height,' +
#           'Max_Swing_Vel,Foot_Attack_Angle,Contact_Time,Step_Time,Cadence\n')
#
# for subj in subjs:
#     for i, v in zip([1, 2], [bda_stats, bdd_stats]):
#         for c in compars:
#             line = f'{int(subj)},{i},{c[0]}'
#             for j in range(8):
#                 line += f',{v[subj][c[0]][j]}'
#             line += '\n'
#             fid.write(line)
#
# fid.close()


# ANOVA

# Repeated measures ANOVA

def avg_stuff(data, subs, wd, m=0):
    savg = np.zeros((len(subs), 4))  # 4 different bag types

    rem = []  # indices to be removed if subj contains no data
    for sub, i in zip(subs, range(len(subs))):
        for bt, j in zip(['hb', 'lb', 'hs', 'ls'], range(4)):
            try:
                savg[i, j] = np.mean(data[sub][bt][wd][:, m])
            except IndexError:
                rem.append((i, j))

    for inds in rem:
        i, j = inds
        savg = np.delete(savg, i, axis=j)
    return savg


def rmANOVA(data):
    """
    Repeated Measures ANOVA

    data : array_like
        MxN array where M is the number of subjects, and N is the number of repeated measures
    """

    # TODO add sphericity check/corrections
    # TODO add normality check
    # TODO add option to not use all columns

    N = data.size  # number of elements
    s, a = data.shape  # number of subjects, number of repeated measures (time points)

    dfb = a-1  # Degrees of freedom (DoF) between
    dfw = N-a  # DoF within
    dfs = s-1  # DoF subjects
    dfe = dfw - dfs  # DoF error
    dft = N-1  # DoF total

    ssb = np.sum(np.sum(data, axis=0)**2)/s - np.sum(data)**2/N  # SS_between
    ssw = np.sum(data**2) - np.sum(np.sum(data, axis=0)**2)/s  # SS_within
    sss = np.sum(np.sum(data, axis=1)**2)/a - np.sum(data)**2/N  # SS_subjects
    sse = ssw - sss  # SS_error
    sst = np.sum(data**2) - np.sum(data)**2/N  # SS_total

    msb = ssb/dfb  # MS_between
    msw = ssw/dfw  # MS_within
    mss = sss/dfs  # MS_subject
    mse = sse/dfe  # MS_error
    mst = sst/dft  # MS_total

    F = msb/mse  # F-statistic

    p = stats.f.sf(F, dfb, dfe)

    eta_sq = ssb/(ssb + sse)

    # return (dfb, msb), (dfw, msw), (dfs, mss), (dfe, mse), (dft, mst)
    return F, p, eta_sq


metrics = ['step length', 'lateral deviation', 'step height', 'max swing velocity', 'foot attack angle', 'contact time',
           'step time', 'cadence']
for mm, met in zip(range(8), metrics):
    uavg = avg_stuff(comb, subjs, 'up', m=mm)
    davg = avg_stuff(comb, subjs, 'down', m=mm)

    Fu, pu, eu = rmANOVA(uavg)
    Fd, pd, ed = rmANOVA(davg)

    print(f'{met:19s}   Up: p={pu:.3f} eta^2={eu:.4f}   Down: p={pd:.3f} eta^2={ed:.4f}')
