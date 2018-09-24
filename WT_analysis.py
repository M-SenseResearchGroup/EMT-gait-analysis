# Walk and Turn.
# Code initially by Chris Petrillo
# Python Version: Lukas Adamowicz
# 7/25/18

import sys
sys.path.append("C:\\Users\\Lukas Adamowicz\\Dropbox\\Masters\\MC10py")
sys.path.append("C:\\Users\\Lukas Adamowicz\\Dropbox\\Masters\\Project - EMT Study\\Gait Analysis\\IMU_Gait_Analysis")
import MC10py
from IMUGaitParameters import GaitParameters


# def load_mc10(study_dir, segment=True, sync=True, save=True, save_loc=None, save_subj=False, return_data=True):
paths = MC10py.LoadMC10('C:\\Users\\Lukas Adamowicz\\Documents\\Study Data\\EMT\\ASYM_OFFICIAL', pre_time=0,
                        segment=True, sync=True, save=True, save_loc='import', save_subj=False, return_data=False)

raw_data = MC10py.OpenMC10('C:\\Users\\Lukas Adamowicz\\Documents\\Study Data\\EMT\\ASYM_OFFICIAL\\data.pickle')

gp = GaitParameters(raw_data)
gp.Run('..\\Data\\wt_results.csv')

