import numpy as np
import matplotlib.pyplot as plt
from Biomechanics import Biomechanics
mass = 73.482
height = 1.64592
speed = 0.8
incline = 'neutral'
fn = 'D:/Alexis_Subject_8/S8 Data/S8 DH 05.txt'
s_7 = Biomechanics(fn, mass, height, speed, incline)
s_7.get_stance()

#s_27_n_08.plot_fz_steps()
s_7.plot_joint_force()
s_7.plot_joint_moment()