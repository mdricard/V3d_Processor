import numpy as np
import matplotlib.pyplot as plt
from Biomechanics import Biomechanics

fn = 'D:/Alexis_Subject_27/LON_Fall.txt'
s_27_n_08 = Biomechanics(fn)
s_27_n_08.get_stance()

#s_27_n_08.plot_fz_steps()
s_27_n_08.plot_joint_force()
s_27_n_08.plot_joint_moment()