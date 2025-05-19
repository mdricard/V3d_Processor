import numpy as np
import matplotlib.pyplot as plt
from Biomechanics import Biomechanics

fn = 'C:/Alexis_Subject_7/S7 new data 2.txt'
s_7 = Biomechanics(fn)
s_7.get_stance()

#s_27_n_08.plot_fz_steps()
s_7.plot_joint_force()
s_7.plot_joint_moment()