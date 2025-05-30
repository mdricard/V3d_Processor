import numpy as np
import matplotlib.pyplot as plt
from Biomechanics import Biomechanics
mass = 73.482
height = 1.64592
speed = 1.2
incline = 'Neutral'
subject = 8
shoe = 'barefoot'
fn = 'D:/Alexis_Subject_8/S8 Data/S8 Neutral 12.txt'
s_7 = Biomechanics(fn, subject, mass, height, speed, incline, shoe)
s_7.get_stance()
#s_7.plot_joint_angle()
#s_7.plot_left_fy()
s_7.analyze_joint_force()
s_7.save_a_step(0)
#s_7.plot_right_fy()
#s_7.plot_shear_force()
#s_27_n_08.plot_fz_steps()
#s_7.plot_joint_force()
#s_7.plot_joint_moment()