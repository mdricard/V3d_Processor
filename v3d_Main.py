import numpy as np
import matplotlib.pyplot as plt
from Biomechanics import Biomechanics

#  --------- Alexis input Subject #, mass and height ------------
mass = 73.482
height = 1.64592
subject = 8
#  --------------------------------------------------------------
path = 'D:/Alexis_Subject_' + str(subject) + '/S' + str(subject) + ' Data/'
for speed in range(3):
    for incline in range(3):
        for shoe in range(2):
            if speed == 0:
                speed_str = '05'
            elif speed == 1:
                speed_str = '08'
            else:
                speed_str = '12'
            if incline == 0:
                incline_str = 'Neutral'
            elif incline == 1:
                incline_str = 'UH'
            else:
                incline_str = 'DH'
            if shoe == 0:
                shoe_str = 'HK'
                filename = 'S' + str(subject) + ' ' + shoe_str + ' ' + incline_str + ' ' + speed_str + '.txt'
            else:
                shoe_str = ''
                filename = 'S' + str(subject) + shoe_str + ' ' + incline_str + ' ' + speed_str + '.txt'
            fn = path + filename
            #print('Speed ' + speed_str + ' Incline ' + incline_str + ' Shoe ' + shoe_str + ' File ' + filename)
            subj = Biomechanics(fn, subject, mass, height, speed, incline, shoe)
            subj.get_stance()
            subj.analyze_joint_force()
            subj.plot_joint_force()
            subj.plot_joint_moment()
            subj.save_stats_long()


"""

fn = 'D:/Alexis_Subject_' + str(subject) + '/S' + str(subject) + ' Data/'
s_7 = Biomechanics(fn, subject, mass, height, speed=0, incline=0, shoe=0)
s_7.get_stance()

s_7.analyze_joint_force()
#s_7.save_a_step(0)
#s_7.plot_right_fy()
#s_7.plot_shear_force()
#s_27_n_08.plot_fz_steps()
#s_7.plot_joint_force()
#s_7.plot_joint_moment()
#s_7.plot_joint_angle()
#s_7.plot_left_fy()
"""