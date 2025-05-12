import numpy as np
import matplotlib.pyplot as plt
from Biomechanics import Biomechanics

fn = 'D:/Alexis_Subject_27/LON_Rise.txt'
s_27_n_08 = Biomechanics(fn)
s_27_n_08.get_stance()

s_27_n_08.plot_fz()