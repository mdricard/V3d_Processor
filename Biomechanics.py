import numpy as np
import matplotlib.pyplot as plt
from math import degrees, asin
from BiomechTools import low_pass, zero_crossing, max_min, simpson_nonuniform, critically_damped, residual_analysis

class Biomechanics:
    RON = np.zeros(50)
    ROFF = np.zeros(50)
    LON = np.zeros(50)
    LOFF = np.zeros(50)

    def __init__(self, filename):
        with open(filename) as f:
            f.readline()  # skip first line
            self.var_name = f.readline().rstrip().split()
            self.n_vars = len(set(self.var_name))  # use set to count # of unique values in a list
            self.var_name.insert(0, 'pt_num')
            f.readline()  # skip next line
            f.readline()  # skip next line
            xyz = f.readline().rstrip().split()
            for i in range(0, len(self.var_name)):
                self.var_name[i] = self.var_name[i] + '_' + xyz[i]
            self.var_name[0] = 'pt_num'

            data = np.genfromtxt(filename, delimiter='\t', skip_header=5)
        self.n_rows = data.shape[0]  # number of rows of array
        self.n_cols = data.shape[1]
        self.FP1_X = data[:, 1]
        self.FP1_Y = data[:, 2]
        self.FP1_Z = data[:, 3]
        self.FP2_X = data[:, 4]
        self.FP2_Y = data[:, 5]
        self.FP2_Z = data[:, 6]
        self.R_HIP_Jt_Angle_X = data[:, 7]
        self.R_HIP_Jt_Angle_Y = data[:, 8]
        self.R_HIP_Jt_Angle_Z = data[:, 9]
        self.Rt_Knee_Jt_Angle_X = data[:, 10]
        self.Rt_Knee_Jt_Angle_Y = data[:, 11]
        self.Rt_Knee_Jt_Angle_Z = data[:, 12]
        self.Rt_Knee_Jt_Force_X = data[:, 13]
        self.Rt_Knee_Jt_Force_Y = data[:, 14]
        self.Rt_Knee_Jt_Force_Z = data[:, 15]
        self.Rt_Knee_Jt_Moment_X = data[:, 16]
        self.Rt_Knee_Jt_Moment_Y = data[:, 17]
        self.Rt_Knee_Jt_Moment_Z = data[:, 18]

        #for i in range(1, self.n_cols):
        #    print(self.var_name[i])
        #self.n_steps = int(
        #    (self.n_cols - 1) / 3 / self.n_vars)  # drop first col, divide by 3 (x,y,z) then divide by # variables in file
        #self.v3d = pd.DataFrame(data, columns=self.var_name)

    def get_stance(self):
        Lf, Lf_rf = zero_crossing(self.FP1_Z, 16, 0, self.n_rows-1)
        Rt, Rt_rf = zero_crossing(self.FP2_Z, 16, 0, self.n_rows-1)
        i = 0
        while (Lf[i] > Rt[i]) or (Lf_rf[i] != 'rising'):
            i = i + 1
        Lf = Lf[i]
        Rt = Rt[i]
        print(
            'LF: ', Lf,
            'RT: ', Rt
        )
        if Lf[0] < Rt[0]:
            print('LON ', Lf[0])
            print('RON', Rt[0])

    def plot_first_step(self):
        plt.plot(self.FP2_Z[Lf[0]:Rt[1]], 'r', label='FP2 Z')
        plt.plot(self.FP1_Z, 'b', label='FP1 Z')
        plt.grid(True)
        plt.legend()
        plt.show()

    def plot_fz(self):
        plt.plot(self.FP2_Z, 'r', label='FP2 Z')
        plt.plot(self.FP1_Z, 'b', label='FP1 Z')
        plt.grid(True)
        plt.legend()
        plt.show()
