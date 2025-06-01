import numpy as np
import matplotlib.pyplot as plt
from math import degrees, asin
from BiomechTools import low_pass, zero_crossing, max_min, simpsons_rule, critically_damped, residual_analysis, get_max_value, get_min_value

class Biomechanics:
    RON = np.zeros(60, dtype=int)
    ROFF = np.zeros(60, dtype=int)
    LON = np.zeros(60, dtype=int)
    LMID = np.zeros(60, dtype=int)  # start of propulsion for left leg
    RMID = np.zeros(60, dtype=int)
    LOFF = np.zeros(60, dtype=int)
    RMID = np.zeros(60, dtype=int)
    LOFF = np.zeros(60, dtype=int)
    speed_str = ''
    incline_str = ''
    shoe_str = ''
    height = 0
    speed = 0
    incline = 0
    shoe = 0
    peak_comp = np.zeros(60)
    peak_comp_pt = np.zeros(60, dtype=int)
    peak_shear = np.zeros(60)
    peak_shear_pt = np.zeros(60, dtype=int)
    comp_impulse = np.zeros(60)
    shear_impulse = np.zeros(60)
    peak_add = np.zeros(60)
    peak_add_pt = np.zeros(60, dtype=int)
    add_impulse = np.zeros(60)
    trail_leg_prop = np.zeros(60)
    lead_leg_braking = np.zeros(60)
    hip_ron = np.zeros(60)
    knee_ron = np.zeros(60)
    knee_flex_range = np.zeros(60)
    n_reps = 0
    subject = ''
    mass = 0
    n_steps = 0
    n_vars = 0
    var_name = []
    v3d = []
    n_rows = 0
    n_cols = 0
    FP1_X = []
    FP1_Y = []

    def __init__(self, filename, subject, mass, height, speed, incline, shoe):
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
            self.subject = subject
            self.mass = mass
            self.height = height
            self.speed = speed
            self.incline = incline
            self.shoe = shoe
        data = np.genfromtxt(filename, delimiter='\t', skip_header=7, skip_footer=2)
        self.n_rows = data.shape[0]  # number of rows of array
        self.n_cols = data.shape[1]
        self.FP1_X = data[:, 1]
        self.FP1_Y = -1.0 * data[:, 2] / (self.mass * 9.8)     # convert to BW
        self.FP1_Z = data[:, 3]
        self.FP2_X = data[:, 4]
        self.FP2_Y = -1.0 * data[:, 5] / (self.mass * 9.8)     # convert to BW
        self.FP2_Z = data[:, 6]
        self.R_HIP_Jt_Angle_X = data[:, 7]
        self.R_HIP_Jt_Angle_Y = data[:, 8]
        self.R_HIP_Jt_Angle_Z = data[:, 9]
        self.Rt_Knee_Jt_Angle_X = data[:, 10]
        self.Rt_Knee_Jt_Angle_Y = data[:, 11]
        self.Rt_Knee_Jt_Angle_Z = data[:, 12]
        self.Rt_Knee_Jt_Force_X = data[:, 13] / (self.mass * 9.8)     # convert to BW
        self.Rt_Knee_Jt_Force_Y = data[:, 14] / (self.mass * 9.8)     # convert to BW
        self.Rt_Knee_Jt_Force_Z = data[:, 15] / (self.mass * 9.8)     # convert to BW
        self.Rt_Knee_Jt_Moment_X = data[:, 16]
        self.Rt_Knee_Jt_Moment_Y = data[:, 17]
        self.Rt_Knee_Jt_Moment_Z = data[:, 18]

        # smooth Forces at 20 Hz
        self.FP1_X = critically_damped(self.FP1_X, 1000, 20)
        self.FP1_Y = critically_damped(self.FP1_Y, 1000, 20)
        self.FP1_Z = critically_damped(self.FP1_Z, 1000, 20)
        self.FP2_X = critically_damped(self.FP2_X, 1000, 20)
        self.FP2_Y = critically_damped(self.FP2_Y, 1000, 20)
        self.FP2_Z = critically_damped(self.FP2_Z, 1000, 20)

        #for i in range(1, self.n_cols):
        #    print(self.var_name[i])
        #self.n_steps = int(
        #    (self.n_cols - 1) / 3 / self.n_vars)
        # drop first col, divide by 3 (x,y,z) then divide by # variables in file
        #self.v3d = pd.DataFrame(data, columns=self.var_name)

    def get_stance(self):
        Lf, Lf_rf = zero_crossing(self.FP1_Z, 20, 0, self.n_rows - 1)
        Rt, Rt_rf = zero_crossing(self.FP2_Z, 20, 0, self.n_rows - 1)
        # Find the First Step  LON must be less than RON
        last_pt = len(Rt) - 1
        if Rt_rf[last_pt] == 'rising':
            last_pt = last_pt - 1
        i = 0
        step = 0
        while (Lf[i] > Rt[i]) and (Lf_rf[i] != 'rising'):
            i = i + 1
        while i < last_pt - 2:
            self.LON[step] = Lf[i]
            self.LOFF[step] = Lf[i + 1]
            if Rt_rf[i + 1] == 'rising':
                self.RON[step] = Rt[i + 1]
                self.ROFF[step] = Rt[i + 2]
            else:
                self.RON[step] = Rt[i + 2]
                self.ROFF[step] = Rt[i + 3]
            i = i + 2
            step = step + 1
            self.n_steps = step

        # Find the start of propulsion for left leg
        for step in range(self.n_steps):
            k = self.LON[step] + 100    # skip first 100 pts
            while self.FP1_Y[k] < 0.0 and k < self.LOFF[step]:
                k = k + 1
            self.LMID[step] = k

        # Find the start of propulsion for right leg
        for step in range(self.n_steps):
            k = self.RON[step] + 100    # skip first 100 pts
            while self.FP2_Y[k] < 0.0 and k < self.ROFF[step]:
                k = k + 1
            self.RMID[step] = k

    def plot_first_step(self):
        plt.plot(self.FP2_Z[Lf[0]:Rt[1]], 'r', label='FP2 Z')
        plt.plot(self.FP1_Z, 'b', label='FP1 Z')
        plt.grid(True)
        plt.legend()
        plt.show()

    def plot_fz_steps(self):
        for i in range(self.n_steps):
            plt.plot(self.FP1_Z[int(self.LON[i]):int(self.LOFF[i])], 'b', label='FP1 Z')
            plt.plot(self.FP2_Z[int(self.RON[i]):int(self.ROFF[i])], 'r', label='FP2 Z')
            plt.grid(True)
            plt.title('Vertical Force Step ' + str(i))
            plt.legend()
            plt.show()

    def plot_left_fy(self):
        for i in range(self.n_steps):
            plt.plot(self.FP1_Y[int(self.LON[i]):int(self.LOFF[i])], 'b', label='FP1 Y')
            plt.vlines(self.LMID[i] - self.LON[i], -.4, .4, colors='red', linestyles='dotted')
            #plt.plot(self.FP2_Z[int(self.RON[i]):int(self.ROFF[i])], 'r', label='FP2 Z')
            plt.grid(True)
            plt.title('Left Leg A/P Force (BW) ' + str(i))
            print(self.LON[i], self.LMID[i], self.LOFF[i])
            plt.legend()
            plt.show()

    def plot_right_fy(self):
        for i in range(self.n_steps):
            plt.plot(self.FP2_Y[int(self.RON[i]):int(self.ROFF[i])], 'b', label='FP2 Y')
            plt.vlines(self.RMID[i] - self.RON[i], -.4, .4, colors='red', linestyles='dotted')
            #plt.plot(self.FP2_Z[int(self.RON[i]):int(self.ROFF[i])], 'r', label='FP2 Z')
            plt.grid(True)
            plt.title('Right Leg A/P Force (BW) ' + str(i))
            print(self.RON[i], self.RMID[i], self.ROFF[i])
            plt.legend()
            plt.show()

    def get_plot_titletext(self):
        if self.speed == 0:
            self.speed_str = '05 '
        elif self.speed == 1:
            self.speed_str = '08 '
        else:
            self.speed_str = '12 '
        if self.incline == 0:
            self.incline_str = 'Neutral '
        elif self.incline == 1:
            self.incline_str = 'UH '
        else:
            self.incline_str = 'DH '
        if self.shoe == 0:
            self.shoe_str = 'HK '
        else:
            self.shoe_str = 'Barefoot '


    def plot_joint_force(self):
        self.get_plot_titletext()
        for i in range(self.n_steps):
            plt.plot(self.Rt_Knee_Jt_Force_Z[int(self.RON[i]):int(self.ROFF[i])], label='Step ' + str(i))
            plt.grid(True)
            plt.title('Subject ' + str(self.subject) + ' ' + self.shoe_str + self.incline_str + ' ' + self.speed_str + ' Rt Knee Comp Force')
            #plt.vlines(self.peak_comp_pt[i], -1.0, 0, colors='red', linestyles='dotted')
            plt.legend()
        plt.show()


    def plot_shear_force(self):
        self.get_plot_titletext()
        for i in range(self.n_steps):
            plt.plot(self.Rt_Knee_Jt_Force_Y[int(self.RON[i]):int(self.ROFF[i])], label='Step ' + str(i))
            plt.grid(True)
            plt.title('Subject ' + str(self.subject) + ' ' + self.shoe_str + self.incline_str + ' ' + self.speed_str + ' Rt Knee Shear Force')
            plt.legend()
        plt.show()


    def plot_joint_moment(self):
        self.get_plot_titletext()
        for i in range(self.n_steps):
            plt.plot(self.Rt_Knee_Jt_Moment_Y[int(self.RON[i]):int(self.ROFF[i])], label='Step ' + str(i))
            plt.grid(True)
            plt.title('Subject ' + str(self.subject) + ' ' + self.shoe_str + self.incline_str + ' ' + self.speed_str  + ' Rt Knee Adduction Moment')
            plt.legend()
        plt.show()

    def plot_joint_angle(self):
        for i in range(self.n_steps):
            plt.plot(self.Rt_Knee_Jt_Angle_X[self.RON[i]:self.ROFF[i]], label='Step ' + str(i))
            plt.grid(True)
            plt.title('Subject ' + str(self.subject) + ' ' + self.incline + ' ' + str(self.speed) + ' Rt Knee Angle')
            plt.legend()
        plt.show()

    def plot_fz(self):
        plt.plot(self.FP2_Z, 'r', label='FP2 Z')
        plt.plot(self.FP1_Z, 'b', label='FP1 Z')
        plt.grid(True)
        plt.legend()
        plt.show()

    def save_a_step(self, step):
        f_step_name = 'd:/step_' + str(step) + '.csv'
        a = self.FP1_Y[self.LON[step] : self.ROFF[step]]
        b = self.FP1_Z[self.LON[step] : self.ROFF[step]]
        c = self.FP2_Y[self.LON[step] : self.ROFF[step]]
        d = self.FP2_Z[self.LON[step] : self.ROFF[step]]
        e = self.Rt_Knee_Jt_Force_Y[self.LON[step] : self.ROFF[step]]
        f = self.Rt_Knee_Jt_Force_Z[self.LON[step] : self.ROFF[step]]
        g = self.Rt_Knee_Jt_Moment_Y[self.LON[step]: self.ROFF[step]]

        np.savetxt(f_step_name, np.column_stack((a, b, c, d, e, f, g)),  fmt='%.6f', delimiter=',', newline='\n', header="fp1 Y, fp1 Z, fp2 Y, fp2 Z, Rt Knee Force Y, Rt Knee Force Z, Rt Knee Adduction Moment Y", comments="")


    def analyze_joint_force(self):
        for i in range(self.n_steps):
            self.peak_comp[i], self.peak_comp_pt[i] = get_min_value(self.Rt_Knee_Jt_Force_Z, self.RON[i], self.RMID[i])
            self.comp_impulse[i] = simpsons_rule(self.Rt_Knee_Jt_Force_Z, self.RON[i], self.peak_comp_pt[i], 1.0)
            self.peak_add[i], self.peak_add_pt[i] = get_min_value(self.Rt_Knee_Jt_Moment_Y, self.RON[i], self.RMID[i])
            self.add_impulse[i] = simpsons_rule(self.Rt_Knee_Jt_Moment_Y, self.RON[i], self.peak_add_pt[i], 0.001)
            self.peak_shear[i], self.peak_shear_pt[i] = get_max_value(self.Rt_Knee_Jt_Force_Y, self.RON[i], self.RMID[i])
            self.shear_impulse[i] = simpsons_rule(self.Rt_Knee_Jt_Force_Y, self.RON[i], self.peak_shear_pt[i], 1.0)
            self.trail_leg_prop[i] = simpsons_rule(self.FP1_Y, self.LMID[i], self.LOFF[i], 1.0)
            self.lead_leg_braking[i] = simpsons_rule(self.FP2_Y, self.RON[i], self.peak_comp_pt[i], 1.0)
            self.hip_ron[i] = self.R_HIP_Jt_Angle_X[self.RON[i]]
            self.knee_ron[i] = self.Rt_Knee_Jt_Angle_X[self.RON[i]]
            self.knee_flex_range[i] = self.Rt_Knee_Jt_Angle_X[self.peak_comp_pt[i]] - self.knee_ron[i] - self.hip_ron[i]

    def save_stats_long(self):
        stat_file_path = 'D:/Alexis_Stats/'
        fn = stat_file_path + 'Alexis_Stats.csv'
        with open(fn, 'a') as stat_file:
            stat_file.write('subject, shoe,speed, incline, step, comp_force,comp_impulse,shear_force,shear_impulse,add_mom, add_impulse, trail_prop, lead_braking, hip_ron, knee_ron\n')
            for step in range(self.n_steps):
                stat_file.write(
                    str(self.subject) + ',' + str(self.shoe)  + ',' + str(self.speed) + ',' + str(self.incline) + ',' + str(step) + ',' + str(self.peak_comp[step]) + ',' + str(self.comp_impulse[step]) + ',' + str(self.peak_shear[step]) + ',' + str(self.shear_impulse[step]) + ',' + str(self.peak_add[step]) + ',' + str(self.add_impulse[step]) + ',' + str(self.trail_leg_prop[step]) + ',' + str(self.lead_leg_braking[step]) + ',' + str(self.hip_ron[step]) + ',' + str(self.knee_ron[step]) + '\n')
        stat_file.close()
        print('Stats saved to ' + fn)
