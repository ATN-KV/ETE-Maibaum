from sre_constants import SUCCESS
from typing import Tuple
import matplotlib as mlt
import matplotlib.pyplot as plt
from brokenaxes import brokenaxes
import numpy as np
import pandas as pd
from matplotlib.colors import LightSource
from mpl_toolkits import mplot3d
import math
import scipy as sp
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from MaibaumPyth_force import func_Force_over_LengthRope as func_force
from MaibaumPyth_force import l_St, lr_min, lr_max, nan, katnr

# global variables
res = 1000

# motor data

motor_type = "G4-M8-encoder"
j_motor = (4.5 / (100**2))
m_1 = 26.5
m_2 = 13
omega_1 = 1750 * (np.pi / 30)
omega_2 = 4700 * (np.pi / 30)

# motor_type = "G4-M6-encoder"
# j_motor = (2.4 / (100**2))
# m_1 = 13.2
# m_2 = 6.0
# omega_1 = 2000 * (np.pi / 30)
# omega_2 = 5800 * (np.pi / 30)

# motor_type = "G5-M2-encoder"
# j_motor = (4.4 / (100**2))
# m_1 = 13.5
# m_2 = 6.0
# omega_1 = 2200 * (np.pi / 30)
# omega_2 = 6800 * (np.pi / 30)

# motor_type = "G5-M4-encoder"
# j_motor = (7.8 / (100**2))
# m_1 = 27
# m_2 = 13.5
# omega_1 = 1350 * (np.pi / 30)
# omega_2 = 4200 * (np.pi / 30)

# setup variables 
j_gearbox = (8 / (100**2))
d_hoistdrum = 0.2
d_hoistdrum_gearshaft = 0.05
i_g_ud = 150
etha_g = 0.83
a_ud = 14
t_max = 75 + (katnr * np.pi)
t_max_def = 150
m_winch = 40


def main():

    fig_1_GridRows = 2
    fig_1_GridCols = 2

    #### 2D Plot (higher resolution)
    a = np.linspace(l_St, 100, 90)
    i_g = np.linspace(144, 148, 9)

    # set global plt.parameters
    parameters = {"xtick.labelsize":8}
    plt.rcParams.update(parameters)

    # create figure
    fig_1 = plt.figure(figsize=(18, 9), constrained_layout=True)
    fig_2 = plt.figure(figsize=(18, 9), constrained_layout=True)
    fig_3 = plt.figure(figsize=(18, 9), constrained_layout=True)


    # calculations iterating
    # time_to_finish_forloop = np.inf
    # a_forloop = nan
    # i_g_forloop = nan
    # for var_a in a:
    #     ti = np.inf
    #     print("#a#", var_a, "#ig#", i_g_forloop, "#time#", time_to_finish_forloop)
    #     for var_i_g in i_g:
    #         dict_forloop = func_calculations(par_a=var_a, par_i_g=var_i_g, par_res=res, par_M1=m_1, par_M2=m_2, 
    #                                             par_Omega1=omega_1, par_Omega2=omega_2, par_t_max=150, 
    #                                             par_etha_g=etha_g, par_j_gearbox=j_gearbox, par_j_motor=j_motor)
    #         if dict_forloop["success"] == True:
    #             if (dict_forloop["time_to_finish"] < time_to_finish_forloop) and not np.isnan(dict_forloop["time_to_finish"]):
    #                 time_to_finish_forloop = dict_forloop["time_to_finish"]
    #                 a_forloop = var_a
    #                 i_g_forloop = var_i_g
    #             if ti < dict_forloop["time_to_finish"]:
    #                 break
    #             ti = dict_forloop["time_to_finish"]

    # a_ud = a_forloop
    # i_g_ud = i_g_forloop

    # # calculations iterating a = const
    # time_to_finish_forloop = np.inf
    # i_g_forloop = nan
    # for var_i_g in i_g:
    #     print(var_i_g," # ", end="")
    #     dict_forloop = func_calculations(par_a=a_ud, par_i_g=var_i_g, par_res=res, par_M1=m_1, par_M2=m_2, 
    #                                         par_Omega1=omega_1, par_Omega2=omega_2, par_t_max=100, 
    #                                         par_etha_g=etha_g, par_j_gearbox=j_gearbox, par_j_motor=j_motor)
    #     if dict_forloop["success"] == False:
    #         continue
    #     else:
    #         print(dict_forloop["time_to_finish"])
    #         if (dict_forloop["time_to_finish"] < time_to_finish_forloop) and not np.isnan(dict_forloop["time_to_finish"]) and (dict_forloop["time_to_finish"] != 0):
    #             time_to_finish_forloop = dict_forloop["time_to_finish"]
    #             i_g_forloop = var_i_g

    # i_g_ud = i_g_forloop

    # calculation 1 condition
    dict_sc_1 = func_calculations(par_a=a_ud, par_i_g=i_g_ud, par_res=2000, par_M1=m_1, par_M2=m_2, 
                                    par_Omega1=omega_1, par_Omega2=omega_2, par_t_max=150, 
                                    par_etha_g=etha_g, par_j_gearbox=j_gearbox, par_j_motor=j_motor)

    # printing results
    if dict_sc_1["success"] == False:
        print("\n\n\n\n\nResults of the calculation:\n")
        print(f"Motor used: {motor_type}")
        print("Motor is to weak")
        print(f"distance a: {a_ud} m\ngear ratio: 1 : {i_g_ud}")
        print("\n\n\n\n\n")
    else:
        print("\n\n\n\n\nResults of the calculation:\n")
        print(f"Motor used: {motor_type}")
        print("time till finish: {}s" .format(dict_sc_1["time_to_finish"]))
        print("a_exp: {}" .format(round(dict_sc_1["a_exp"], 6)))
        print("b_exp: {}" .format(round(dict_sc_1["b_exp"], 6)))
        print(f"distanc a: {round(a_ud, 3)} m")
        print(f"gearbox ratio: 1 : {i_g_ud}")
        print(f"gearbox efficiency: {etha_g}")
        print(f"diameter hoistdrum: {d_hoistdrum} m")
        print(f"J gearbox: {j_gearbox} kg×m²")
        print("J collective: {} kg×m²" .format(round(dict_sc_1["j"], 6)))
        print("\n\n\n\n\n")

        # plotting results

        # plotting phi
        ax_1 = fig_1.add_subplot(fig_1_GridRows, fig_1_GridCols, 1)
        # axes labeling and title
        ax_1.set_xlabel("time in [s]")
        ax_1.set_ylabel("phi motor in [turns]")
        ax_1.set_title(f"phi motor at (a = {round(a_ud, 3)}m)")
        #### Set the X - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_1.xaxis.set_major_locator(loc)
        #### Set the Y - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=250.0)
        ax_1.yaxis.set_major_locator(loc)
        # setting grid options
        ax_1.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
        # plot
        ax_1.plot(dict_sc_1["t"], (dict_sc_1["phi_ot"]/(2*np.pi)), color="g")

        # plotting omega
        ax_2 = fig_1.add_subplot(fig_1_GridRows, fig_1_GridCols, 2)
        # axes labeling and title
        ax_2.set_xlabel("time in [s]")
        ax_2.set_ylabel("omega motor in [rad/s]")
        ax_2.set_title(f"omega motor at (a = {round(a_ud, 3)}m)")
        #### Set the X - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_2.xaxis.set_major_locator(loc)
        #### Set the Y - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=25.0)
        ax_2.yaxis.set_major_locator(loc)
        # setting grid options
        ax_2.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
        # plot
        ax_2.plot(dict_sc_1["t"], dict_sc_1["omega_ot"], color="r")

        # plotting force over time
        ax_3 = fig_1.add_subplot(fig_1_GridRows, fig_1_GridCols, 3)
        # axes labeling and title
        ax_3.set_xlabel("time in [s]")
        ax_3.set_ylabel("force at winch in [kN]")
        ax_3.set_title(f"force at winch at (a = {round(a_ud, 3)}m)")
        #### Set the X - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_3.xaxis.set_major_locator(loc)
        #### Set the Y - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_3.yaxis.set_major_locator(loc)
        # setting grid options
        ax_3.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
        # plot
        ax_3.plot(dict_sc_1["t"], dict_sc_1["f_winch"]/1000, color="#EA8000")

        # rope length over time
        ax_4 = fig_1.add_subplot(fig_1_GridRows, fig_1_GridCols, 4)
        # axes labeling and title
        ax_4.set_xlabel("time in [s]")
        ax_4.set_ylabel("length of rope in [m]")
        ax_4.set_title(f"rope length at (a = {round(a_ud, 3)}m)")
        #### Set the X - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_4.xaxis.set_major_locator(loc)
        #### Set the Y - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_4.yaxis.set_major_locator(loc)
        # setting grid options
        ax_4.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
        # plot
        ax_4.plot(dict_sc_1["t"], dict_sc_1["rope_length"], color="#E400FF")

        # plotting alpha
        ax_1 = fig_2.add_subplot(fig_1_GridRows, fig_1_GridCols, 3)
        # axes labeling and title
        ax_1.set_xlabel("time in [s]")
        ax_1.set_ylabel("alpha motor in [rad/s²]")
        #### Set the X - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_1.xaxis.set_major_locator(loc)
        #### Set the Y - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_1.yaxis.set_major_locator(loc)
        # setting grid options
        ax_1.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
        # set y lim
        ax_1.set_ylim(-20, 25)
        # hide spines
        # plot
        ax_1.plot(dict_sc_1["t"], dict_sc_1["alpha"], color="#00F3FF")
        ax_11 = fig_2.add_subplot(fig_1_GridRows, fig_1_GridCols, 1)
        # axes labeling and title
        ax_11.set_ylabel("alpha motor in [rad/s²]")
        ax_11.set_title(f"alpha motor at (a = {round(a_ud, 3)}m)")
        #### Set the X - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_11.xaxis.set_major_locator(loc)
        #### Set the Y - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_11.yaxis.set_major_locator(loc)
        # setting grid options
        ax_11.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
        # set y lim
        ax_11.set_ylim(2030, 2075)
        # plot
        ax_11.plot(dict_sc_1["t"], dict_sc_1["alpha"], color="#00F3FF")

        # plotting accel. torque
        ax_3 = fig_2.add_subplot(fig_1_GridRows, fig_1_GridCols, 2)
        # axes labeling and title
        ax_3.set_xlabel("time in [s]")
        ax_3.set_ylabel("accel torque motor in [Nm]")
        ax_3.set_title(f"accel torque at motor at (a = {round(a_ud, 3)}m)")
        #### Set the X - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_3.xaxis.set_major_locator(loc)
        #### Set the Y - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=0.0025)
        ax_3.yaxis.set_major_locator(loc)
        # setting grid options
        ax_3.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
        # set y lim
        ax_3.set_ylim(2.545, 2.595)
        # hide spines
        # plot
        ax_3.plot(dict_sc_1["t"], (dict_sc_1["alpha"] * dict_sc_1["j"]), color="#00FF2A")
        ax_33 = fig_2.add_subplot(fig_1_GridRows, fig_1_GridCols, 4)
        # axes labeling and title
        ax_33.set_xlabel("time in [s]")
        ax_33.set_ylabel("accel torque motor in [Nm]")
        #### Set the X - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_33.xaxis.set_major_locator(loc)
        #### Set the Y - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=0.0025)
        ax_33.yaxis.set_major_locator(loc)
        # setting grid options
        ax_33.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
        # set y lim
        ax_33.set_ylim(-0.025, 0.03)
        # hide spines
        # plot
        ax_33.plot(dict_sc_1["t"], (dict_sc_1["alpha"] * dict_sc_1["j"]), color="#00FF2A")

        # plotting Mfmot
        ax_1 = fig_3.add_subplot(1, 2, 1)
        # axes labeling and title
        ax_1.set_xlabel("time in [s]")
        ax_1.set_ylabel("torque tree at motor in [Nm] (M_F_Mot)")
        ax_1.set_title(f"torque tree at motor at (a = {round(a_ud, 3)}m)")
        #### Set the X - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_1.xaxis.set_major_locator(loc)
        #### Set the Y - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=0.5)
        ax_1.yaxis.set_major_locator(loc)
        # setting grid options
        ax_1.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
        # set y lim
        #ax_1.set_ylim(-20, 25)
        # hide spines
        # plot
        ax_1.plot(dict_sc_1["t"], dict_sc_1["torque_f_mot"], color="#00F3FF")

        # plotting Mmot
        ax_2 = fig_3.add_subplot(1, 2, 2)
        # axes labeling and title
        ax_2.set_xlabel("time in [s]")
        ax_2.set_ylabel("torque at motor in [Nm] (M_Mot)")
        ax_2.set_title(f"torque at motor at (a = {round(a_ud, 3)}m)")
        #### Set the X - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=2.0)
        ax_2.xaxis.set_major_locator(loc)
        #### Set the Y - marker stepsize
        loc = mlt.ticker.MultipleLocator(base=0.5)
        ax_2.yaxis.set_major_locator(loc)
        # setting grid options
        ax_2.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
        # set y lim
        #ax_2.set_ylim(-20, 25)
        # hide spines
        # plot
        ax_2.plot(dict_sc_1["t"], (dict_sc_1["torque_f_mot"] + dict_sc_1["alpha"] * dict_sc_1["j"]), color="#FF0000")


        plt.show()




# ==========================
    # ==========================
    # functions
    # ==========================
    # ==========================
def func_dSdx(t, S, par_a=0.0, par_i_g=0.0, par_j=0.0, par_a_exp=0.0, par_b_exp=0.0, 
                par_lr_start=0.0, par_M1=0.0, par_Omega1=0.0, par_etha_g=0.0):
    x, v = S
    func = 0
    term_1 = 0
    term_2 = 0

    term_1 = (-((d_hoistdrum/2) * (1/(par_i_g*par_etha_g)) * func_force(par_a, par_lr_start-(d_hoistdrum*(x/(2*par_i_g)))))) / par_j
    term_2 = ( ( ((par_a_exp*(par_b_exp**v))-par_M1) * (np.sign(v-par_Omega1)+1) * (1/2) ) + par_M1) / par_j
    func = term_1 + term_2

    return [v, func]

def func_calculations(par_a=0.0, par_i_g=0.0, par_res=0.0, par_M1=0.0, par_M2=0.0, par_Omega1=0.0, par_Omega2=0.0, 
                        par_t_max=0.0, par_etha_g=0.0, par_j_gearbox=0, par_j_motor=0):
    # check for max force: is the motor strong enough with the i_g and etha_g?
    # creating values for the axis
    force_max_motor = (par_M1*par_i_g*par_etha_g*2)/d_hoistdrum
    lr = np.linspace(lr_min, lr_max, par_res)
    for var_lr in lr:
        x = func_force(par_a, var_lr)
        if x > force_max_motor:
            dict_return = {
                "success" : False,
                "time_to_finish" : nan
            }
            return dict_return

    time_to_finish = nan
    # initializing force
    F_winch = np.linspace(0, 1, par_res)
    F_winch[:] = nan
    # initializing rope length
    rope_length = np.linspace(0, 1, par_res)
    rope_length[:] = nan

    # calc rotation inertia of the whole setup
    j_winch = m_winch * ((((d_hoistdrum/2)**2) + ((d_hoistdrum_gearshaft/2)**2)) / 2)
    j_red = j_winch * ((1/par_i_g)**2) + par_j_gearbox
    j = par_j_motor + j_red

    # calc lr_start
    lr_start = np.sqrt((par_a-np.cos(np.radians(60))*l_St)**2+(np.sin(np.radians(60))*l_St)**2) + 2 * l_St

    # calc of exp-funct.
    b_exp = np.e**((np.log(par_M1) - np.log(par_M2)) / (par_Omega1 - par_Omega2))
    a_exp = par_M1 / (b_exp**par_Omega1)

    # defining initial conditions
    x_0 = 0
    v_0 = 0
    S_0 = (x_0, v_0)

    # solving of ODE
    t = np.linspace(0, par_t_max, par_res)
    sol = odeint(func_dSdx, y0=S_0, t=t, tfirst=True, args=(par_a, par_i_g, j, a_exp, b_exp, lr_start, par_M1, par_Omega1, par_etha_g))

    # assigning results of ODE
    phi_ot = sol.T[0]
    omega_ot = sol.T[1]

    # initializing alpha
    alpha = np.linspace(0,1,len(omega_ot))
    alpha[:] = nan
    
    # calculation of terms for alpha
    term_1 = ((d_hoistdrum/2) * (1/(par_i_g*par_etha_g)) * func_force(par_a, lr_start-(d_hoistdrum*(phi_ot/(2*par_i_g))))) / j
    term_2 = ( ( ((a_exp*(b_exp**omega_ot))-par_M1) * (np.sign(omega_ot-par_Omega1)+1) * (1/2) ) + par_M1) / j

    # changing ODE results due to max omega
    i = 0
    for i in range(len(phi_ot)):
        # if: if omega is bigger than omega max (omega = const; alpha = 0)
        if (omega_ot[i] > par_Omega2):
            omega_ot[i] = par_Omega2
            phi_ot[i] = phi_ot[i-1] + par_Omega2 * (par_t_max / par_res)
            alpha[i] = 0
        # elif: if force is not zero but omega reached end
        elif np.isnan(phi_ot[i]) and not np.isnan(func_force(par_a, lr_start-(d_hoistdrum*((phi_ot[i-1] + par_Omega2 * (par_t_max / par_res))/(2*par_i_g))))):
            omega_ot[i] = par_Omega2
            phi_ot[i] = phi_ot[i-1] + par_Omega2 * (par_t_max / par_res)
            alpha[i] = 0
            if np.isnan(func_force(par_a, lr_start-(d_hoistdrum*((phi_ot[i] + par_Omega2 * (par_t_max / par_res))/(2*par_i_g))))):
                omega_ot[i] = nan
                phi_ot[i] = nan
        elif (omega_ot[i] <= par_Omega2):
            alpha = term_2 - term_1


    # calculating motor torque
    torque_f_mot = (d_hoistdrum/2) * (1/(par_i_g*par_etha_g)) * func_force(par_a, lr_start-(d_hoistdrum*(phi_ot/(2*par_i_g))))

    # calc force over time
    i = 0
    for x in phi_ot:
        F_winch[i] = func_force(par_a, lr_start-(d_hoistdrum*(x/(2*par_i_g))))
        i += 1

    # calc rope length over time and time to finish
    i = 0
    toggle_timefinder = True
    for x in phi_ot:
        rope_length[i] = lr_start-(d_hoistdrum*(x/(2*par_i_g)))
        if toggle_timefinder and np.isnan(phi_ot[i]):
            toggle_timefinder = False
            time_to_finish = i * (par_t_max/par_res)
        i += 1

    # dictionary for returning values
    dict_return = {
        "success" : True,
        "f_winch": F_winch,
        "rope_length" : rope_length,
        "lr_start" : lr_start,
        "a_exp" : a_exp,
        "b_exp" : b_exp,
        "t" : t,
        "phi_ot" : phi_ot,
        "omega_ot" : omega_ot,
        "omega_ot_theory" : sol.T[1],
        "time_to_finish" : time_to_finish,
        "j" : j,
        "alpha" : alpha,
        "torque_f_mot": torque_f_mot

    }
    return dict_return


# ==========================
    # ==========================
    # MAIN
    # ==========================
    # ==========================
if __name__ == "__main__":
    main()