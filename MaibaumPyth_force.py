from typing import Tuple
import matplotlib as mlt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LightSource
from mpl_toolkits import mplot3d
import math

# Variables
l_St = 40/3
nan = np.nan
katnr = 8
l_centmass = 20
f_gravity = 20000-(300*katnr)
res = 500
res_Optimizing = 50
res_3D = 100
letter_phi = "φ"
letterC_phi = "Φ"

# creating values for the axis
lr_min = math.sqrt((math.cos(math.radians(60))*l_St)**2+(l_St)**2)
lr_max = math.sqrt((100-math.cos(math.radians(60))*l_St)**2+(math.sin(math.radians(60))*l_St)**2)+(80/3)

def main():

    fig_1_GridRows = 1
    fig_1_GridCols = 1
    fig_2_GridRows = 2
    fig_2_GridCols = 2
    

    
    # =============
        # =============
        # DATA
        # =============
        # =============

    #### 2D Plot (higher resolution)
    a = np.linspace(l_St, 100, res)
    lr = np.linspace(lr_min, lr_max, res)
    # # initializing force
    # F_winch = np.array([[nan for i in range(res) for i in range(res)]])

    #### Optimizing (less resolution)
    a_Opt = np.linspace(l_St, 100, res_Optimizing)
    lr_Opt = np.linspace(lr_min, lr_max, res_Optimizing)

    #### 3D Plot (less resolution)
    a_3D = np.linspace(l_St, 100, res_3D)
    lr_3D = np.linspace(lr_min, lr_max, res_3D)
    # creating mesh grid for the axis
    A_3D, LR_3D = np.meshgrid(a_3D, lr_3D)
    # initializing force
    F_winch_3D = np.array([[nan for i in range(res_3D) for i in range(res_3D)]])



    # =============
        # =============
        # calculations
        # =============
        # =============


    # calc max-force optimized a
    # min = np.inf
    # a_Fmin = 0
    # for var_a in a_Opt:
    #     max = 0
    #     for var_lr in lr_Opt:
    #         x = func_Force_over_LengthRope(var_a, var_lr)
    #         if x > max:
    #             max = x
    #     if min > max:
    #         min = max
    #         a_Fmin = var_a

    # print(f"The distanc a at which the maximum force ({min/1000}kN) is the lowest is: {a_Fmin}m")

    # =============
        # =============
        # 1 FIGURE
        # =============
        # =============
    # create a figure
    fig_1 = plt.figure(figsize=(18, 9), constrained_layout=True)

    # =============
        # 1 subplot
        # =============
    # add first axes: axs_1
    ax_1 = fig_1.add_subplot(fig_1_GridRows, fig_1_GridCols, 1, projection="3d")
    # title and label of: axs_1
    ax_1.set_xlabel("a in [m]")
    ax_1.set_ylabel("length of rope in [m]")
    ax_1.set_zlabel("force at winch in [kN]")
    ax_1.set_title(f"force at winch, Ktlg.Nr.: {katnr}")
    # getting Z - Data for graph
    F_winch_3D = func_Force_over_LengthRope(A_3D, LR_3D)/1000
    # creating colors for colormap
    ax_1_colors = np.linspace(1, 100, 7718)
    # marker ticker setting of axes
    loc = mlt.ticker.MultipleLocator(base=10.0)
    ax_1.xaxis.set_major_locator(loc)
    ax_1.yaxis.set_major_locator(loc)
    loc = mlt.ticker.MultipleLocator(base=5.0)
    ax_1.zaxis.set_major_locator(loc)
    # size setting of axes
    ax_1.set_zlim([0,34])
    # creating a color map
    ax_1_cmap = plt.get_cmap("hsv")
    # add gridlines
    ax_1.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.2)
    # Change viewing angle to reorient camera x degrees and rotate y
    ax_1.view_init(25,30)
    # Ploting the 3d Graph
    #sctt = ax_1.scatter3D(A,LR,F_winch, alpha = 0.8, c = ax_1_colors, cmap = ax_1_cmap)#s=1.5, cmap = ax_1_cmap)
    #fig_1.colorbar(sctt, ax = ax_1, shrink = 0.5, aspect = 5)
    ax_1.scatter3D(A_3D,LR_3D,F_winch_3D, c = "#EA8000", s=1, marker="o")
    #ax_1.plot_surface(A_3D,LR_3D,F_winch_3D, cmap=ax_1_cmap)



    # =============
        # =============
        # 2 FIGURE
        # =============
        # =============
    # create a figure
    fig_2 = plt.figure(figsize=(18, 9), constrained_layout=True)

# =============
        # 2 FIGURE
        # 1 subplot
        # =============
    # add axes ax_2
    ax_1 = fig_2.add_subplot(fig_2_GridRows, fig_2_GridCols, 1)
    # axes labeling and title
    ax_1.set_xlabel("length of rope in [m]")
    ax_1.set_ylabel("force at winch in [kN]")
    # ax_1.set_title(f"force at winch (a = {round(a_Fmin, 3)}m)")
    ax_1.set_title(f"force at winch (a = {23.233}m)")
    #### Set the X - marker stepsize
    loc = mlt.ticker.MultipleLocator(base=2.0)
    ax_1.xaxis.set_major_locator(loc)
    #### Set the Y - marker stepsize
    loc = mlt.ticker.MultipleLocator(base=2.0)
    ax_1.yaxis.set_major_locator(loc)
    # setting grid options
    ax_1.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
    # size setting of axes
    #ax_2.set_xlim([5,21])
    #ax_2.set_ylim([80,140])
    # plotting axis 2
    ax_1.plot(lr, func_Force_over_LengthRope(23.233, lr)/1000)

    # =============
        # 2 FIGURE
        # 2 subplot
        # =============
    # add axes ax_2
    ax_2 = fig_2.add_subplot(fig_2_GridRows, fig_2_GridCols, 2)
    # axes labeling and title
    ax_2.set_xlabel("length of rope in [m]")
    ax_2.set_ylabel("force at winch in [kN]")
    ax_2.set_title("force at winch (a = 40/3m)")
    #### Set the X - marker stepsize
    loc = mlt.ticker.MultipleLocator(base=2.0)
    ax_2.xaxis.set_major_locator(loc)
    #### Set the Y - marker stepsize
    loc = mlt.ticker.MultipleLocator(base=2.0)
    ax_2.yaxis.set_major_locator(loc)
    # setting grid options
    ax_2.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
    # size setting of axes
    #ax_2.set_xlim([5,21])
    #ax_2.set_ylim([80,140])
    # plotting axis 2
    ax_2.plot(lr, func_Force_over_LengthRope(l_St, lr)/1000)

    # =============
        # 2 FIGURE
        # 3 subplot
        # =============
    # add axes ax_3
    ax_3 = fig_2.add_subplot(fig_2_GridRows, fig_2_GridCols, 3)
    # axes labeling and title
    ax_3.set_xlabel("length of rope in [m]")
    ax_3.set_ylabel("force at winch in [kN]")
    ax_3.set_title("force at winch (a = 14m)")
    #### Set the X - marker stepsize
    loc = mlt.ticker.MultipleLocator(base=2.0)
    ax_3.xaxis.set_major_locator(loc)
    #### Set the Y - marker stepsize
    loc = mlt.ticker.MultipleLocator(base=2.0)
    ax_3.yaxis.set_major_locator(loc)
    # setting grid options
    ax_3.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
    # size setting of axes
    #ax_3.set_xlim([5,21])
    #ax_3.set_ylim([80,140])
    # plotting axis 3
    ax_3.plot(lr, func_Force_over_LengthRope(14, lr)/1000)

    # =============
        # 2 FIGURE
        # 4 subplot
        # =============
    # add axes ax_4
    ax_4 = fig_2.add_subplot(fig_2_GridRows, fig_2_GridCols, 4)
    # axes labeling and title
    ax_4.set_xlabel("length of rope in [m]")
    ax_4.set_ylabel("force at winch in [kN]")
    ax_4.set_title("force at winch (a = 100m)")
    #### Set the X - marker stepsize
    loc = mlt.ticker.MultipleLocator(base=2.0)
    ax_4.xaxis.set_major_locator(loc)
    #### Set the Y - marker stepsize
    loc = mlt.ticker.MultipleLocator(base=2.0)
    ax_4.yaxis.set_major_locator(loc)
    # setting grid options
    ax_4.grid(visible = True, color ='grey', linestyle ='-.', linewidth = 0.3, alpha = 0.8)
    # size setting of axes
    #ax_4.set_xlim([5,21])
    #ax_4.set_ylim([80,140])
    # plotting axis 4
    ax_4.plot(lr, func_Force_over_LengthRope(100, lr)/1000)


    # =============
        # =============
        # display Plots
        # =============
        # =============
    plt.show()


# ==========================
    # ==========================
    # functions
    # ==========================
    # ==========================
@np.vectorize
def func_phi(fu_z, fu_a):
    phi_i = ((fu_z**2)+(l_St**2)-(fu_a**2))/(2*fu_z*l_St)
    if (phi_i > 1) or (phi_i < -1):
        return nan
    elif (np.degrees(np.arccos(phi_i))+60) > 180:
        return nan
    else:
        return np.degrees(np.arccos(phi_i))

@np.vectorize
def func_gamma(fu_a, fu_z):
    gamma_i = ((fu_a**2)+(l_St**2)-(fu_z**2))/(2*fu_a*l_St)
    if (gamma_i > 1) or (gamma_i < -1):
        return nan
    else:
        return np.degrees(np.arccos(gamma_i))

@np.vectorize
def func_Force_over_LengthRope(fu_a, fu_lr):
    phi_1 = 0.0
    phi_2 = 0.0
    phi_3 = 0.0
    phi_3_max = 0.0
    gamma = 0.0
    Z1 = False
    Z2 = False
    Z3 = False
    f_rope_tree = 0.0
    f_rope_winch = nan
    z1_max = np.sqrt((fu_a-np.cos(np.radians(60))*l_St)**2+(np.sin(np.radians(60))*l_St)**2)
    z1_min = 0
    z1 = fu_lr - (80/3)
    z2 = fu_lr - (40/3)
    z3 = fu_lr

    #### max length of rope for Z1
    for i in np.linspace(0, 130, 250):
        z = 130 - i
        if z == 0:
            a = nan
        else:
            a = (l_St**2 + z**2 - fu_a**2) / (2*l_St*z)
        if (a >= -1) and (a <= 1):
            if np.degrees(np.arccos(a)) >= 120:
                z1_min = z
                break

    #### max length of rope
    if fu_lr > (z1_max + 80/3):
        return nan

    #### min length of rope, else more than 90 degrees of the tree
    if fu_lr < (np.sqrt((fu_a**2)+(l_St**2))):
        return nan

    #### calc phi3 max 
    if (fu_a/z3 > 1) or (fu_a/z3 < -1):
        phi_3_max = nan
    else: phi_3_max = np.degrees(np.arcsin(fu_a/z3))

    #### calc phi_1 start
    phi_1_start = func_phi(z1_max, fu_a)

    #### calc phi
    if (fu_lr - (80/3)) <= z1_max:
        phi_1 = func_phi(z1, fu_a)

    if (phi_1 <= 120) and not(np.isnan(phi_1)) and (z1 > z1_min): # and (phi_1 > phi_1_start)
        Z1 = True
    else:
        phi_2 = func_phi(z2, fu_a)
        if (phi_2 <= 120) and not(np.isnan(phi_2)):
            Z2 = True
        else:
            phi_3 = func_phi(z3, fu_a)
            if (phi_3 <= phi_3_max) and not(np.isnan(phi_3)):
                Z3 = True

    #### calc gamma and force
    if Z1 and not(Z2 or Z3):
        gamma = func_gamma(fu_a, z1)
        #### calc rope force at tree
        f_rope_tree = (f_gravity*np.cos(np.radians(60-gamma))*l_centmass) / (np.sin(np.radians(60))*l_St)
        #### calc rope force at winch
        if (2*np.sin(np.radians(phi_1))) == 0:
            f_rope_winch = 0
        else:
            f_rope_winch = (f_rope_tree * np.sqrt(3)) / (2*np.sin(np.radians(phi_1)))
    elif Z2 and not(Z1 or Z3):
        gamma = func_gamma(fu_a, z2)
        if gamma < 30:      # tree already stands
            return nan
        #### calc rope force at tree
        f_rope_tree = (f_gravity*np.cos(np.radians(120-gamma))*l_centmass) / (np.sin(np.radians(60))*l_St)
        #### calc rope force at winch
        f_rope_winch = (f_rope_tree * np.sqrt(3)) / (2*np.sin(np.radians(phi_2)))
    elif Z3 and not(Z2 or Z1):
        gamma = func_gamma(fu_a, z3)
        if gamma < 90:
            return nan
        #### calc rope force at winch
        f_rope_winch = (f_gravity*np.cos(np.radians(180-gamma))*l_centmass) / (np.sin(np.radians(phi_3))*l_St)

    ## DEBUGING START
    if fu_lr < 24.2 and fu_lr > 20:
        i = 0
    ## DEBUGING END

    #### return f_rope_winch
    if f_rope_winch < 0:# or np.inf:
        return nan
    else: return f_rope_winch


# ==========================
    # ==========================
    # MAIN
    # ==========================
    # ==========================
if __name__ == "__main__":
    main()
