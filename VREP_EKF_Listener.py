#!/usr/bin/env python

import csv
import tf as transf
import message_filters
import rospy
import sympy
from sympy import *
from itertools import izip
#from functionHandle import writerf
from euler2mat import RTmtx
from math import sqrt
import numpy as np
from threading import Thread
from matplotlib import pyplot as plt
from geometry_msgs.msg import PoseStamped
from nav_msgs.msg import Odometry
from ps_msgs.msg import ps_platform_motion_msg
from ps_msgs.msg import ps_platform_steering_report_msg
from ps_msgs.msg import ps_platform_wheel_speed_report_msg

from time import sleep
from mpl_toolkits.mplot3d import Axes3D

import scipy.linalg
import sys
import matplotlib.dates as mdates
import math
from scipy.stats import norm
from sympy import Symbol, symbols, Matrix, sin, cos, sqrt, atan, im, pprint
from sympy import init_printing
_FLOAT_EPS_4 = np.finfo(float).eps * 4.0
init_printing(use_latex=True)



#START OF Constant Params

dt = 0.130

m = 941.20 	# Mass of vehicle in kgg
h = 0.84 	# m
tf = 1.28	# m
tr = tf
#d = tf/2
a = 0.968	# m
b = 1.392	# m
g = 9.81

Iz = 1125.67

Ca_FL = 40000.0#200520.0    # cornering stiffness / lateral stiffness
Ca_FR = 40000.0
Ca_RL = 40000.0#76426.0
Ca_RR = 40000.0

Ck_FL = 150000.0  #30.0        # braking stiffness / longitudinal stiffness
Ck_FR = 150000.0  #30.0
Ck_RL = 150000.0  #24.0
Ck_RR = 150000.0 #24.0


wheelr = (0.4/2.0) + 0.8

#END OF Constant Params



def x_handle(x, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR):
    # print('printing x: ')
    # print(x)
    xdash = x
    xdash[0] = (x[0] + x[2] * dt + 0.5 * x[4] * dt ** 2)  # x coord
    xdash[1] = (x[1] + x[3] * dt + 0.5 * x[5] * dt ** 2)  # y coord
    xdash[2] = (x[2] + x[4] * dt)  # x vel
    xdash[3] = (x[3] + x[5] * dt)  # y vel
    xdash[4] = (((1 / m) * (
    Ck_FL * x[24] * cos(delta) - Ca_FL * x[15] * sin(delta) + Ck_FR * x[25] * cos(delta) - Ca_FR * x[16] * sin(
        delta) + Ck_RL * x[26] + Ck_RR * x[27])) * cos(x[6]) - ((1 / m) * (
    Ck_FL * x[24] * sin(delta) + Ca_FL * x[15] * cos(delta) + Ck_FR * x[25] * sin(delta) + Ca_FR * x[16] * cos(
        delta) + Ca_RL * x[17] + Ca_RR * x[18])) * sin(x[6]))  # x accel
    xdash[5] = (((1 / m) * (
    Ck_FL * x[24] * sin(delta) + Ca_FL * x[15] * cos(delta) + Ck_FR * x[25] * sin(delta) + Ca_FR * x[16] * cos(
        delta) + Ca_RL * x[17] + Ca_RR * x[18])) * cos(x[6]) + ((1 / m) * (
    Ck_FL * x[24] * cos(delta) - Ca_FL * x[15] * sin(delta) + Ck_FR * x[25] * cos(delta) - Ca_FR * x[16] * sin(
        delta) + Ck_RL * x[26] + Ck_RR * x[27])) * sin(x[6]))  # y accel
    xdashYaw = (x[6] + x[9] * dt + 0.5 * x[12] * dt ** 2)  # yaw

    if xdashYaw <= np.pi and xdashYaw >= -np.pi:
        xdash[6] = xdashYaw
    elif xdashYaw > np.pi:
        xdash[6] = (-(2 * np.pi - xdashYaw))
    elif xdashYaw < -np.pi:
        xdash[6] = (-(-2 * np.pi - xdashYaw))

    xdashPitch = x[7]  # pitch
    if xdashPitch <= np.pi and xdashPitch >= -np.pi:
        xdash[7] = (x[7])
    elif xdashPitch > np.pi:
        xdash[7] = -(2 * np.pi - x[7])
    elif xdashPitch < -np.pi:
        xdash[7] = -(-2 * np.pi - x[7])

    xdashRoll = x[8]  # roll

    if xdashRoll <= np.pi and xdashRoll >= -np.pi:
        xdash[8] = x[8]
    elif xdashRoll > np.pi:
        xdash[8] = -(2 * np.pi - x[8])
    elif xdashRoll < -np.pi:
        xdash[8] = -(-2 * np.pi - x[8])

    xdash[9] = x[9] + x[12] * dt  # yaw rate
    xdash[10] = x[10]  # pitch rate
    xdash[11] = x[11]  # roll rate
    xdash[12] = (1 / Iz) * (-tf * 0.5 * (Ck_FL * x[24] * cos(delta) - Ca_FL * x[15] * sin(delta)) + tf * 0.5 * (
    Ck_FR * x[25] * cos(delta) - Ca_FR * x[16] * sin(delta)) - tr * 0.5 * Ck_RL * x[26] + tr * 0.5 * Ck_RR * x[27] +
                            a * (Ca_FL * x[15] * cos(delta) + Ck_FL * x[24] * sin(delta)) + a * (
                            Ca_FR * x[16] * cos(delta) + Ck_FR * x[25] * sin(delta)) - b * Ca_RL * x[17] - b * Ca_RR *
                            x[18])  # yaww accel
    xdash[13] = x[13]  # pitch accel
    xdash[14] = x[14]  # roll accel
    xdash[15] = delta - atan((x[3] * cos(x[6]) - x[2] * sin(x[6]) + a * x[9]) / (
    x[3] * sin(x[6]) + x[2] * cos(x[6]) - tf * 0.5 * x[9]))  # alpha slip angle FL + ve
    xdash[16] = delta - atan((x[3] * cos(x[6]) - x[2] * sin(x[6]) + a * x[9]) / (
    x[3] * sin(x[6]) + x[2] * cos(x[6]) + tf * 0.5 * x[9]))  # alpha slip angle FR + ve
    xdash[17] = atan((-(x[3] * cos(x[6]) - x[2] * sin(x[6])) + b * x[9]) / (
    x[3] * sin(x[6]) + x[2] * cos(x[6]) - tr * 0.5 * x[9]))  # alpha slip angle RL -ve
    xdash[18] = atan((-(x[3] * cos(x[6]) - x[2] * sin(x[6])) + b * x[9]) / (
    x[3] * sin(x[6]) + x[2] * cos(x[6]) + tr * 0.5 * x[9]))  # alpha slip angle RR -ve
    xdash[19] = atan((x[3] * cos(x[6]) - x[2] * sin(x[6])) / (x[3] * sin(x[6]) + x[2] * cos(x[6])))  # Beta slip angle
    xdash[20] = sqrt((x[3] * cos(x[6]) - x[2] * sin(x[6]) + a * x[9]) ** 2 + (
    x[3] * sin(x[6]) + x[2] * cos(x[6]) - tf * 0.5 * x[9]) ** 2) * cos(x[15])  # Actual wheel velocity FL
    xdash[21] = sqrt((x[3] * cos(x[6]) - x[2] * sin(x[6]) + a * x[9]) ** 2 + (
    x[3] * sin(x[6]) + x[2] * cos(x[6]) + tf * 0.5 * x[9]) ** 2) * cos(x[16])  # Actual wheel velocity FR
    xdash[22] = sqrt((x[3] * cos(x[6]) - x[2] * sin(x[6]) - b * x[9]) ** 2 + (
    x[3] * sin(x[6]) + x[2] * cos(x[6]) - tr * 0.5 * x[9]) ** 2) * cos(x[17])  # Actual wheel velocity RL
    xdash[23] = sqrt((x[3] * cos(x[6]) - x[2] * sin(x[6]) - b * x[9]) ** 2 + (
    x[3] * sin(x[6]) + x[2] * cos(x[6]) + tr * 0.5 * x[9]) ** 2) * cos(x[18])  # Actual wheel velocity RR



    if x[20] == 0.0 and angular_whlSpeed_FL == 0.0:
        xdash[24] = 0.0
    elif x[20] >= (angular_whlSpeed_FL * wheelr):  # Slip ratio FL
        xdash[24] = (angular_whlSpeed_FL * wheelr / x[20]) - 1.0
    else:
        xdash[24] = 1.0 - (x[20] / (angular_whlSpeed_FL * wheelr))

    if x[21] == 0.0 and angular_whlSpeed_FR == 0.0:
        xdash[25] = 0.0
    elif x[21] >= angular_whlSpeed_FR * wheelr:  # Slip ratio FR
        xdash[25] = (angular_whlSpeed_FR * wheelr / x[21]) - 1.0
    else:
        xdash[25] = 1.0 - (x[21] / (angular_whlSpeed_FR * wheelr))

    if x[22] == 0.0 and angular_whlSpeed_RL == 0.0:
        xdash[26] = 0.0
    elif x[22] >= angular_whlSpeed_RL * wheelr:  # Slip ratio RL
        xdash[26] = (angular_whlSpeed_RL * wheelr / x[22]) - 1.0
    else:
        xdash[26] = 1.0 - (x[22] / (angular_whlSpeed_RL * wheelr))

    if x[23] == 0.0 and angular_whlSpeed_RR == 0.0:
        xdash[27] = 0.0
    elif x[23] >= angular_whlSpeed_RR * wheelr:  # Slip ratio RR
        xdash[27] = (angular_whlSpeed_RR * wheelr / x[23]) - 1.0
    else:
        xdash[27] = 1.0 - (x[23] / (angular_whlSpeed_RR * wheelr))
    xdash[28] = (0.5 * m * g - m * x[5] * h / tf) * b / (a + b) - 0.5 * m * x[4] * h / (a + b)  # Fz FL
    xdash[29] = (0.5 * m * g + m * x[5] * h / tf) * b / (a + b) - 0.5 * m * x[4] * h / (a + b)  # Fz FR
    xdash[30] = (0.5 * m * g - m * x[5] * h / tr) * a / (a + b) + 0.5 * m * x[4] * h / (a + b)  # Fz RL
    xdash[31] = (0.5 * m * g + m * x[5] * h / tr) * a / (a + b) + 0.5 * m * x[4] * h / (a + b)  # Fz RR

    return xdash



def y_handle(x):

    y = []
    y.append(x[0])
    y.append(x[1])
    y.append(x[6])
    y.append(x[2])
    y.append(x[3])
    y.append(x[4])
    y.append(x[5])

    return y

n_states = 32
n_meas_states = 7

# EKF
def x_Jacob(x, delta_Jc , angular_whlSpeed_FL_Jc, angular_whlSpeed_FR_Jc, angular_whlSpeed_RL_Jc, angular_whlSpeed_RR_Jc):
    delta = delta_Jc
    awsFL = angular_whlSpeed_FL_Jc
    awsFR = angular_whlSpeed_FR_Jc
    awsRL = angular_whlSpeed_RL_Jc
    awsRR = angular_whlSpeed_RR_Jc
    ret_x_jacob = np.matrix([[1.0, 0.0, dt, 0.0, 0.5*dt**2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 1.0, 0.0, dt, 0.0, 0.5*dt**2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 1.0, 0.0, dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 1.0, 0.0, dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -(-Ca_FL*x[15]*sin(delta) - Ca_FR*x[16]*sin(delta) + Ck_FL*x[24]*cos(delta) + Ck_FR*x[25]*cos(delta) + Ck_RL*x[26] + Ck_RR*x[27])*sin(x[6])/m - (Ca_FL*x[15]*cos(delta) + Ca_FR*x[16]*cos(delta) + Ca_RL*x[17] + Ca_RR*x[18] + Ck_FL*x[24]*sin(delta) + Ck_FR*x[25]*sin(delta))*cos(x[6])/m, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -Ca_FL*sin(x[6])*cos(delta)/m - Ca_FL*sin(delta)*cos(x[6])/m, -Ca_FR*sin(x[6])*cos(delta)/m - Ca_FR*sin(delta)*cos(x[6])/m, -Ca_RL*sin(x[6])/m, -Ca_RR*sin(x[6])/m, 0.0, 0.0, 0.0, 0.0, 0.0, -Ck_FL*sin(x[6])*sin(delta)/m + Ck_FL*cos(x[6])*cos(delta)/m, -Ck_FR*sin(x[6])*sin(delta)/m + Ck_FR*cos(x[6])*cos(delta)/m, Ck_RL*cos(x[6])/m, Ck_RR*cos(x[6])/m, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (-Ca_FL*x[15]*sin(delta) - Ca_FR*x[16]*sin(delta) + Ck_FL*x[24]*cos(delta) + Ck_FR*x[25]*cos(delta) + Ck_RL*x[26] + Ck_RR*x[27])*cos(x[6])/m - (Ca_FL*x[15]*cos(delta) + Ca_FR*x[16]*cos(delta) + Ca_RL*x[17] + Ca_RR*x[18] + Ck_FL*x[24]*sin(delta) + Ck_FR*x[25]*sin(delta))*sin(x[6])/m, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -Ca_FL*sin(x[6])*sin(delta)/m + Ca_FL*cos(x[6])*cos(delta)/m, -Ca_FR*sin(x[6])*sin(delta)/m + Ca_FR*cos(x[6])*cos(delta)/m, Ca_RL*cos(x[6])/m, Ca_RR*cos(x[6])/m, 0.0, 0.0, 0.0, 0.0, 0.0, Ck_FL*sin(x[6])*cos(delta)/m + Ck_FL*sin(delta)*cos(x[6])/m, Ck_FR*sin(x[6])*cos(delta)/m + Ck_FR*sin(delta)*cos(x[6])/m, Ck_RL*sin(x[6])/m, Ck_RR*sin(x[6])/m, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, dt, 0.0, 0.0, 0.5*dt**2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (Ca_FL*a*cos(delta) + 0.5*Ca_FL*tf*sin(delta))/Iz, (Ca_FR*a*cos(delta) - 0.5*Ca_FR*tf*sin(delta))/Iz, -Ca_RL*b/Iz, -Ca_RR*b/Iz, 0.0, 0.0, 0.0, 0.0, 0.0, (Ck_FL*a*sin(delta) - 0.5*Ck_FL*tf*cos(delta))/Iz, (Ck_FR*a*sin(delta) + 0.5*Ck_FR*tf*cos(delta))/Iz, -0.5*Ck_RL*tr/Iz, 0.5*Ck_RR*tr/Iz, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, -(-(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)*cos(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2 - sin(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf))/((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2 + 1), -(-(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)*sin(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2 + cos(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf))/((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2 + 1), 0.0, 0.0, -((x[2]*sin(x[6]) - x[3]*cos(x[6]))*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2 + (-x[2]*cos(x[6]) - x[3]*sin(x[6]))/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf))/((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2 + 1), 0.0, 0.0, -(a/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf) + 0.5*tf*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2)/((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2 + 1), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, -(-(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)*cos(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2 - sin(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf))/((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2 + 1), -(-(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)*sin(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2 + cos(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf))/((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2 + 1), 0.0, 0.0, -((x[2]*sin(x[6]) - x[3]*cos(x[6]))*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2 + (-x[2]*cos(x[6]) - x[3]*sin(x[6]))/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf))/((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2 + 1), 0.0, 0.0, -(a/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf) - 0.5*tf*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2)/((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2 + 1), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, (-(x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)*cos(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2 + sin(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr))/((x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2 + 1), (-(x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)*sin(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2 - cos(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr))/((x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2 + 1), 0.0, 0.0, ((x[2]*sin(x[6]) - x[3]*cos(x[6]))*(x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]))/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr))/((x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2 + 1), 0.0, 0.0, (b/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr) + 0.5*tr*(x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2)/((x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2 + 1), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, (-(x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)*cos(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2 + sin(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr))/((x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2 + 1), (-(x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)*sin(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2 - cos(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr))/((x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2 + 1), 0.0, 0.0, ((x[2]*sin(x[6]) - x[3]*cos(x[6]))*(x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]))/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr))/((x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2 + 1), 0.0, 0.0, (b/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr) - 0.5*tr*(x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2)/((x[2]*sin(x[6]) - x[3]*cos(x[6]) + x[9]*b)**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2 + 1), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, (-(-x[2]*sin(x[6]) + x[3]*cos(x[6]))*cos(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]))**2 - sin(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6])))/((-x[2]*sin(x[6]) + x[3]*cos(x[6]))**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]))**2 + 1), (-(-x[2]*sin(x[6]) + x[3]*cos(x[6]))*sin(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6]))**2 + cos(x[6])/(x[2]*cos(x[6]) + x[3]*sin(x[6])))/((-x[2]*sin(x[6]) + x[3]*cos(x[6]))**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]))**2 + 1), 0.0, 0.0, ((-x[2]*sin(x[6]) + x[3]*cos(x[6]))*(x[2]*sin(x[6]) - x[3]*cos(x[6]))/(x[2]*cos(x[6]) + x[3]*sin(x[6]))**2 + (-x[2]*cos(x[6]) - x[3]*sin(x[6]))/(x[2]*cos(x[6]) + x[3]*sin(x[6])))/((-x[2]*sin(x[6]) + x[3]*cos(x[6]))**2/(x[2]*cos(x[6]) + x[3]*sin(x[6]))**2 + 1), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, (-(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)*sin(x[6]) + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)*cos(x[6]))*cos(x[15])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2), ((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)*cos(x[6]) + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)*sin(x[6]))*cos(x[15])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2), 0.0, 0.0, ((-2*x[2]*sin(x[6]) + 2*x[3]*cos(x[6]))*(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)/2 + (-2*x[2]*cos(x[6]) - 2*x[3]*sin(x[6]))*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)/2)*cos(x[15])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2), 0.0, 0.0, (a*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a) - 0.5*tf*(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf))*cos(x[15])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2), 0.0, 0.0, 0.0, 0.0, 0.0, -sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tf)**2)*sin(x[15]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, (-(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)*sin(x[6]) + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)*cos(x[6]))*cos(x[16])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2), ((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)*cos(x[6]) + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)*sin(x[6]))*cos(x[16])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2), 0.0, 0.0, ((-2*x[2]*sin(x[6]) + 2*x[3]*cos(x[6]))*(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)/2 + (-2*x[2]*cos(x[6]) - 2*x[3]*sin(x[6]))*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)/2)*cos(x[16])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2), 0.0, 0.0, (a*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a) + 0.5*tf*(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf))*cos(x[16])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) + x[9]*a)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tf)**2)*sin(x[16]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, (-(-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)*sin(x[6]) + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)*cos(x[6]))*cos(x[17])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2), ((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)*cos(x[6]) + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)*sin(x[6]))*cos(x[17])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2), 0.0, 0.0, ((-2*x[2]*sin(x[6]) + 2*x[3]*cos(x[6]))*(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)/2 + (-2*x[2]*cos(x[6]) - 2*x[3]*sin(x[6]))*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)/2)*cos(x[17])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2), 0.0, 0.0, (-b*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b) - 0.5*tr*(x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr))*cos(x[17])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) - 0.5*x[9]*tr)**2)*sin(x[17]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, (-(-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)*sin(x[6]) + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)*cos(x[6]))*cos(x[18])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2), ((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)*cos(x[6]) + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)*sin(x[6]))*cos(x[18])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2), 0.0, 0.0, ((-2*x[2]*sin(x[6]) + 2*x[3]*cos(x[6]))*(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)/2 + (-2*x[2]*cos(x[6]) - 2*x[3]*sin(x[6]))*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)/2)*cos(x[18])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2), 0.0, 0.0, (-b*(-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b) + 0.5*tr*(x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr))*cos(x[18])/sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -sqrt((-x[2]*sin(x[6]) + x[3]*cos(x[6]) - x[9]*b)**2 + (x[2]*cos(x[6]) + x[3]*sin(x[6]) + 0.5*x[9]*tr)**2)*sin(x[18]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -awsFL*wheelr/x[20]**2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -awsFR*wheelr/x[21]**2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -awsRL*wheelr/x[22]**2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -awsRR*wheelr/x[23]**2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, -0.5*h*m/(a + b), -b*h*m/(tf*(a + b)), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, -0.5*h*m/(a + b), b*h*m/(tf*(a + b)), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.5*h*m/(a + b), -a*h*m/(tr*(a + b)), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.5*h*m/(a + b), a*h*m/(tr*(a + b)), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                             ])
    return ret_x_jacob

def y_Jacob():
    ret_y_jacob = np.matrix([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    return ret_y_jacob


def ekf(fstate,x,P,hmeas,z,Q,R, delta_ekf , angular_whlSpeed_FL_ekf, angular_whlSpeed_FR_ekf, angular_whlSpeed_RL_ekf, angular_whlSpeed_RR_ekf):

    #print('>>>>>>>>>>>>>>>>>>>>>>>>X before fstate :')
    #print(type(x))
    x1 = fstate(x, delta_ekf , angular_whlSpeed_FL_ekf, angular_whlSpeed_FR_ekf, angular_whlSpeed_RL_ekf, angular_whlSpeed_RR_ekf)
    A = x_Jacob(x1, delta_ekf , angular_whlSpeed_FL_ekf, angular_whlSpeed_FR_ekf, angular_whlSpeed_RL_ekf, angular_whlSpeed_RR_ekf)                     #nonlinear update and linearization at current state
    #print('X after fstate <<<<<<<<<<<<<<<<<<<<<<<<:')
    #print(type(x1))
    p = Q
    P = np.add((np.dot(A,P)).dot(A.T),Q)                               #partial update)
    z1 = hmeas(x1)
    H = y_Jacob()                    #nonlinear measurement and linearization
    P12=P.dot(np.transpose(H))                                   #cross covariance
    invHolder = sympy.Matrix((np.add(H.dot(P12),R)).tolist())
    invoHolder = invHolder.inv()
    K=P12.dot(invoHolder.tolist())                       #Kalman filter gain
    kdotHolder = K.dot([a_i - b_i for a_i, b_i in zip(z, z1)])
    KdotHolder = kdotHolder.tolist()
    x=np.add(x1,KdotHolder[0])                            #state estimate

    P=P-K.dot(np.transpose(P12))                               #state covariance matrix
    return x, P


# Code starts

q_proc = 0.1
r_meas = 0.1

Q_proc = (q_proc**2)*np.identity(n_states)

Q_proc_e, Q_proc_E = np.linalg.eig(Q_proc)  # value vector

R_meas = (r_meas**2)*np.identity(n_meas_states)

Q_proc[6,6] = 0.01**2
Q_proc[9,9] = 0.01**2
Q_proc[12,12] = 0.01**2


R_meas_e, R_meas_E = np.linalg.eig(R_meas)


s_state = np.zeros(n_states)
#s_state = np.asmatrix(s_state).T

x = s_state + (q_proc**2)*np.random.normal(0,0.01,n_states)

P = 1*np.identity(n_states)

# Preallocation for Storage
x0  = [] # x
x1  = [] # x vel
x2  = [] # z
x3  = [] # z vel
x4  = [] # yaw
x5  = [] # yaw rate
x6  = [] # roll
x7  = [] # roll rate
x8  = [] # pitch
x9  = [] # pitch rate
x10 = [] # y
x11 = [] # y vel
x12 = [] # Vxy
x13 = [] #ax
x14 = [] # az
x15 = [] # ay
x16 = [] # Axy
x17 = [] # fxFL
x18 = [] # fxFR
x19 = [] # fxRL
x20 = [] # fxRR
x21 = [] # fyFL
x22 = [] # fyFR
x23 = [] # fyRL
x24 = [] # fyRR
x25 = [] # fzFL
x26 = [] # fzFR
x27 = [] # fzRL


iterI = 0


# Transformation mtx




print("Processed through constants")

def listener():
    vo_pose = message_filters.Subscriber('/mono_odometer/pose', PoseStamped)
    whl_speed = message_filters.Subscriber('/wheel_speeds', ps_platform_wheel_speed_report_msg)
    motion = message_filters.Subscriber('/motion_data', ps_platform_motion_msg)
    steer = message_filters.Subscriber('/vehicle_report/steering', ps_platform_steering_report_msg)
    #sensOut4 = message_filters.Subscriber('/odom', Odometry)
    ts = message_filters.ApproximateTimeSynchronizer([vo_pose, whl_speed, motion, steer], 2, 0.1)
    ts.registerCallback(EKF)

    rospy.spin()



def EKF(vop, whlSpeed, motn, strAng):
    print('EKF got launched')
    rospy.loginfo(rospy.get_caller_id() + '  I start %s', iterI)
    global x, P, iterI
    global x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27
    voTmtx = RTmtx(-0.32753, -0.057194, -0.62872, 0, -89.99, -90)
    voTmtx_rotation = RTmtx(0, 0, 0, 0, -89.99, -90)
    senTmtx = RTmtx(+1.4059, -0.03949, +0.4788, 0, 0, 0)

    voCoordinates = np.dot(voTmtx, [[vop.pose.position.x],[vop.pose.position.y],[vop.pose.position.z], [1.0]])
    voOrientation = np.dot(voTmtx_rotation, [[vop.pose.orientation.x],[vop.pose.orientation.y],[vop.pose.orientation.z], [1.0]])
    voEulers = transf.transformations.euler_from_quaternion([voOrientation[0], voOrientation[1], voOrientation[2], vop.pose.orientation.w], 'szyx')

    senCoordinates = np.dot(senTmtx, [[motn.position[0]], [motn.position[1]], [motn.position[2]], [1.0]])
    senCoordinates = senCoordinates.tolist()
    senOrientation = [motn.orientation[0], motn.orientation[1], motn.orientation[2], 1.0] #np.dot(senTmtx, [[motn.orientation[0]], [motn.orientation[1]], [motn.orientation[2]], [1.0]])

    senEulers = transf.transformations.euler_from_quaternion([senOrientation[0], senOrientation[1], senOrientation[2], motn.orientation[3]], 'szyx')


    senRotationRate = [motn.rotation_rate[0], motn.rotation_rate[1], motn.rotation_rate[2], 1.0]  # np.dot(senTmtx, [[motn.rotation_rate[0]], [motn.rotation_rate[1]], [motn.rotation_rate[2]], [1.0]])
    senVelocity = [motn.velocity[0], motn.velocity[1], motn.velocity[2], 1.0]  # np.dot(senTmtx, [[motn.velocity[0]], [motn.velocity[1]], [motn.velocity[2]], [1.0]])
    senAcceleration = [motn.acceleration[0], motn.acceleration[1], motn.acceleration[2], 1.0]  #np.dot(senTmtx, [[motn.acceleration[0]], [motn.acceleration[1]], [motn.acceleration[2]], [1.0]])

    # Steering angle "strAng.steering_wheel_angle" in radians
    delta_cin = strAng.steering_wheel_angle
    # WheelSpeed "whlSpeed.front/rear_left/right" in radians/s
    angular_whlSpeed_FL_cin = whlSpeed.front_left
    angular_whlSpeed_FR_cin = whlSpeed.front_right
    angular_whlSpeed_RL_cin = whlSpeed.rear_left
    angular_whlSpeed_RR_cin = whlSpeed.rear_right

    meas = [senCoordinates[0][0], senCoordinates[1][0], senEulers[0], senVelocity[0], senVelocity[1], senAcceleration[0], senAcceleration[1]]
    #print('x : ')
    #print(senCoordinates[0])
    #print('y : ')
    #print(senCoordinates[1])
    #print('yaw : ')
    #print(senEulers[0])
    #print('x vel : ')
    #print(senVelocity[0])
    #print('y vel : ')
    #print(senVelocity[1])
    #print('x accel : ')
    #print(senAcceleration[0])
    #print('y accel : ')
    #print(senAcceleration[1])

    # fxfl = Ck_FL*x[20] fxfr = Ck_FR*x[21] fxrl = Ck_RL*x[22] fxrr = Ck_RR*x[23]
    # fyfl = Ca_FL*x[15] fyfr = Ca_FR*x[16] fyrl = Ca_RL*x[17] fyrr = Ca_RR*x[18]

    x, P = ekf(x_handle,x,P,y_handle,meas,Q_proc,R_meas, delta_cin, angular_whlSpeed_FL_cin, angular_whlSpeed_FR_cin, angular_whlSpeed_RL_cin, angular_whlSpeed_RR_cin)
    #print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Printing x at end<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
    #print(np.array(x))
    iterI += 1


    # Save states for Plotting
    '''
    with open('/home/equatrace/catkin_ws/src/viso2/viso2_ros/rata/Statdata.csv', 'wb') as fi:
        write = csv.writer(fi)
        write.writerows(izip(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27))# ,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49))

    #writerf(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,datafile)
    '''
    rospy.loginfo(rospy.get_caller_id() + '  I finished %s', iterI-1)
    return

if __name__ == "__main__":
    rospy.init_node("dissertacao")
    listener()

