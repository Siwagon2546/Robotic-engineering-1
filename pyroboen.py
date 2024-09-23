import numpy as np 
import math as m 


pi = np.pi

def d2r(val):
    return val*pi/180
    
def r2d(val):
    return val*180/pi
  

def Rotx(val):
    out = np.zeros((4, 4))
    out[0,:] = [1, 0, 0, 0]
    out[1,:] = [0, np.cos(val), -np.sin(val), 0]
    out[2,:] = [0, np.sin(val), np.cos(val), 0]
    out[3,:] = [0, 0, 0, 1]
    return out    
    
def Rotz(val):
    out = np.zeros((4, 4))
    out[0,:] = [np.cos(val), -np.sin(val), 0, 0]
    out[1,:] = [np.sin(val), np.cos(val), 0, 0]
    out[2,:] = [0, 0, 1, 0]
    out[3,:] = [0, 0, 0, 1]
    return out  
def Tranx(val):
    out = np.zeros((4, 4))
    out[0,:] = [1, 0, 0, val]
    out[1,:] = [0, 1, 0, 0]
    out[2,:] = [0, 0, 1, 0]
    out[3,:] = [0, 0, 0, 1]
    return out

def Tranz(val):
    out = np.zeros((4, 4))
    out[0,:] = [1, 0, 0, 0]
    out[1,:] = [0, 1, 0, 0]
    out[2,:] = [0, 0, 1, val]
    out[3,:] = [0, 0, 0, 1]
    return out

def DH_table(theta):
    DH = np.zeros((6, 4))
    DH[0,:] = [ d2r(0),    0,      0.0892,      d2r(theta[0])]
    DH[1,:] = [ d2r(90),   0,      0.110,       d2r(90+theta[1])]
    DH[2,:] = [ d2r(0),    0.425,     0,        d2r(theta[2])]
    DH[3,:] = [ d2r(0),    0.392,   0,          d2r(-90+theta[3])]
    DH[4,:] = [ d2r(-90),  0,      0.09475,     d2r(theta[4])]
    DH[5,:] = [ d2r(90),   0,      0.0825,      d2r(theta[5])]
    return DH

def T_pos_euler(matrix):
    position = matrix[:3, 3]
    r = matrix[:3, :3]
    psi = np.arctan2(-r[1][2],r[2][2])
    phi = np.arcsin(r[0][2])
    theta = np.arctan2(-r[0][1],r[0][0])
    euler = np.array([(psi), (phi), (theta)])
    return position, euler

def Homogeneous(DH,F,L):
    H = np.eye(4)
    for i in range(F,L):
        Rx = Rotx(DH[i][0])
        Tx = Tranx(DH[i][1])
        Tz = Tranz(DH[i][2])
        Rz = Rotz(DH[i][3])
        H = H @ Rx @ Tx  @ Tz @ Rz
    return H
def find_J(DH):
    r_6E_6 = np.append(T_pos_euler(Homogeneous(DH, 0, 6))[0],1)
    J = np.zeros((6, len(DH)))
    for i in range(len(DH)):
        dum = Homogeneous(DH, i+1,6) @ r_6E_6
        r_iE_6 = dum[:-1]
        T_i_0 = Homogeneous(DH, 0, i+1)
        R_i_0 = T_i_0[:3, :3]
        r_iE_0 = R_i_0 @ r_iE_6
        k_i_0 = R_i_0[:, 2]
        J[0:3, i] = np.cross(k_i_0,r_iE_0)
        J[3:6, i] = k_i_0
    return J

def FK():
    print(np.dot(Homogeneous(DH,6),[0,0,0,1]))

    
th = [0,0,0,0,0,0]
DH = DH_table(th)
print(np.dot(Homogeneous(DH,0,6),[0,0,0,1]))
print(T_pos_euler(Homogeneous(DH,0,6))[1])
print(find_J(DH))