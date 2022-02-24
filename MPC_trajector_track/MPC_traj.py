
import numpy as np
import time
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
import math
from qpsolvers import solve_qp

t=960
time = np.linspace(0, t, 30000,endpoint=True)
i=0

x_state_current=np.zeros((3,30000))
u_current=np.zeros((2,30000))

x_ref=np.zeros((3,30000))
u_ref=np.zeros((2,30000))

#robot  trajectory initial coordinate and velocity
x_state_current[:,0]=np.array([0,3,0])
u_current[:,0]=np.array([0.0,0.0])

#reference trajectory initial coordinate and velocity
x_ref[:,0]=np.array([0.0,0.0,0])
u_ref[:,0]=np.array([1.0,0.0])

radius_of_wheel = 0.033
dist_between_wheels = 0.16

x_cord = []
y_cord = []
x_cord.append(0.0)
y_cord.append(3.0)

x_ref_traj = []
y_ref_traj = []
x_ref_traj.append(0.0)
y_ref_traj.append(0.0)



def mpc_controller(x_state_current, u_current, x_ref,u_ref ,i):

    global T
    global N

    A = np.array([[1, 0, -u_ref[0,i]*math.sin(x_ref[2,i])*T],
                  [0, 1, u_ref[0,i]*math.cos(x_ref[2,i])*T],
                  [0, 0, 1]])

    B = np.array([[math.cos(x_ref[2,i])*T, 0], 
                  [math.sin(x_ref[2,i])*T, 0],
                  [0, T]])

    Q = 0.5*np.eye(np.size(A,0)*N)
    R = 0.5*np.eye(np.size(B,1)*N)

    x_traj = np.zeros((3, N))
    u_traj = np.zeros((2, N))

    x_traj[:,0] = x_ref[:,i]
    u_traj[:,0] = u_ref[:,i]

    for j in range(1, N):
        delta_x = np.array([[u_traj[0,j-1]*math.cos(x_traj[2,j-1])*T],
                            [u_traj[0,j-1]*math.sin(x_traj[2,j-1])*T],
                            [u_traj[1,j-1]*T]])

        x_traj[:,j] = x_traj[:,j-1] + delta_x[:,0]
        u_traj[:,j] = u_traj[:,j-1]

    A_dash = np.zeros((np.size(A,0)*N, np.size(A,1)))
    

    for j in range(N):
        A_i = np.array([[1, 0, -u_traj[0,j]*math.sin(x_traj[2,j])*T],
            [0, 1, u_traj[0,j]*math.cos(x_traj[2,j])*T],
            [0, 0, 1]])

        if(j == 0):
            ind = j*np.size(A,0)   
            A_dash[ind:ind+np.size(A,0), :] = A_i
        else:
            ind = j*np.size(A,0)   
            A_dash[ind:ind+np.size(A,0), :] = np.dot(A_dash[ind-np.size(A,0): ind , :], A_i)

    B_dash = np.zeros((np.size(B,0)*N, np.size(B,1)*N))
    for j in range(N):
        ind_i = j*np.size(B,0) 
        for l in range(N):
            ind_j = l*np.size(B,1) 
            alph = np.zeros((np.size(A,0), np.size(A,1)))
            if(l == j):
                alph = np.array([[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]])

                B_i = np.array([[math.cos(x_traj[2,j])*T, 0], 
                    [math.sin(x_traj[2,j])*T, 0],
                    [0, T]])
                
                B_dash[ind_i:ind_i+np.size(B,0), ind_j:ind_j+np.size(B,1)] = np.dot(alph, B_i)
                break
            else:
                alph = np.array([[1, 0, -u_traj[0,l+1]*math.sin(x_traj[2,l+1])*T],
                    [0, 1, u_traj[0,l+1]*math.cos(x_traj[2,l+1])*T],
                    [0, 0, 1]])

                for z in range(l+2,j):
                    alph = alph * np.array([[1, 0, -u_traj[0,z]*math.sin(x_traj[2,z])*T],
                                            [0, 1, u_traj[0,z]*math.cos(x_traj[2,z])*T],
                                            [0, 0, 1]])
                
                B_i = np.array([[math.cos(x_traj[2,l])*T, 0],
                                [math.sin(x_traj[2,l])*T, 0],
                                [0, T]])
                
                B_dash[ind_i:ind_i+np.size(B,0), ind_j:ind_j+np.size(B,1)] = np.dot(alph, B_i)
    #print(B_dash)
    
    x_tilda = x_state_current[:,i] - x_ref[:,i]
    H = 2 * (np.dot(np.dot(np.transpose(B_dash), Q), B_dash) + R)
    f = 2 * np.dot(np.dot(np.dot(np.transpose(B_dash), Q), A_dash), x_tilda)

    f = np.reshape(f, (2*N))
    #objective_func = 1/2*(np.transpose(u_current[:,i]-u_ref[:,i]))*H*(u_current[:,i]-u_ref[:,i]) + (np.transpose(f)*(u_current[:,i]-u_ref[:,i]))

    D = np.zeros((4*N,2*N))
    #print(D.shape)
    for j in range(N):
        ind = j*np.size(u_ref,0)   
        D[2*ind:2*ind+2*np.size(u_ref,0), ind:ind+np.size(u_ref,0)] = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])
    
    C = np.zeros((4*N, 1))
    
    for j in range(N):
        ind = 2*j*np.size(u_ref,0)
        C[ind:ind+2*np.size(u_ref,0),:]=np.array([[3-u_ref[0,i]],[0.3-u_ref[1,i]],[-(-3-u_ref[0,i])],[-(-0.3-u_ref[1,i])]])
    
    u_dash = solve_qp(H, f, D, C, solver="cvxopt")
    #u_dash = solve_qp(H, f, D, C)

    u_current[:,i]= u_dash[0:2] + u_ref[:,i]
    velocity_rw = (2*u_current[0][i] + dist_between_wheels*u_current[1][i]) / (2*radius_of_wheel)
    velocity_lw = ((2*u_current[0][i])/radius_of_wheel) - velocity_rw
    #print(u_current[:,i] , velocity_rw , velocity_lw)

    delta_x = np.array([[u_current[0,i]*math.cos(x_state_current[2,i])*T],
                        [u_current[0,i]*math.sin(x_state_current[2,i])*T],
                        [u_current[1,i]*T]])

    x_state_current[0,i+1] = x_state_current[0,i] + delta_x[0]
    x_cord.append(x_state_current[0,i+1])
    
    x_state_current[1,i+1] = x_state_current[1,i] + delta_x[1]
    y_cord.append(x_state_current[1,i+1])

    x_state_current[2,i+1] = x_state_current[2,i] + delta_x[2]
    #print(x_state_current[0,i+1], x_state_current[1,i+1], x_state_current[2,i+1])
    #print(delta_x)
    delta_x_ref = np.array([[u_ref[0,i]*math.cos(x_ref[2,i])*T],
                            [u_ref[0,i]*math.sin(x_ref[2,i])*T],
                            [u_ref[1,i]*T]])

    x_ref[0,i+1] = x_ref[0,i] + delta_x_ref[0]
    x_ref[1,i+1] = x_ref[1,i] + delta_x_ref[1]
    x_ref[2,i+1] = x_ref[2,i] + delta_x_ref[2]

    x_ref_traj.append(x_ref[0,i+1])
    y_ref_traj.append(x_ref[0,i+1])

    u_ref[:,i+1] = u_ref[:,i]


while time[i] < 960:

    T=0.032
    N=5
    
    mpc_controller(x_state_current, u_current, x_ref,u_ref ,i)

    i=i+1
    print(i)

print(x_cord)
print(y_cord)

#plot robot trajectory and reference trajectory
plt.plot(x_state_current[0],x_state_current[1],label="Robot Trajectory")
plt.plot(x_ref[0],x_ref[1],label="Reference Trajectory")
plt.xlabel("X axis")
plt.ylabel("Y axis")
plt.legend(loc = "upper left")
plt.show()