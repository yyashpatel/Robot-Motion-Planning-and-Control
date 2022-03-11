import numpy as np
import time
from matplotlib import pyplot as plt
import math

t=20
time = np.linspace(0, t, 200,endpoint=True)
i=0

#time step
T=0.1

#starting point robot 1
x1 = 1
y1 = 1
x1_list = [] 
y1_list = []
x1_list.append(x1)
y1_list.append(y1)

#starting point robot 2
x2 = 10
y2 = 1
x2_list= []
y2_list = []
x2_list.append(x2)
y2_list.append(y2)

#starting point robot 3
x3 = 5
y3 = 8
x3_list= []
y3_list = []
x3_list.append(x3)
y3_list.append(y3)


while time[i] <= 19:

	#consensus equation
	x1_dot =  -((x1 - x2) + (x1 - x3))
	y1_dot = -((y1 - y2) + (y1 - y3))

	x2_dot = -((x2 - x1) + (x2 - x3))
	y2_dot = -((y2 - y1) + (y2 - x3))

	x3_dot = -((x3 - x1) + (x3 - x2))
	y3_dot = -((y3 - y1) + (y3 - x2))

	#robot position update using kinematics
	x1 = x1 + x1_dot*T
	x1_list.append(x1)
	y1 = y1 + y1_dot*T
	y1_list.append(y1)

	x2 = x2 + x2_dot*T
	x2_list.append(x2)
	y2 = y2 + y2_dot*T
	y2_list.append(y2)

	x3 = x3 + x3_dot*T
	x3_list.append(x3)
	y3 = y3 + y3_dot*T
	y3_list.append(y3)

	print(x1_dot,y1_dot)
	i=i+1
	#print(i)

print(x1,y1,x2,y2,x3,y3)
plt.plot(x1_list, y1_list)
plt.plot(x2_list, y2_list)
plt.plot(x3_list, y3_list)
plt.xlabel("X axis")
plt.ylabel("Y axis")
plt.show()