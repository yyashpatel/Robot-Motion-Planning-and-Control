import numpy as np
import time
from matplotlib import pyplot as plt
import math

class BOT():
	def __init__(self,n, s, S, m, time, x_in,  y_in, theta_in, x_f, y_f, theta_f, T, k, i, s_dot, alpha_y, alpha_x, beta_y, beta_x, x_dot_s, y_dot_s, x_double_dot_s, y_double_dot_s, velocity):

		self.n=n
		self.m=m
		self.s=s
		self.S=S
		self.time=time
		self.T=T

		self.x_in = x_in
		self.y_in = y_in
		self.theta_in = theta_in

		self.x_f = x_f
		self.y_f = y_f
		self.theta_f = theta_f

		self.k=k
		self.i=i
		self.s_dot=s_dot

		self.alpha_x=alpha_x
		self.alpha_y=alpha_y
		self.beta_x=beta_x
		self.beta_y=beta_y

		self.x_traj=[]
		self.y_traj=[]
		self.theta_traj=[]

		self.v=[]
		self.w=[]

		self.x_path=[]
		self.y_path=[]
		self.theta_path=[]

		self.x_s=x_in
		self.y_s=y_in
		self.theta_s=theta_in

		self.x_bot=x_in
		self.y_bot=y_in
		self.theta_bot=theta_in

	#compute the geometric path of the robot
	def robot_path(self):
		#computes x,y coordinates that robot traverses
		for s in self.S:
			#equation x(s) computes x coordinate
			x = ((s**3)*self.x_f) - (((s-1)**3)*self.x_in) + (self.alpha_x*(s**2)*(s-1)) + ((self.beta_x)*(s)*((s-1)**2))
			#equation y(s) computes y coordinate
			y = ((s**3)*self.y_f) - (((s-1)**3)*self.y_in) + (self.alpha_y*(s**2)*(s-1)) + ((self.beta_y)*(s)*((s-1)**2))

			#derivative of x(s)
			x_dot_s = 3*(s**2)*self.x_f - 3*((s-1)**2)*self.x_in + self.alpha_x*(3*(s**2)-2*s) + self.beta_x*(3*(s**2)-4*s+1)
			#derivative of y(s)
			y_dot_s = 3*(s**2)*self.y_f - 3*((s-1)**2)*self.y_in + self.alpha_y*(3*(s**2)-2*s) + self.beta_y*(3*(s**2)-4*s+1)

			x_double_dot_s = 6*s*self.x_f - 6*(s-1)*self.x_in + self.alpha_x*(6*s-2) + self.beta_x*(6*s-4)
			y_double_dot_s = 6*s*self.y_f - 6*(s-1)*self.y_in + self.alpha_y*(6*s-2) + self.beta_y*(6*s-4)

			#calculated velocity that the robot requires
			self.velocity = math.sqrt(pow(x_dot_s,2) + pow(y_dot_s,2))
			self.v.append(self.velocity)
			
			#calculated angular velocity that robot requires
			omega = ((y_double_dot_s*x_dot_s) - (x_double_dot_s*y_dot_s)) / self.velocity**2
			self.w.append(omega)

			# unicycle robot kinematics, to update the robot position
			self.theta_s = math.atan2(y_dot_s,x_dot_s)
			self.theta_path.append(self.theta_s)

			self.x_s = self.x_s + self.velocity*math.cos(self.theta_s)*self.T
			self.x_path.append(self.x_s)

			self.y_s = self.y_s + self.velocity*math.sin(self.theta_s)*self.T
			self.y_path.append(self.y_s)
        
	#compute the trajectory of the robot with the timing law
	def robot_traj(self):

		for t in self.time:
			vel_t = self.v[self.i]*self.s_dot     #s_dot is the timing law, self.v and self.w are the velocities calculated for the geometric path
			omega_t = self.w[self.i]*self.s_dot
			#print(vel_t, omega_t)
			
			#use differential kinematics to compute x,y,theta
			self.theta_bot = self.theta_bot + omega_t*0.032
			self.theta_traj.append(self.theta_bot)

			self.x_bot = self.x_bot + vel_t*math.cos(self.theta_bot)*0.032
			self.x_traj.append(self.x_bot)

			self.y_bot = self.y_bot + vel_t*math.sin(self.theta_bot)*0.032 
			self.y_traj.append(self.y_bot)
			#print(t,i)
			self.i = self.i+1

	def plot(self):
		
		plt.plot(self.x_traj , self.y_traj)
		plt.xlabel("X axis")
		plt.ylabel("Y axis")
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()

def main():

	#S=[0,1], geometric path parameter
	n = 1
	S = np.linspace(0,n,300)
	print(S)

	#define timing law(s_dot = 1/t,here t=m)
	#here we take 9.6 because 9.6/0.032=300, where 300 is the totalnumber of steps
	m= 9.6
	s_dot = 1/m

	#total time taken to complete the trajectory
	time = np.linspace(0,m,300)

	#start position of the robot
	x_in = 0.0
	y_in = 0.0
	theta_in = 0

	#goal position of the robot
	x_f = 0.0
	y_f = 1.0
	theta_f = 0

	#time step
	T = 1/300

	#parameter to scale robot path
	k=10
	i=0
	
	#for x(s) cubic polynomial equation
	alpha_x = k*math.cos(theta_f) - 3*x_f
	beta_x = k*math.cos(theta_in) + 3*x_in

	#for y(s) cubic polynomial equation
	alpha_y = k*math.sin(theta_f) - 3*y_f
	beta_y = k*math.sin(theta_in) + 3*y_in

	#just defining dertivative of x(s) and y(s)
	s=0
	x_dot_s = 3*(s**2)*x_f - 3*((s-1)**2)*x_in + alpha_x*(3*(s**2)-2*s) + beta_x*(3*(s**2)-4*s+1)
	y_dot_s = 3*(s**2)*y_f - 3*((s-1)**2)*y_in + alpha_y*(3*(s**2)-2*s) + beta_y*(3*(s**2)-4*s+1)

	# second derivative of x(s) and y(s)
	x_double_dot_s = 6*s*x_f - 6*(s-1)*x_in + alpha_x*(6*s-2) + beta_x*(6*s-4)
	y_double_dot_s = 6*s*y_f - 6*(s-1)*y_in + alpha_y*(6*s-2) + beta_y*(6*s-4)

	#compute velocity according to given equation in the attached references
	velocity = math.sqrt(pow(x_dot_s,2) + pow(y_dot_s,2))

	bot = BOT(n, s, S, m, time, x_in,  y_in, theta_in, x_f, y_f, theta_f, T, k, i, s_dot, alpha_y, alpha_x, beta_y, beta_x, x_dot_s, y_dot_s, x_double_dot_s, y_double_dot_s, velocity)
	bot.robot_path()
	bot.robot_traj()
	bot.plot()

if __name__ == "__main__":
    main()
