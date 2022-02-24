import math
import matplotlib.pyplot as plt
import numpy as np

class Robot():
    def __init__(self, x_ref, y_ref, theta_ref, velocity_ref, robot_x, robot_y, robot_theta, robot_x_goal, robot_y_goal, kP = 2, kI = 0.04, kD = 0.2):
        self.x_ref = x_ref
        self.y_ref = y_ref
        self.theta_ref = theta_ref
        self.x_ref_list = []
        self.y_ref_list = []
        self.velocity_ref = velocity_ref
        self.t_ref = 0.2
        self.omega_ref = 0   # for straight line

        self.robot_x = robot_x
        self.robot_y = robot_y
        self.robot_theta = robot_theta
        self.robot_vel = velocity_ref
        self.robot_omega = None
        self.robot_x_goal = robot_x_goal
        self.robot_y_goal = robot_y_goal

        self.robot_x_list = []
        self.robot_y_list = []

        self.Kp = kP
        self.Ki = kI
        self.Kd = kD

        self.cumm_error = 0
        self.prev_error = 0

    def makeTrajectory(self):
        self.x_ref += self.velocity_ref * math.cos(self.theta_ref) * self.t_ref
        self.y_ref += self.velocity_ref * math.sin(self.theta_ref) * self.t_ref
        self.theta_ref +=  self.omega_ref * self.t_ref
        self.x_ref_list.append(self.x_ref)
        self.y_ref_list.append(self.y_ref)
        return

    def trajectory(self, x_goal, y_goal, theta_goal):
        while True:
            self.makeTrajectory()
            if self.x_ref >= x_goal:
                break
        return

    def isArrived(self):
        if abs(self.robot_x - self.robot_x_goal) < 0.05 or abs(self.robot_y - self.robot_y_goal) < 0.05:
            return True
        return False

    def makeAction(self, w):
        theta_dt = w

        self.robot_x += self.robot_vel * math.cos(self.robot_theta) * self.t_ref
        self.robot_y += self.robot_vel * math.sin(self.robot_theta) * self.t_ref
        self.robot_theta += theta_dt*self.t_ref
        return

    def iteratePID(self):
        g_theta = math.atan2(self.robot_y_goal - self.robot_y, self.robot_x_goal - self.robot_x)
        e = g_theta - self.robot_theta

        e_P = e
        e_I = self.cumm_error + e
        e_D = e - self.prev_error

        w = self.Kp*e_P + self.Ki*e_I + self.Kd*e_D

        self.cumm_error += e
        self.prev_error = e

        return w

    def runPID(self):
        #print("dddddddddd")
        while(self.isArrived() == False):
            w = self.iteratePID()
            self.makeAction(w)
            self.robot_x_list.append(self.robot_x)
            self.robot_y_list.append(self.robot_y)
        return

    def plotReference(self): 
        plt.plot(self.x_ref_list, self.y_ref_list, label="Reference Trajectory")
        plt.plot(self.robot_x_list, self.robot_y_list, label="Robot Trajectory")
        plt.legend(loc = "upper left")
        plt.show()

def main():
    ref_start_x = 2
    ref_start_y = 2
    ref_lin_velocity = 1
    
    ref_goal_x = 8
    ref_goal_y = 4
    ref_goal_theta = math.atan2(ref_goal_y - ref_start_y , ref_goal_x - ref_start_x)

    robot_x = ref_start_x
    robot_y = ref_start_y
    robot_theta = 0

    robot_x_goal = ref_goal_x
    robot_y_goal = ref_goal_y
    
    robot = Robot(ref_start_x, ref_start_y, ref_goal_theta, ref_lin_velocity, robot_x, robot_y, robot_theta, robot_x_goal, robot_y_goal) 
    robot.trajectory(ref_goal_x, ref_goal_y, ref_goal_theta)
    robot.runPID()
    #print(robot.robot_x_list)
    robot.plotReference()

if __name__ == "__main__":
    main()