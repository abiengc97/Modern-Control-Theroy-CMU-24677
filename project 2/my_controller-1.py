# Fill in the respective functions to implement the controller

# Import libraries
import numpy as np
from base_controller import BaseController
from scipy import signal, linalg
from util import *

# CustomController class (inherits from BaseController)
class CustomController(BaseController):

    def __init__(self, trajectory):

        super().__init__(trajectory)

        # Define constants
        # These can be ignored in P1
        self.lr = 1.39
        self.lf = 1.55
        self.Ca = 20000
        self.Iz = 25854
        self.m = 1888.6
        self.g = 9.81

        # Add additional member variables according to your need here.
        self.previous_error=0
        self.integral_error=0
        self.previous_e1=0
        self.previous_e2=0
        
    def PID(self,current_error,kp,ki,kd):
        delt_time=0.01
        if self.previous_error is None:
            self.previous_error=current_error
        
        self.integral_error+=current_error*delt_time
        derivative_error=(current_error-self.previous_error)/delt_time
        self.previous_error=current_error
        return kp*current_error+ki*self.integral_error+kd*derivative_error
    

    def update(self, timestep):
        delt_time=0.01
        trajectory = self.trajectory

        lr = self.lr
        lf = self.lf
        Ca = self.Ca
        Iz = self.Iz
        m = self.m
        g = self.g

        # Fetch the states from the BaseController method
        delT, X, Y, xdot, ydot, psi, psidot = super().getStates(timestep)

        # Design your controllers in the spaces below. 
        A = np.array([[0, 1, 0, 0], [0, -4*Ca / (m * xdot), 4*Ca / m, -(2*Ca*(lf - lr))/(m*xdot)], [0, 0, 0, 1], [0, -(2*Ca*(lf - lr)) / (Iz * xdot), (2*Ca*(lf - lr)) / Iz, (-2*Ca*(np.power(lf, 2) + np.power(lr, 2))) / (Iz * xdot)]])
        B = np.array([[0], [2*Ca / m], [0], [(2 * Ca* lf) / Iz]])
        
        
        # ---------------|Lateral Controller|-------------------------
        # Remember, your controllers will need to use the states
        
        time_horizon=10
        _,nearidx=closestNode(X,Y,trajectory)
        [X_d,Y_d]=trajectory[nearidx+time_horizon]  #calculates the desired trajectory point after a time horzion of 10 seconds,
        psi_d=np.arctan2(Y_d-Y,X_d-X)               #calculates the desired heading angle
        
        # to calculate control inputs (F, delta). 
        
        

        # ---------------|Lateral Controller|-------------------------
        P = np.array([-20, -10, -1, 0])
        K = signal.place_poles(A, B, P)
        K=K.gain_matrix
        
        psi_error=wrapToPi(psi-psi_d)
         #calculates the steering angle using PID controller
        e1=np.sqrt(np.power(X_d-X,2)+np.power(Y_d-Y,2))
        e2=psi_error
        e1dot=((X-X_d)*(xdot*np.cos(psi)-ydot*np.sin(psi))+((Y-Y_d)*(xdot*np.sin(psi)+ydot*np.cos(psi))))/e1
        e2dot=psidot
        
        
        


        e = np.hstack((e1, e1dot, e2, e2dot))

        delta = -np.matmul(K,e)
        delta=float(delta)
        # ---------------|Longitudinal Controller|-------------------------
        velocity_error = 30-xdot
        F = self.PID(velocity_error, 50, 0.0001, 0.0001)


        
        

        # Return all states and calculated control inputs (F, delta)
        return X_d, Y_d, xdot, ydot, psi, psidot, F, delta
