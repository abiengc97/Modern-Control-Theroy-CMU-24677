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

        # ---------------|Lateral Controller|-------------------------
        # Remember, your controllers will need to use the states
        
        time_horizon=10
        _,nearidx=closestNode(X,Y,trajectory)
        [X_d,Y_d]=trajectory[nearidx+time_horizon]  #calculates the desired trajectory point after a time horzion of 10 seconds,
        psi_d=np.arctan2(Y_d-Y,X_d-X)               #calculates the desired heading angle
        
        # to calculate control inputs (F, delta). 


        # ---------------|Lateral Controller|-------------------------
        psi_error=wrapToPi(psi_d-psi)
        delta=self.PID(psi_error,5,0,0) #calculates the steering angle using PID controller
        delta=clamp(delta,-np.pi/6,np.pi/6)   #clamps the steering angle between -30 and 30 degrees

        # ---------------|Longitudinal Controller|-------------------------
        position_error=np.sqrt((X_d-X)**2+(Y_d-Y)**2)
        F=self.PID(position_error,500,0,0)
        F=clamp(F,0,15736) #clamps the throttle input between 0 and 15736 N


        
        

        # Return all states and calculated control inputs (F, delta)
        return X_d, Y_d, xdot, ydot, psi, psidot, F, delta
