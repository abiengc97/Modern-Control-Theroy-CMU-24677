# Fill in the respective functions to implement the controller

# Import libraries
import numpy as np
from base_controller import BaseController
from scipy import signal, linalg
from scipy.signal import StateSpace, lsim, dlsim
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
        self.integralPsiError = 0
        self.previousPsiError = 0
        self.previousXdotError = 0

        # Add additional member variables according to your need here.
        
    def dlqr(self,A,B,Q,R):
        S = np.matrix(linalg.solve_discrete_are(A, B, Q, R))
        K = -np.matrix(linalg.inv(B.T@S@B+R)@(B.T@S@A))
        return K
    def update(self, timestep):

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
        A = np.array([[0, 1, 0, 0], [0, -4*Ca / (m * xdot), 4*Ca / m, -(2*Ca*(lf - lr))/(m*xdot)], [0, 0, 0, 1], [0, -(2*Ca*(lf - lr)) / (Iz * xdot), (2*Ca*(lf - lr)) / Iz, (-2*Ca*(np.power(lf, 2) + np.power(lr,   2))) / (Iz * xdot)]])
        B = np.array([[0], [2*Ca / m], [0], [(2 * Ca* lf) / Iz]])
        C = np.eye(4)
        D = np.zeros((4,1))
        
        # Find the closest node to the vehicle
        _, node = closestNode(X, Y, trajectory)
# Choose a node that is ahead of our current node based on index
        forwardIndex = 100
        
        
        # Remember, your controllers will need to use the states
        # to calculate control inputs (F, delta). 

        # ---------------|Lateral Controller|-------------------------
        try:
            psiDesired=np.arctan2(trajectory[node+forwardIndex,1]-trajectory[node,1],trajectory[node+forwardIndex,0]-trajectory[node,0])
            e1 = (Y - trajectory[node+forwardIndex,1])*np.cos(psiDesired) -(X - trajectory[node+forwardIndex,0])*np.sin(psiDesired)
        except:
            psiDesired = np.arctan2(trajectory[-1,1]-trajectory[node,1],trajectory[-1,0]-trajectory[node,0])
            e1 = (Y - trajectory[node,1])*np.cos(psiDesired) -(X - trajectory[node,0])*np.sin(psiDesired)

        e1dot=ydot+xdot*wrapToPi(psi-psiDesired)
        e2=wrapToPi(psi-psiDesired)
        e2dot=psidot
        e=np.hstack((e1,e1dot,e2,e2dot))
        sys_ct=StateSpace(A,B,C,D)
        sys_dt=sys_ct.to_discrete(delT)
        dt_A = sys_dt.A
        dt_B = sys_dt.B
        R=100
        Q=np.array([[10,0,0,0], [0,10,0,0], [0,0,300,0], [0,0,0,400]])
        K= self.dlqr(dt_A,dt_B,Q,R)
        delta=float(np.matmul(K,e))

        delta=clamp(delta,-np.pi/6,np.pi/6)

        # ---------------|Longitudinal Controller|-------------------------
        kp = 1000000
        ki = 9
        kd = 20
        if abs(delta) > 0.1:
            desiredVelocity = 9
        else:
            desiredVelocity = 11
        xdotError = (desiredVelocity - xdot)
        self.integralXdotError += xdotError
        derivativeXdotError = xdotError - self.previousXdotError
        self.previousXdotError = xdotError
        F = kp*xdotError + ki*self.integralXdotError*delT +kd*derivativeXdotError/delT
        F=clamp(F,0,15736)
        # Return all states and calculated control inputs (F, delta)
        return X, Y, xdot, ydot, psi, psidot, F, delta

