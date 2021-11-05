#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
Class definition for physical atribute
"""
import numpy as np

globals()['G'] = 6.67e-11 #Gravitational constant in SI units
globals()['Ms'] = 2e30 #Solar mass in kg
globals()['au'] = 1.5e11 #Astronomical unit in m

class Body:

    def __init__(self, mass, position, velocity):
        self.m = mass
        self.q = position
        self.v = velocity
        self.p = velocity*mass

    def __repr__(self): # Called upon "print(body)"
        return "Body of mass: {0:.2f}kg, position: {1}, velocity: {2}".format(self.m, self.p, self.v)

    def __str__(self): # Called upon "str(body)"
        return "Body of mass: {0:.2f}kg, position: {1}, velocity: {2}".format(self.m, self.p, self.v)

class System:

    def __init__(self, bodylist):
        self.bodylist = bodylist
    
    def get_masses(self): #return the masses of each object
        return np.array([body.m for body in self.bodylist])
    
    def get_positions(self): #return the positions of the bodies
        return np.array([body.q for body in self.bodylist])
    
    def get_velocities(self): #return the positions of the bodies
        return np.array([body.v for body in self.bodylist])
    
    def get_momenta(self): #return the momenta of the bodies
        return np.array([body.p for body in self.bodylist])

    def Mass(self): #return total system mass
        mass = 0
        for body in self.bodylist:
            mass = mass + body.m
        return mass

    def COM(self): #return center of mass in cartesian np_array
        coord = np.zeros(3)
        for body in self.bodylist:
            coord = coord + body.m*body.q
        coord = coord/self.Mass()
        return coord

    def COMV(self): #return center of mass velocity in cartesian np_array
        coord = np.zeros(3)
        for body in self.bodylist:
            coord = coord + body.p
        coord = coord/self.Mass()
        return coord

    def COMShift(self): #Shift coordinates of bodies in system to COM frame and set COM at rest
        for body in self.bodylist:
            body.q = body.q - self.COM()
            body.p = body.p - self.COMV()
        return 0

    def Lval(self): #return angular momentum of bodies in system
        L = np.zeros(3)
        for body in self.bodylist:
            L = L + np.cross(body.q,body.p)
        return L

    def Eval(self): #return total energy of bodies in system
        T = 0
        W = 0
        for body in self.bodylist:
            T = T + 1./2.*body.m*np.linalg.norm(body.v)**2
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.q-otherbody.q)
                    W = W - G*body.m*otherbody.m/rij
        return T + W
    
    def __repr__(self): # Called upon "print(system)"
        return str([print(body) for body in self.bodylist])

    def __str__(self): # Called upon "str(system)"
        return str([str(body) for body in self.bodylist])
