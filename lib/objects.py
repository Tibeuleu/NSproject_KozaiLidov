#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
Class definition for physical attribute
"""
from os import system
import numpy as np
from lib.plots import DynamicUpdate
from lib.units import *

class Body:

    def __init__(self, mass, position, velocity):
        self.m = mass
        self.q = position
        self.v = velocity
        self.p = velocity * mass
        self.a = np.zeros(3)
        self.ap = np.zeros(3)
        self.j = np.zeros(3)
        self.jp = np.zeros(3)
        self.qp = np.zeros(3)
        self.vp = np.zeros(3)

    def __repr__(self): # Called upon "print(body)"
        return r"Body of mass: {0:.2f} $M_\odot$, position: {1}, velocity: {2}".format(self.m/Ms, self.q, self.v)

    def __str__(self): # Called upon "str(body)"
        return r"Body of mass: {0:.2f} $M_\odot$".format(self.m/Ms)

class System(Body):

    def __init__(self, bodylist, blackstyle=True):
        self.blackstyle = blackstyle #for dark mode in plot
        self.bodylist = np.array(bodylist)
        self.time = 0 #lifetime of system
        self.m = self.M
        self.q = self.COM
        self.v = self.COMV

    def __repr__(self):  # Called upon "print(system)"
        return str([print(body) for body in self.bodylist])

    def __str__(self):  # Called upon "str(system)"
        return str([str(body) for body in self.bodylist])

    @property
    def get_masses(self): #return the masses of each object
        return np.array([body.m for body in self.bodylist])
    
    @property
    def get_positions(self): #return the positions of the bodies
        xdata = np.array([body.q[0] for body in self.bodylist])
        ydata = np.array([body.q[1] for body in self.bodylist])
        zdata = np.array([body.q[2] for body in self.bodylist])
        return xdata, ydata, zdata
    
    @property
    def get_velocities(self): #return the positions of the bodies
        vxdata = np.array([body.v[0] for body in self.bodylist])
        vydata = np.array([body.v[1] for body in self.bodylist])
        vzdata = np.array([body.v[2] for body in self.bodylist])
        return vxdata, vydata, vzdata
    
    @property
    def get_momenta(self): #return the momenta of the bodies
        pxdata = np.array([body.p[0] for body in self.bodylist])
        pydata = np.array([body.p[1] for body in self.bodylist])
        pzdata = np.array([body.p[2] for body in self.bodylist])
        return pxdata, pydata, pzdata

    @property
    def M(self): #return total system mass
        mass = 0
        for body in self.bodylist:
            mass = mass + body.m
        return mass

    @property
    def COM(self): #return center of mass in cartesian np_array
        coord = np.zeros(3)
        for body in self.bodylist:
            coord = coord + body.m*body.q
        coord = coord/self.M
        return coord

    @property
    def COMV(self): #return center of mass velocity in cartesian np_array
        coord = np.zeros(3)
        for body in self.bodylist:
            coord = coord + body.p
        coord = coord/self.M
        return coord

    def COMShift(self): #Shift coordinates of bodies in system to COM frame and set COM at rest
        for body in self.bodylist:
            body.q = body.q - self.COM
            body.p = body.p - self.COMV

    @property
    def L(self): #return angular momentum of bodies in system
        L = np.zeros(3)
        for body in self.bodylist:
            L = L + np.cross(body.q,body.p)
        return L

    @property
    def E(self): #return total energy of bodies in system
        T = 0
        W = 0
        for body in self.bodylist:
            T = T + 1./2.*body.m*np.linalg.norm(body.v)**2
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.q-otherbody.q)
                    W = W - G*body.m*otherbody.m/rij
        E = T + W
        return E
         
    @property
    def mu(self):
        sum = 0
        prod = 1
        for body in self.bodylist:
            prod = prod * body.m
        mu = prod/self.M
        return mu

    @property
    def ex(self): #exentricity of system (if composed of 2 bodies)
        if len(self.bodylist) != 2 :
            return np.nan
        else:
            k = (2.*self.E*(np.linalg.norm(self.L)**2))/((G**2)*(self.M**2)*(self.mu**3)) + 1.
            return k

    @property
    def sma(self): #semi major axis of system (if composed of 2 bodies)
        if len(self.bodylist) != 2 :
            return np.nan
        else:
            sma = -G*self.M*self.mu/(2.*self.E)
            return sma
