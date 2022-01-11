#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
Class definition for physical attribute
"""
from os import system
import numpy as np
from astropy.coordinates import Angle
from astropy import units as u
from lib.plots import DynamicUpdate
from lib.units import *

class Body:

    def __init__(self, mass, position, velocity):
        self.m = mass
        self.q = position
        self.v = velocity
        self.qb = position
        self.vb = velocity
        self.a = np.zeros(3,dtype=np.longdouble)
        self.ap = np.zeros(3,dtype=np.longdouble)
        self.j = np.zeros(3,dtype=np.longdouble)
        self.jp = np.zeros(3,dtype=np.longdouble)
        self.qp = np.zeros(3,dtype=np.longdouble)
        self.vp = np.zeros(3,dtype=np.longdouble)

    def __repr__(self): # Called upon "print(body)"
        return r"Body of mass: {0:.1e} $M_\odot$, position: {1}, velocity: {2}".format(self.m, self.q, self.v)

    def __str__(self): # Called upon "str(body)"
        return r"Body of mass: {0:.1e} $M_\odot$".format(self.m)

    @property
    def p(self):
        return self.v*self.m

    @property
    def pb(self):
        return self.vb*self.m

class System(Body):

    def __init__(self, bodylist, main = False, blackstyle=True):
        self.blackstyle = blackstyle #for dark mode in plot
        self.bodylist = np.array(bodylist)
        if main == True :
            self.COMShift()
        self.time = 0 #lifetime of system
        self.m = self.M
        self.q = self.COM
        self.v = self.COMV
        self.coordarray = []

    def __repr__(self):  # Called upon "print(system)"
        return str([print(body) for body in self.bodylist])

    def __str__(self):  # Called upon "str(system)"
        return str([str(body) for body in self.bodylist])


    def get_masses(self): #return the masses of each object
        return np.array([body.m for body in self.bodylist],dtype=np.longdouble)


    def get_positions(self): #return the positions of the bodies
        xdata = np.array([body.q[0] for body in self.bodylist],dtype=np.longdouble)
        ydata = np.array([body.q[1] for body in self.bodylist],dtype=np.longdouble)
        zdata = np.array([body.q[2] for body in self.bodylist],dtype=np.longdouble)
        return xdata, ydata, zdata


    def get_velocities(self): #return the positions of the bodies
        vxdata = np.array([body.v[0] for body in self.bodylist],dtype=np.longdouble)
        vydata = np.array([body.v[1] for body in self.bodylist],dtype=np.longdouble)
        vzdata = np.array([body.v[2] for body in self.bodylist],dtype=np.longdouble)
        return vxdata, vydata, vzdata


    def get_momenta(self): #return the momenta of the bodies
        pxdata = np.array([body.p[0] for body in self.bodylist],dtype=np.longdouble)
        pydata = np.array([body.p[1] for body in self.bodylist],dtype=np.longdouble)
        pzdata = np.array([body.p[2] for body in self.bodylist],dtype=np.longdouble)
        return pxdata, pydata, pzdata

    @property
    def M(self): #return total system mass
        mass = 0
        for body in self.bodylist:
            mass = mass + body.m
        return mass

    @property
    def mu(self):
        prod = 1
        for body in self.bodylist:
            prod = prod * body.m
        mu = prod/self.M
        return mu

    @property
    def COM(self): #return center of mass in cartesian np_array
        coord = np.zeros(3,dtype=np.longdouble)
        for body in self.bodylist:
            coord = coord + body.m*body.q
        coord = coord/self.M
        return coord

    @property
    def COMV(self): #return center of mass velocity in cartesian np_array
        coord = np.zeros(3,dtype=np.longdouble)
        for body in self.bodylist:
            coord = coord + body.m*body.v
        coord = coord/self.M
        return coord

    def COMShift(self): #Shift coordinates of bodies in system to COM frame and set COM at rest
        COM = self.COM
        COMV = self.COMV
        for body in self.bodylist:
            body.q = body.q - COM
            body.v = body.v - COMV

    def COMShiftBin(self): #Shift coordinates of inner binary system to COM frame and set COM at rest
        COM = self.COM
        COMV = self.COMV
        for body in self.bodylist:
            body.qb = body.qb - COM
            body.vb = body.vb - COMV

    @property
    def LBIN(self): #return angular momentum of inner binary
        self.COMShiftBin()
        L = np.zeros(3,dtype=np.longdouble)
        for body in self.bodylist:
            L = L + np.cross(body.qb,body.pb)
        return L

    @property
    def EBIN(self): #return total energy of inner binary
        self.COMShiftBin()
        T = 0
        W = 0
        for body in self.bodylist:
            T = T + 1./2.*body.m*np.linalg.norm(body.vb)**2
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.qb-otherbody.qb)
                    W = W - Ga*body.m*otherbody.m/rij
        E = T + W
        return E

    @property
    def LCOM(self): #return angular momentum of the center of mass
        LCOM = np.zeros(3,dtype=np.longdouble)
        dr = self.bodylist[0].m/self.mu*self.bodylist[0].q
        dv = self.bodylist[0].m/self.mu*self.bodylist[0].v
        LCOM = self.mu*np.cross(dr,dv)

        LCOM = self.L

        return LCOM

    @property
    def ECOM(self): #return mechanical energy of the center of mass
        dr = self.bodylist[0].m/self.mu*self.bodylist[0].q
        dv = self.bodylist[0].m/self.mu*self.bodylist[0].v
        ECOM = self.mu/2.*np.linalg.norm(dv)**2 - Ga*self.M*self.mu/np.linalg.norm(dr)

        ECOM = self.E

        return ECOM

    @property
    def L(self): #return angular momentum of bodies in system
        L = np.zeros(3,dtype=np.longdouble)
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
                    W = W - Ga*body.m*otherbody.m/rij
        E = T + W
        return E





    @property
    def ecc(self): #exentricity of two body sub system
        if len(self.bodylist) == 2 :
            ecc = (2.*self.EBIN*(np.linalg.norm(self.LBIN)**2))/((Ga**2)*(self.M**2)*(self.mu**3)) + 1.
        else :
            ecc = np.nan
        return ecc

    @property
    def sma(self): #semi major axis of two body sub system
        if len(self.bodylist) == 2 :
            sma = -Ga*self.M*self.mu/(2.*self.EBIN)
        else :
            sma = np.nan
        return sma

    @property
    def phi(self,body1,body2): #return angle in degree between plans formed by body1 and body2 (perurbator) trajectories
        if len(self.bodylist) == 3 :
            body1 = self.bodylist[0]
            body2 = self.bodylist[2]
            n1 = np.cross(body1.q, body1.v)
            n2 = np.cross(body2.q, body2.v)
            phi = np.arccos(np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2)))
            phi = Angle(phi, u.radian)
            phi = phi.dec
        else :
            phi = np.nan
        return phi
"""""
    def update_coordarray(self): #add current positions of bodies in system in coordarray array.
        sub_array = []
        for body in self.bodylist:
            sub_array.append(body.q)
        self.coordarray.append(sub_array)

    def orbital_analysis(self): #derive semi major axis and eccentricity evolution.
"""""