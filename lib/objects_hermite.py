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
        self.p = velocity*mass
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

    def Update_a(self): #update acceleration of bodies in system
        for body in self.bodylist:
            body.a = np.zeros(3)
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.q-otherbody.q)
                    body.a = body.a - (body.q-otherbody.q)*G*otherbody.m/(rij**3)

    def Update_j(self): #update jerk of bodies in system
        for body in self.bodylist:
            body.j = np.zeros(3)
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.q-otherbody.q)
                    deltav = (body.v-otherbody.v)
                    deltar = (body.q-otherbody.q)
                    vr = deltav + 3.*deltar*np.inner(deltav,deltar)/(rij**2)
                    body.j = body.j - G*otherbody.m/(rij**3)*vr

    def Predict(self,dt):  # update predicted position and velocities of bodies in system
        for body in self.bodylist:
            body.qp = body.q +dt*body.v+((dt**2)*body.a/2.)+((dt**3)*body.j/6.)
            body.vp = body.v + dt*body.a + ((dt**2)*body.j/2.)

    def Update_ap(self): #update acceleration of bodies in system
        for body in self.bodylist:
            body.ap = np.zeros(3)
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.qp-otherbody.qp)
                    body.ap = body.ap - (body.qp-otherbody.qp)*G*otherbody.m/(rij**3)

    def Update_jp(self): #update jerk of bodies in system
        for body in self.bodylist:
            body.jp = np.zeros(3)
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.qp-otherbody.qp)
                    deltav = (body.vp-otherbody.vp)
                    deltar = (body.qp-otherbody.qp)
                    vr = deltav + 3.*deltar*np.inner(deltav,deltar)/(rij**2)
                    body.jp = body.jp - G*otherbody.m/(rij**3)*vr

    def Correct(self,dt):  # correct position and velocities of bodies in system
        for body in self.bodylist:
            a2 = (6.*(body.a-body.ap)+dt*(4*body.j+2*body.jp))/(dt**2)
            a3 = (12. * (body.a - body.ap) + dt * 6. * (body.j + body.jp)) / (dt ** 3)

            body.q = body.qp +((dt**4)*a2/24.) + ((dt**5)*a3/120.)
            body.v = body.vp +((dt**3)*a2/6.) + ((dt**4)*a3/24.)

    def HPC(self, dt):  # update position and velocities of bodies in system with hermite predictor corrector
        self.COMShift()
        self.Update_a()
        self.Update_j()
        self.Predict(dt)
        self.Update_ap()
        self.Update_jp()
        self.Correct(dt)
        self.time = self.time + dt
        for body in self.bodylist:
            body.p = body.v*body.m
    
    def hermite(self, duration, dt, recover_param=False, display=False, savename=None):
        if display:
            try:
                system("mkdir tmp")
            except IOError:
                system("rm tmp/*")
            d = DynamicUpdate(self)
            d.launch(self.blackstyle)

        N = np.ceil(duration/dt).astype(int)
        E = np.zeros(N)
        L = np.zeros((N,3))
        for j in range(N):
            self.HPC(dt)

            E[j] = self.E
            L[j] = self.L

            if display and j%100==0:
                # display progression
                if len(self.bodylist) == 1:
                    d.on_running(self, step=j, label="step {0:d}/{1:d}".format(j,N))
                else:
                    d.on_running(self, step=j, label="step {0:d}/{1:d}".format(j,N))
        if display:
            d.close()
            if not savename is None:
                system("convert -delay 5 -loop 0 tmp/??????.png tmp/temp.gif && rm tmp/??????.png")
                system("convert tmp/temp.gif -fuzz 10% -layers Optimize plots/{0:s}_dynsyst.gif".format(savename))
     
        if recover_param:
            return E, L
 
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

    def __repr__(self): # Called upon "print(system)"
        return str([print(body) for body in self.bodylist])

    def __str__(self): # Called upon "str(system)"
        return str([str(body) for body in self.bodylist])
