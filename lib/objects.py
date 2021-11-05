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
        self.a = np.zeros(3)
        self.ap = np.zeros(3)
        self.j = np.zeros(3)
        self.jp = np.zeros(3)
        self.qp = np.zeros(3)
        self.vp = np.zeros(3)

    def __repr__(self): # Called upon "print(body)"
        return "Body of mass: {0:.2f}kg, position: {1}, velocity: {2}".format(self.m, self.p, self.v)

    def __str__(self): # Called upon "str(body)"
        return "Body of mass: {0:.2f}kg, position: {1}, velocity: {2}".format(self.m, self.p, self.v)

class System:

    def __init__(self, bodylist):
        self.bodylist = np.array(bodylist)
        self.time = 0
    
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

    def Update_a(self): #update acceleration of bodies in system
        for body in self.bodylist:
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.q-otherbody.q)
                    body.a = body.a - (body.q-otherbody.q)*G*otherbody.m/(rij**2)
        return 1

    def Update_j(self): #update jerk of bodies in system
        for body in self.bodylist:
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.q-otherbody.q)
                    deltav = (body.v-otherbody.v)
                    deltar = (body.q-otherbody.q)
                    vr = deltav + 3.*deltar*np.inner(deltav,deltar)/(rij**2)
                    body.j = body.j - G*otherbody.m/(rij**3)*vr
        return 1

    def Predict(self,dt):  # update predicted position and velocities of bodies in system
        for body in self.bodylist:
            body.qp = body.q +dt*body.v+((dt**2)*body.a/2.)+((dt**3)*body.j/6.)
            body.vp = body.v + dt*body.a + ((dt**2)*body.j/2.)
        return 1

    def Update_ap(self): #update acceleration of bodies in system
        for body in self.bodylist:
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.qp-otherbody.qp)
                    body.ap = body.ap - (body.qp-otherbody.qp)*G*otherbody.m/(rij**2)
        return 1

    def Update_jp(self): #update jerk of bodies in system
        for body in self.bodylist:
            for otherbody in self.bodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.qp-otherbody.qp)
                    deltav = (body.vp-otherbody.vp)
                    deltar = (body.qp-otherbody.qp)
                    vr = deltav + 3.*deltar*np.inner(deltav,deltar)/(rij**2)
                    body.jp = body.jp - G*otherbody.m/(rij**3)*vr
        return 1

    def Correct(self,dt):  # correct position and velocities of bodies in system
        for body in self.bodylist:
            a2 = (6.*(body.a-body.ap)+dt*(4*body.j+2*body.jp))/(dt**2)
            a3 = (12. * (body.a - body.ap) + dt * 6. * (body.j + body.jp)) / (dt ** 3)

            body.q = body.qp +((dt**4)*a2/24.) + ((dt**5)*a3/120.)
            body.v = body.vp +((dt**3)*a2/6.) + ((dt**4)*a3/24.)
        return 1

    def HPC(self, dt):  # update position and velocities of bodies in system with hermite predictor corrector
        self.update_a()
        self.update_j()
        self.predict(dt)
        self.update_ap()
        self.update_jp()
        self.update(dt)
        self.time = self.time + dt

        return 1