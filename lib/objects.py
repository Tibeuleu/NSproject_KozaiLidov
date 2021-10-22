#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
Class definition for physical atribute
"""
import numpy as np


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
        comcoord = self.COM
        for body in self.bodylist:
            body.q = body.q-self.COM()
            body.p = body.p - self.COMV()
        return 0

    def Lval(self,Lbodylist): #return angular momentum of bodies in bodylist
        comcoord = np.zeros(3)
        for body in Lbodylist:
            comcoord = comcoord + body.m*body.q
        comcoord = comcoord/self.Mass()
        comq = np.zeros((len(Lbodylist),3))
        i = 0
        L = np.zeros(3)
        for body in Lbodylist:
            comq[i] = body.q-comcoord
            L = L + np.cross(comq[i],body.p)
            i = i+1
        return L

    def Eval(self,Lbodylist): #return total energy of bodies in bodylist
        G = 1. #Gravitational constant (here normalized)
        T = 0
        W = 0
        for body in Lbodylist:
            T = T + 1./2.*body.m*np.linalg.norm(body.v)**2
            for otherbody in Lbodylist:
                if body != otherbody:
                    rij = np.linalg.norm(body.q-otherbody.q)
                    W = W - G*body.m*otherbody.m/rij
        return T + W


if __name__ == "__main__":
    # initialisation mass
    m1 = 10
    m2 = 1
    m3 = 1

    # initialisation position
    q1 = np.array([0, 0, 0])
    q2 = np.array([1, 0, 0])
    q3 = np.array([2, 0, 0])

    # initialisation velocity
    v1 = np.array([0, 0, 0])
    v2 = np.array([1, 1, 0])
    v3 = np.array([2, 0, 0])


    star1 = Body(m1,q1,v1)
    star2 = Body(m2,q2,v2)
    star3 = Body(m3,q3,v3)

    Lbodylist = [star1,star2]

    array = np.zeros((len(Lbodylist),3))
    array[0]=star3.q


    tribody = System([star1,star2,star3])

    print("list=",Lbodylist)

    print(tribody.Lval(Lbodylist))
