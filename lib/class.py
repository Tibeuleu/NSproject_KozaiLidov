import numpy as np

"""
Class definition for physical atribute
"""


class Body:


    
    def __init__(self, mass, position, velocity):
        self.m = mass
        self.q = position
        self.v = velocity
        self.p = velocity*mass


class System:

    def __init__(self, bodylist):
        self.bodylist = bodylist

    def Mass(self): #return total system mass
        mass = 0
        for body in self.bodylist:
            mass = mass + body.m
        return mass

    def COM(self): #return center of mass in cartesian np_array
        coord = np.zeros(3)
        for body in self.bodylist:
            coord = coord + body.m*body.q
        coord = coord/self.Mass
        return coord

    def Lval(self,Lbodylist): #return angular momentum of body in bodylist
        comcoord = np.zeros(3)
        for body in Lbodylist:
            comcoord = comcoord + body.m*body.q
        comcoord = comcoord/self.Mass
        comq = np.zeros((len(Lbodylist),3))
        i = 0
        L = 0
        for body in Lbodylist:
            comq[i] = body.q-comcoord
            L = L + np.cross(comq[i],body.p)
            i = i+1
        return L


    #def initialize(self):

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
v2 = np.array([1, 0, 0])
v3 = np.array([2, 0, 0])


star1 = Body(m1,q1,v1)
print('test')
star2 = Body(m2,q2,v2)
star3 = Body(m3,q3,v3)

Lbodylist = [star1,star2]

array = np.zeros((len(Lbodylist),3))
array[0]=star3.q
print(array[0])

star2 = Body

