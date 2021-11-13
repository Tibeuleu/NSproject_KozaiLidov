#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
Implementation of the plotting and visualization functions.
"""
import numpy as np
import time
import matplotlib.pyplot as plt
from lib.units import *

class DynamicUpdate():
    #Suppose we know the x range
    min_x = -10
    max_x = 10

    plt.ion()

    def __init__(self, dyn_syst):
        self.dyn_syst = dyn_syst

    def set_lims(self, factor=1.5):
        self.ax.set_xlim(factor*self.min_x, factor*self.max_x)
        self.ax.set_ylim(factor*self.min_x, factor*self.max_x)
        self.ax.set_zlim(factor*self.min_x, factor*self.max_x)
    
    def set_blackstyle(self):
        self.fig = plt.figure(figsize=(10,10), facecolor='k')
        self.ax = self.fig.add_subplot(projection='3d')
        self.ax.set_facecolor('k')
        self.ax.xaxis.label.set_color('w')
        self.ax.yaxis.label.set_color('w')
        self.ax.zaxis.label.set_color('w')
        self.ax.tick_params(axis='x',colors='w')
        self.ax.tick_params(axis='y',colors='w')
        self.ax.tick_params(axis='z',colors='w')
        self.ax.w_xaxis.line.set_color('w')
        self.ax.w_yaxis.line.set_color('w')
        self.ax.w_zaxis.line.set_color('w')
        self.ax.w_xaxis.set_pane_color((0,0,0,0))
        self.ax.w_yaxis.set_pane_color((0,0,0,0))
        self.ax.w_zaxis.set_pane_color((0,0,0,0))
        

    def launch(self, blackstyle=True):
        #Set up plot
        if blackstyle:
            self.blackstyle = True
            self.set_blackstyle()
        else:
            self.blackstyle = False
            self.fig = plt.figure(figsize=(10,10))
            self.ax = self.fig.add_subplot(projection='3d')

        self.lines = []
        for i,body in enumerate(self.dyn_syst.bodylist):
            x, y, z = body.q
            lines, = self.ax.plot([x],[y],[z],'o',color="C{0:d}".format(i),label="{0:s}".format(str(body)))
            self.lines.append(lines)
        self.lines = np.array(self.lines)
        #Autoscale on unknown axis and known lims on the other
        self.ax.set_autoscaley_on(True)
        self.set_lims()
        #Other stuff
        self.ax.grid()
        if self.blackstyle:
            self.ax.legend(labelcolor='w', frameon=True, framealpha=0.2)
        else:
            self.ax.legend()

    def on_running(self, dyn_syst, step=None, label=None):
        xdata, ydata, zdata = dyn_syst.get_positions()
        values = np.sqrt(np.sum((np.array((xdata,ydata,zdata))**2).T,axis=1))
        self.min_x, self.max_x = -np.max([np.abs(values).max(),self.max_x]), np.max([np.abs(values).max(),self.max_x])
        self.set_lims()
        #Update data (with the new _and_ the old points)
        for i,body in enumerate(dyn_syst.bodylist):
            x, y, z = body.q
            self.lines[i].set_data_3d([x], [y], [z])
        if not label is None:
            if self.blackstyle:
                self.ax.set_title(label,color='w')
            else:
                self.ax.set_title(label)

        #Need both of these in order to rescale
        self.ax.relim()
        self.ax.autoscale_view()
        #We need to draw *and* flush
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        if not step is None and step%1000==0:
            self.fig.savefig("tmp/{0:06d}.png".format(step),bbox_inches="tight")
    
    def close(self):
        plt.close()

    #Example
    def __call__(self):
        import numpy as np
        import time
        self.on_launch()
        xdata = []
        ydata = []
        for x in np.arange(0,10,0.5):
            xdata.append(x)
            ydata.append(np.exp(-x**2)+10*np.exp(-(x-7)**2))
            self.on_running(xdata, ydata)
            time.sleep(1)
        return xdata, ydata

def display_parameters(E,L,parameters,savename=""):
    """
    """
    if savename != "":
        savename += "_"
    duration, step, dyn_syst, integrator = parameters
    bodies = ""
    for body in dyn_syst.bodylist:
        bodies += str(body)+" ; "
    title = "Relative difference of the {} "+"for a system composed of {0:s}\n integrated with {1:s} for a duration of {2:.2f} years with a step of {3:.2e} years.".format(bodies, integrator, duration/yr, step/yr)
    fig1 = plt.figure(figsize=(15,7))
    ax1 = fig1.add_subplot(111)
    ax1.plot(np.arange(E.shape[0])*step/yr, np.abs((E-E[0])/E[0]), label=r"$\left|\frac{\delta E_m}{E_m(t=0)}\right|$")
    ax1.set(xlabel=r"$t (yr)$", ylabel=r"$\left|\frac{\delta E_m}{E_m(t=0)}\right|$", yscale='log')
    ax1.legend()
    fig1.suptitle(title.format("mechanical energy"))
    fig1.savefig("plots/{0:s}dEm.png".format(savename),bbox_inches="tight")
    fig2 = plt.figure(figsize=(15,7))
    ax2 = fig2.add_subplot(111)
    dL = ((L-L[0])/L[0])
    dL[np.isnan(dL)] = 0.
    ax2.plot(np.arange(L.shape[0])*step/yr, np.abs(np.sum(dL,axis=1)), label=r"$\left|\frac{\delta \vec{L}}{\vec{L}(t=0)}\right|$")
    ax2.set(xlabel=r"$t (yr)$", ylabel=r"$\left|\frac{\delta \vec{L}}{\vec{L}(t=0)}\right|$",yscale='log')
    ax2.legend()
    fig2.suptitle(title.format("kinetic moment"))
    fig2.savefig("plots/{0:s}dL2.png".format(savename),bbox_inches="tight")
    plt.show(block=True)
  