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
    """
    Class for dynamic display of the integrated system.
    Initialise with the System of Bodies that will be integrated.
    Launch the display, then call on_running method the update it.
    -----
    dyn_syst = System(bodylist)
    d = DynamicUpdate(dyn_syst)
    [Start integration procedure]
    #Init
    d.launch()
    [in intgration loop]
    #integration
    for step in range(duration):
        #update attributes of dyn_syst
        d.on_running(step=step, label="will be displayed on current update")
    [when integration reach end]
    d.close()
    -----
    Additionnal parameters:
    launch(blackstyle): boolean
        If the display should have black background.
        Default to True.
    launch(lim_factor): float
        Limits of the 3D display are dynamically updated to max_value*lim_factor.
        Should always be >1. to have all bodies in the display.
        Default to 1.5
    """
    #Suppose we know the x range
    min_x = -1
    max_x = 1

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

    def launch(self, blackstyle=True, lim_factor=1.5):
        #Set up plot
        if blackstyle:
            self.blackstyle = True
            self.set_blackstyle()
        else:
            self.blackstyle = False
            self.fig = plt.figure(figsize=(10,10))
            self.ax = self.fig.add_subplot(projection='3d')
        self.lim_factor = 1.5

        self.lines = []
        for i,body in enumerate(self.dyn_syst.bodylist):
            x, y, z = body.q/au-self.dyn_syst.COM/au
            lines, = self.ax.plot([x],[y],[z],'o',color="C{0:d}".format(i),label="{0:s}".format(str(body)))
            self.lines.append(lines)
        self.lines = np.array(self.lines)
        #Autoscale on unknown axis and known lims on the other
        self.ax.set_autoscaley_on(True)
        self.set_lims(factor=self.lim_factor)
        #Other stuff
        self.ax.grid()
        if self.blackstyle:
            self.ax.legend(labelcolor='w', frameon=True, framealpha=0.2)
            self.ax.set_xlabel('AU', color='w')
            self.ax.set_ylabel('AU', color='w')
            self.ax.set_zlabel('AU', color='w')
        else:
            self.ax.legend()
            self.ax.set_xlabel('AU')
            self.ax.set_ylabel('AU')
            self.ax.set_zlabel('AU')

    def on_running(self, step=None, label=None):
        xdata, ydata, zdata = self.dyn_syst.get_positionsCOM()
        values = np.sqrt(np.sum((np.array((xdata,ydata,zdata))**2).T,axis=1))/au
        self.min_x, self.max_x = -np.max([np.abs(values).max(),self.max_x]), np.max([np.abs(values).max(),self.max_x])
        self.set_lims(factor=self.lim_factor)
        #Update data (with the new _and_ the old points)
        for i,body in enumerate(self.dyn_syst.bodylist):
            x, y, z = body.q/au
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
        if not step is None and step%50==0:
            self.fig.savefig("tmp/{0:06d}.png".format(step),bbox_inches="tight")
    
    def close(self):
        plt.close()


def display_parameters(E,L,sma,ecc,phi,parameters,savename="",display_param=True):
    """
    Save integrated parameters plots to multiple png files.
    -----
    Inputs:
    E, L, sma, ecc, phi: list of np.ndarray
        list of integrated parameters value computed for each step length in
        parameters[step] list.
    parameters : list
        list of simulation parameters : duration, steps, system, integrator,
        and initialisation parameters.        
    savename : str
        default savename that will be prepend to each saved file path.
        Default to empty string.
    display_param : boolean
        Set to True if the user wants to display the plots before saving,
        False otherwise. Default to True.
    """
    if savename != "":
        savename += "_"
    duration, step, dyn_syst, integrator, init = parameters
    a, e, psi = init
    bodies = ""
    init_str = ""
    for i, body in enumerate(dyn_syst.bodylist):
        bodies += str(body)+" ; "
        init_str += r"a{0:d} = {1:.2f} au, e{0:d} = {2:.2f}, $\psi${0:d} = {3:.2f}° ; ".format(i+1,a[i]/au,e[i],psi[i]*180./np.pi)
    title1, title2 = "Relative difference of the {0:s} ","for a system composed of {0:s}\n integrated with {1:s} for a duration of {2:.2f} years with initial parameters\n {3:s}".format(bodies, integrator, duration/yr, init_str)

    fig1 = plt.figure(figsize=(15,7))
    ax1 = fig1.add_subplot(111)
    for i in range(len(E)):
        ax1.plot(np.arange(E[i].shape[0]-1)*step[i]/yr, np.abs((E[i][1:]-E[i][0])/E[i][0]), label="step of {0:.2e}s".format(step[i]))
    ax1.set(xlabel=r"$t \, [yr]$", ylabel=r"$\left|\frac{\delta E_m}{E_m(t=0)}\right|$", yscale='log')
    ax1.legend()
    fig1.suptitle(title1.format("mechanical energy")+title2)
    fig1.savefig("plots/{0:s}dEm.png".format(savename),bbox_inches="tight")

    fig2 = plt.figure(figsize=(15,7))
    ax2 = fig2.add_subplot(111)
    for i in range(len(L)):
        dL = ((L[i]-L[i][0])/L[i][0])
        dL[np.isnan(dL)] = 0.
        ax2.plot(np.arange(L[i].shape[0]-1)*step[i]/yr, np.abs(np.sum(dL[1:],axis=1)), label="step of {0:.2e}s".format(step[i]))
    ax2.set(xlabel=r"$t \, [yr]$", ylabel=r"$\left|\frac{\delta \vec{L}}{\vec{L}(t=0)}\right|$",yscale='log')
    ax2.legend()
    fig2.suptitle(title1.format("kinetic moment")+title2)
    fig2.savefig("plots/{0:s}dL2.png".format(savename),bbox_inches="tight")

    fig3 = plt.figure(figsize=(15,7))
    ax3 = fig3.add_subplot(111)
    for i in range(len(E)):
        ax3.plot(np.arange(E[i].shape[0])*step[-1]/yr, E[i], label="step of {0:.2e}s".format(step[i]))
    ax3.set(xlabel=r"$t \, [yr]$", ylabel=r"$E \, [J]$")
    ax3.legend()
    fig3.suptitle("Mechanical energy of the whole system "+title2)
    fig3.savefig("plots/{0:s}E.png".format(savename),bbox_inches="tight")
    
    fig4 = plt.figure(figsize=(15,7))
    ax4 = fig4.add_subplot(111)
    for i in range(len(L)):
        L2 = np.array([np.linalg.norm(Li)**2 for Li in L[i]])
        ax4.plot(np.arange(L[i].shape[0])*step[i]/yr, L2, label=r"$L^2$ for step of {0:.2e}s".format(step[i]))
    ax4.set(xlabel=r"$t \, [yr]$", ylabel=r"$\left|\vec{L}\right|^2 \, [kg^2 \cdot m^4 \cdot s^{-2}]$",yscale='log')
    ax4.legend()
    fig4.suptitle("Squared norm of the kinetic moment of the whole system "+title2)
    fig4.savefig("plots/{0:s}L.png".format(savename),bbox_inches="tight")

    fig5 = plt.figure(figsize=(15,7))
    ax5 = fig5.add_subplot(111)
    ax5.plot(np.arange(sma[-1].shape[0])*step[-1]/yr, sma[-1]/au, label="a (step of {0:.2e}s)".format(step[-1]))
    ax5.plot(np.arange(ecc[-1].shape[0])*step[-1]/yr, ecc[-1], label="e (step of {0:.2e}s)".format(step[-1]))
    ax5.set(xlabel=r"$t \, [yr]$", ylabel=r"$a \, [au] \, or \, e$")
    ax5.legend()
    fig5.suptitle("Semi major axis and eccentricity "+title2)
    fig5.savefig("plots/{0:s}a_e.png".format(savename),bbox_inches="tight")

    fig6 = plt.figure(figsize=(15,7))
    ax6 = fig6.add_subplot(111)
    for i in range(len(phi)):
        ax6.plot(np.arange(phi[i].shape[0])*step[-1]/yr, phi[i], label="step of {0:.2e}s".format(step[i]))
    ax6.set(xlabel=r"$t \, [yr]$", ylabel=r"$\phi \, [^{\circ}]$")
    ax6.legend()
    fig6.suptitle("Inclination angle of the perturbator's orbital plane "+title2)
    fig6.savefig("plots/{0:s}phi.png".format(savename),bbox_inches="tight")
    
    if display_param:
        plt.show(block=True)
  