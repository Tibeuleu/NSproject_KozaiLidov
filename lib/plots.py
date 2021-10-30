#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
Implementation of the plotting and visualization functions.
"""
import numpy as np
import time
import matplotlib.pyplot as plt

class DynamicUpdate():
    #Suppose we know the x range
    min_x = -10
    max_x = 10

    plt.ion()

    def on_launch(self):
        #Set up plot
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(projection='3d')
        self.lines, = self.ax.plot([],[],[],'o')
        #Autoscale on unknown axis and known lims on the other
        self.ax.set_autoscaley_on(True)
        self.ax.set_xlim(self.min_x, self.max_x)
        self.ax.set_ylim(self.min_x, self.max_x)
        self.ax.set_zlim(self.min_x, self.max_x)
        #Other stuff
        self.ax.grid()
        #self.ax.set_aspect('equal')

    def on_running(self, xdata, ydata, zdata, step=None, label=None):
        #Update data (with the new _and_ the old points)
        self.lines.set_data_3d(xdata, ydata, zdata)
        if not label is None:
            self.ax.set_title(label)
        #Need both of these in order to rescale
        self.ax.relim()
        self.ax.autoscale_view()
        #We need to draw *and* flush
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        if not step is None and step%50==0:
            self.fig.savefig("tmp/{0:05d}.png".format(step),bbox_inches="tight")

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
