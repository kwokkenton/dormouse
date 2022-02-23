#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 18:17:26 2022

@author: go34660
"""
"""
actual code
"""

import meep as mp
import numpy as np
from meep.materials import Au


def ampfunc(v3,k): #amplitude function used to for straight wave propagation
    return np.exp(2*np.pi*1j*(k*v3))

class meep_sim():
    def __init__(self,x,y,theta,geometry,resolution,wavelength,wvl_list,k):
        self.xs = x
        self.ys = y
        self.resolution = resolution
        self.theta = theta
        self.wavelength = wavelength
        self.freq = 1/wavelength
        self.wvl_min = wvl_list[0]
        self.wvl_max = wvl_list[1]
        self.fmax = 1/self.wvl_min
        self.fmin = 1/self.wvl_max
        self.fcen = (self.fmax + self.fmin)/2
        self.df = self.fmax - self.fmin
        self.cell = mp.Vector3(x,y*2)
        self.pml = [mp.PML(1.0)]
        self.nfreq = 100
        self.k = k
        self.geometry = geometry
        def amplitude(v):
            return ampfunc(v,self.k)
        self.source = [mp.Source(mp.ContinuousSource(frequency=self.freq,is_integrated=True), #0.633 = 633nm?
                            component=mp.Ez, #in ez direction not really sure what polarization has on these stuff
                            size = mp.Vector3(self.xs,0), #needs to encompass the whole cell for it to be straight
                            center=mp.Vector3(0,0), 
                            amp_func = amplitude #amp function
                            )]
        self.sim = mp.Simulation(resolution=self.resolution,
                            cell_size=self.cell,
                            boundary_layers=self.pml,
                            sources=self.source,
                            k_point=k,
                            #symmetries=symmetries,
                            geometry=self.geometry)
        self.flob = mp.FluxRegion(center = mp.Vector3(0,-0.8*self.ys),
                            size = mp.Vector3(self.xs,0)
                            )
        self.flux = self.sim.add_flux(self.fcen, self.df, self.nfreq,self.flob)
    def run(self,stop,photos = True,time = 100):
        if photos == True:
            self.sim.run(mp.at_every(time, mp.output_png(mp.Ez, "-Zc dkbluered")),until=stop)
        else:
            self.sim.run(until=stop)
        self.tot_flux = mp.get_fluxes(self.flux)
        self.freq_list = mp.get_flux_freqs(self.flux)
        return [self.freq_list,self.tot_flux]
    def image(self):
        ez_data = self.sim.get_array(center=mp.Vector3(), size=self.cell, component=mp.Ez)
        eps_data = self.sim.get_array(center=mp.Vector3(), size=self.cell, component=mp.Dielectric)
        return ez_data



resolution = 10 #resolution per unit length
scale = 10 #scale to change simulation size
ys = 50*scale #size of simulation (y direction)
xs = 30*scale#size of simulation (x direction)

a = 1e-7 #the scale the simulation works in 
theta = np.pi*43/180 #rotation angle
wavelength = 6.63*scale #10^-7
freq = 1/wavelength
eps = 2.31 #glass

circle_dist = 0.1 #distance between ball and surface

wvl_min = 1*scale
wvl_max = 15*scale

circle_rad = 1

k = mp.Vector3(0,-1,0).rotate(mp.Vector3(0,0,1),theta).scale(freq) * np.sqrt(eps) #k vector which decides the dicrection the wave travels in

offset = 0

extra_length = 5

geometry_scatter = [mp.Block(mp.Vector3(mp.inf,xs + extra_length*scale,mp.inf), #infinite block
                     center = mp.Vector3(0,0 + offset*scale + extra_length*scale/2), #center should be centered
                     material = mp.Medium(epsilon=eps)) #glass
                      ,
            mp.Cylinder(radius = circle_rad*scale,
                    height=0,
                    center = mp.Vector3(0,xs/2 + (circle_rad/2 + circle_dist)*scale + offset*scale + extra_length*scale/2),
                    #length of block + length of block (extra) + offset + radius of circle + extra so they aren't touching
                    material = Au)
                    ]

geometry_empty = [mp.Block(mp.Vector3(mp.inf,xs,mp.inf), #infinite block
                     #e1 = mp.Vector3(np.cos(theta),np.sin(theta)), #decides 2 axis which the block is rotated in
                     #e2 = mp.Vector3(-np.sin(theta),np.cos(theta)),
                     center = mp.Vector3(0,0), #center should be high
                     material = mp.Medium(epsilon=eps)) #glass
             ]

sim = meep_sim(xs,ys,theta,geometry_scatter,resolution,wavelength,[wvl_min,wvl_max],k)

for i in range(6):
    data = sim.run(300,False)
    imag = sim.image()
    np.savetxt('gold_{}nm_{}_{}_extra_{}.txt'.format(circle_rad*100,300*i,300*(i+1),extra_length*scale),data)
    np.savetxt('im_gold_{}nm_{}_{}_extra_{}.txt'.format(circle_rad*100,300*i,300*(i+1),extra_length*scale),imag)


