#Sim. procedure

from YABGC_M import YoungAdultBornGranuleCell

import neuron
from neuron import h
import numpy as np
from neuron.units import ms, mV
from neuron import h, gui

my_cell = YoungAdultBornGranuleCell()

#simulation setup
stim = h.IClamp(my_cell.soma(0.5))
stim.delay = 100
stim.dur = 1000
stim.amp = 0.015 #nA

###Record of the soma
soma_v = h.Vector().record(my_cell.soma(0.5)._ref_v)

###Record of time

t = h.Vector().record(h._ref_t)

#Basic properties of the simulation. dt, temperature, sim duration and initial voltage

h.celsius = 32
h.tstop = 500
h.v_init = -62


#plotting the result

import bokeh
from bokeh.io import output_notebook
import bokeh.plotting as plt
output_notebook()

f = plt.figure(x_axis_label='t (ms)', y_axis_label='v (mV)')
f.line(t, soma_v, line_width=2)
plt.show(f)