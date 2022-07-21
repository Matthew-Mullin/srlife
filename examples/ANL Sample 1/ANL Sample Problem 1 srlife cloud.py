#%%
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 10:58:36 2022

@author: matthew
"""


import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../..')

from srlife import receiver

def pressure(pmax,t,a,b,day): #generates an array of dimension: time
    aArray = (pmax/a)*times[np.where(times[:] < a)]
    bArray = pmax+(-pmax/(day-b))*(times[np.where((times[:] > b) & (times[:] < day))]-b)
    dayArray = pmax*np.ones(len(np.where((times[:] <= b) & (times[:] >= a))[0]))
    nightArray = np.zeros(len(np.where(times[:] >= day)[0]))
    pArray = np.concatenate([aArray, dayArray, bArray, nightArray], axis=None)
    return pArray
    
def innerTemp(c, d, day, times, T_in_bot, T_in_top, T_base, nt, nz): #generates an array of dimension: time x theta x z
    #generate Time Variation
    cArray = (1/c)*times[np.where(times[:] < c)]
    dArray = 1+(-1/(day-d))*(times[np.where((times[:] > d) & (times[:] <= day))]-d)
    dayArray = 1*np.ones(len(np.where((times[:] <= d) & (times[:] >= c))[0]))
    nightArray = np.zeros(len(np.where(times[:] > day)[0]))
    tArray = np.concatenate([cArray, dayArray, dArray, nightArray], axis=None)  
    
    T_in_z = np.linspace((T_in_bot-T_base),(T_in_top-T_base), nz)
    thetaArray = np.ones(nt)
    innerArray = (tArray[:, None, None]*thetaArray[None,:,None]*T_in_z[None, None, :]+T_base) 
    return innerArray

def outerTemp(a, b, c, d, day, times, T_in_bot, T_in_top, T_out_bot, T_out_top, T_base, innerTempArray, nt, nz): #generates an array of dimension: time x theta x z
    cArray = np.zeros(len(np.where(times[:] < c)[0]))    
    aArray = ((1-0)/(a-c))*(times[np.where((times[:] <= a) & (times[:] >= c))]-c)+0
    dayArray = 1*np.ones(len(np.where((times[:] < b) & (times[:] > a))[0]))
    bArray = ((0-1)/(d-b))*(times[np.where((times[:] <= d) & (times[:] >= b))]-b)+1
    dArray = np.zeros(len(np.where((times[:] > d) & (times[:] < day))[0]))  
    nightArray = np.zeros(len(np.where(times[:] >= day)[0]))
    tArray = np.concatenate([cArray, aArray, dayArray, bArray, dArray, nightArray], axis=None)
    
    T_out_z = np.linspace((T_out_bot-T_in_bot),(T_out_top-T_in_top), nz)
    thetaArray = np.ones(nt)
    outerArray = (tArray[:, None, None]*thetaArray[None,:,None]*T_out_z[None, None, :]) 
    
    returnArray = innerTempArray+outerArray    
    return returnArray


if __name__ == "__main__":
  # Setup the base receiver
  period = 24.0 # Loading cycle period, hours
  days = 1 # Number of cycles represented in the problem 
  panel_stiffness = "disconnect" # Panels are disconnected from one another

  model = receiver.Receiver(period, days, panel_stiffness)

  # Setup each of the two panels
  tube_stiffness = "disconnect"
  panel_0 = receiver.Panel(tube_stiffness)

  # Basic receiver geometry
  r_outer = 20 # mm
  thickness = 2.0 # mm
  height = 500.0 # mm

  # Tube discretization
  nr = 6#12
  nt = 20#40
  nz = 15#30

  # Inner/Outer Surface Temperature Parameters
  delta_T = 100 #C #how does the inner and outer temps differ with height #Estimated
  #T_out_top= 770 #C #outer at max z
  T_out_top = 825 #C
  T_out_bot = T_out_top-delta_T #C #Estimated
  T_in_bot = 640 #C  #inner at min z
  T_in_top = T_in_bot + delta_T #C #Estimated
  
  #startup/shutdown hours
  day = 12 #hours in day
  a = 2 #hours
  b = 10 #hours
  c = 1 #hours
  d = 11 #hours
  
  cloudTimeStart=np.array((4,5,6,7,8))
  offset = 1/60 #hours
  cloudTimeBefore = cloudTimeStart-offset
  cloudTimeStop=cloudTimeStart+(8/60)
  cloudTimeAfter = cloudTimeStop+offset
  
  evenlySpaced = False
  if evenlySpaced:   
      # Time increments throughout the 24 hour day
      #min increment is 8 mins
      times=np.linspace(0,24,360+1)
  else:
      times = np.sort(np.concatenate((np.asarray((0, c, a, b, d, day, 24)),cloudTimeStart,cloudTimeStop)))
      res = 3
      resizedTimes = np.zeros((len(times)-1)*res+1)
      for i in range(len(times)-1):
          resizedTimes[i*res:(i+1)*res] = np.linspace(times[i], times[i+1], res, endpoint=False)
      resizedTimes[-1] = times[-1]
      times = resizedTimes
  times = np.sort(np.concatenate((times,cloudTimeBefore,cloudTimeAfter)))
   
  T_base = 30 # C
  innerTempArray = innerTemp(c, d, day, times, T_in_bot, T_in_top, T_base, nt, nz)
  outerTempArray = outerTemp(a, b, c, d, day, times, T_in_bot, T_in_top, T_out_bot, T_out_top, T_base, innerTempArray, nt, nz)
  
  # ID pressure history
  p_max = 2.0 # MPa
  pressureArray = pressure(p_max, times, a, b, day)
  
  #cloud Events
  cloudPressure = 0.5 #MPa
  cloudTemp = 550 #C #T_base
  #cloudTime=np.array((4,5,6,7,8));
  
  index =  np.zeros(len(times),dtype=bool)
  for i in range(len(cloudTimeStart)):
      indexTemp = (times >= cloudTimeStart[i]) & (times<cloudTimeStop[i])
      index[indexTemp] = True
  innerTempArray[index,:,:] = cloudTemp;
  outerTempArray[index,:,:] = cloudTemp;    
  pressureArray[index] = cloudPressure;
  
  #plotting
  fig, ax1 = plt.subplots()

  color = 'red'
  ax1.set_xlabel('Time [hr]')
  ax1.set_ylabel('Temperature [C]', color=color)
  ax1.plot(times, innerTempArray[:,0,0], color='blue', label='Tinner(t,z=0mm)')
  ax1.plot(times, innerTempArray[:,0,-1], color='green', label='Tinner(t,z=500mm)')
  ax1.plot(times, outerTempArray[:,0,0], color='red', label='Touter(t,z=0mm)')
  ax1.plot(times, outerTempArray[:,0,-1], color='yellow', label='Touter(t,z=500mm)')
  ax1.tick_params(axis='y', labelcolor=color)
  ax1.legend()

  ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

  color = 'black'
  ax2.set_ylabel('Pressure (MPa)', color=color)  # we already handled the x-label with ax1
  ax2.plot(times, pressureArray, color=color, label='pressure')
  ax2.tick_params(axis='y', labelcolor=color)
  ax2.set_ylim(0, 6)
  ax2.legend(loc='right')

  fig.tight_layout()  # otherwise the right y-label is slightly clipped
  plt.show()
  
  #Convert all temps to Kelvin
  innerTempArray += 273.15 #K
  outerTempArray += 273.15 #K
  T_base += 273.15 #K

  # Setup each tube in turn and assign it to the correct panel
  # Tube 0
  tube_0 = receiver.Tube(r_outer, thickness, height, nr, nt, nz, T0 = T_base)
  tube_0.set_times(times)
  
  tube_0.set_bc(receiver.FixedTempBC(r_outer-thickness,
    height, nt, nz, times, innerTempArray), "inner")
  
  tube_0.set_bc(receiver.FixedTempBC(r_outer, height,
    nt, nz, times, outerTempArray), "outer")
  
  tube_0.set_pressure_bc(receiver.PressureBC(times, pressureArray))

  # Assign to panel 0
  panel_0.add_tube(tube_0)
  
  model.add_panel(panel_0)

  # Save the receiver to an HDF5 file
  model.save("sample1_Cloud.hdf5")
# %%
