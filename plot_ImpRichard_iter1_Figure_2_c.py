import os as os
import os.path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plotpath='Model_plots/'
path='../data/data_ADI_NOTgcr_EXP3_Dp_1M10_pseut_res05'

grid_zonal=64
grid_merid=32

zonal_ext=2
zonal_ext_o=0
merid_ext=2


validation_set=[10440]
for i in range(10440,12240,1):
  validation_set.append(i)
for i in range(18000,19800,1):
  validation_set.append(i)
for i in range(25560,27360,1):
  validation_set.append(i)
for i in range(33120,34920,1):
  validation_set.append(i)
for i in range(40680,42839,1):
  validation_set.append(i)



#inputs=max(zonal_ext,merid_ext*2+1)*3+ ((merid_ext-1)*2+1)*2+ ((merid_ext-2)*2+1)*2  #zonal bands for 3 inputs+stencil

latitudes=np.zeros(32)



file_grid = path+'/Precon5_H_exp3_time0_codes_FT_bits52.txt'
xcoord, ycoord=np.loadtxt( file_grid, usecols=(0,1), unpack=True)

ncols, nrows = len(set(xcoord)), len(set(ycoord)) 
latitudes= sorted(set(ycoord))


file_H0_in = path+'/Topo/'+'Timestep'+str(0)+'.txt'
H0_in=np.loadtxt( file_H0_in, usecols=(0), unpack=True)
grid_H0_in = np.flipud(H0_in.reshape((nrows, ncols), order='F'))

orig_error=np.zeros(((len(validation_set))*grid_zonal))
iter_error=np.zeros(((len(validation_set))*grid_zonal))




for lat in range(31,nrows):
  print(lat)
  orig_error[:]=0.0
  timecounter=0
  for time in validation_set:

    file_h_in = path+'/H_in/'+'Timestep'+str(time)+'.txt'
    file_h_out = path+'/H_in/'+'Timestep'+str(time+1)+'.txt'


    h_in=np.loadtxt( file_h_in, usecols=(0), unpack=True)
    h_out=np.loadtxt( file_h_out, usecols=(0), unpack=True)

    grid_h_in = np.flipud(h_in.reshape((nrows, ncols), order='F'))
    grid_h_out = np.flipud(h_out.reshape((nrows, ncols), order='F'))

    grid_h_in = grid_h_in /9.80616
    grid_h_out=grid_h_out /9.80616


    grid_h_in= grid_h_in-grid_h_out   # I want to predict the error

    for lon in range(0,ncols):
      orig_error[(timecounter)*grid_zonal+lon]=grid_h_in[lat,lon]
    timecounter+=1

  for  iteration in range(1,2):
   print(iteration)
   iter_error[:]=0.0
   timecounter=0
   for time in validation_set:
    print(time)

    file_h_iter = path+'/H_iter/'+'Timestep'+str(time)+'_iter'+str(iteration)+'.txt'
    if (os.path.exists(file_h_iter)):
      file_h_out = path+'/H_in/'+'Timestep'+str(time+1)+'.txt'

      h_iter=np.loadtxt( file_h_iter, usecols=(0), unpack=True)
      h_out=np.loadtxt( file_h_out, usecols=(0), unpack=True)

      grid_h_iter = np.flipud(h_iter.reshape((nrows, ncols), order='F'))
      grid_h_out = np.flipud(h_out.reshape((nrows, ncols), order='F'))

      grid_h_iter = grid_h_iter /9.80616
      grid_h_out=grid_h_out /9.80616

      grid_h_iter= grid_h_iter-grid_h_out
  
      for lon in range(0,ncols):
        iter_error[(timecounter)*grid_zonal+lon]=grid_h_iter[lat,lon]
        #print(iter_error[(timecounter)*grid_zonal+lon]/orig_error[(timecounter)*grid_zonal+lon])
        #raw_input('lala')
      timecounter+=1
    else:
      for lon in range(0,ncols):
        iter_error[(timecounter)*grid_zonal+lon]=0.0
        #print(iter_error[(timecounter)*grid_zonal+lon]/orig_error[(timecounter)*grid_zonal+lon])
        #raw_input('lala')
      timecounter+=1


   Relative_dec=np.zeros((len(iter_error))) 
   Relative_dec[:]= np.absolute((iter_error[:])/orig_error[:])
   plt.xlabel(r'$\Delta \Phi$', fontsize=16)
   plt.ylabel("Relative Decrease in Absolute Error", fontsize=16)
   plt.ylim(10.0**(-9.0), 10.0**(4.0))
   #Relative_dec[:]= np.absolute((y_ML_model[:,0]-y_val[:,0])/y_val[:,0])
   plt.ylim(10.0**(-9.0), 10.0**(4.0))
   plt.yscale('log')
   plt.ylim(10.0**(-9.0), 10.0**(4.0))
   plt.axhline(y=10.0**(-1.0), color='r', linestyle='--')
   plt.scatter(orig_error, Relative_dec, vmin=10.0**(-9.0), vmax=10.0**(4.0))
   plt.ylim(10.0**(-9.0), 10.0**(4.0))
   plt.savefig(plotpath +'Scatter_ADI_PRECON_iteration'+str(iteration)+'lat_'+str(lat)+'_pseut_ext.jpg', bbox_inches=0)
   plt.close()



  
