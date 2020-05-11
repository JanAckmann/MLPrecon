import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

path='../data/data_ADI_NOTgcr_EXP3_Dp_1M10_pseut_res05'

grid_zonal=64
grid_merid=32

zonal_ext=2
zonal_ext_o=0
merid_ext=2


validation_set=[]
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

iteration=1

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




for lat in range(0,nrows):
  print(lat)
  orig_error[:]=0.0
  iter_error[:]=0.0
  timecounter=0
  for time in validation_set:

    file_h_in = path+'/H_in/'+'Timestep'+str(time)+'.txt'
    file_R_in = path+'/RHS/'+'Timestep'+str(time)+'.txt'
    file_h_out = path+'/H_in/'+'Timestep'+str(time+1)+'.txt'
    file_h_iter = path+'/H_iter/'+'Timestep'+str(time)+'_iter'+str(iteration)+'.txt'

    h_in=np.loadtxt( file_h_in, usecols=(0), unpack=True)
    R_in=np.loadtxt( file_R_in, usecols=(0), unpack=True)
    h_out=np.loadtxt( file_h_out, usecols=(0), unpack=True)
    h_iter=np.loadtxt( file_h_iter, usecols=(0), unpack=True)

    grid_h_in = np.flipud(h_in.reshape((nrows, ncols), order='F'))
    grid_R_in = np.flipud(R_in.reshape((nrows, ncols), order='F'))
    grid_h_iter = np.flipud(h_iter.reshape((nrows, ncols), order='F'))

    grid_h_out = np.flipud(h_out.reshape((nrows, ncols), order='F'))

    grid_h_in = grid_h_in /9.80616

    grid_h_out=grid_h_out /9.80616
    grid_R_in = grid_R_in /9.80616

    grid_h_iter = grid_h_iter /9.80616
    #print(time)
    #print(grid_h_iter)
    #print(grid_h_in)
    #print(grid_h_out)

    grid_h_in= grid_h_in-grid_h_out   # I want to predict the error
    grid_h_iter= grid_h_iter-grid_h_out

    grid_R_in= grid_h_in-grid_R_in 
  
    for lon in range(0,ncols):
      orig_error[(timecounter)*grid_zonal+lon]=grid_h_in[lat,lon]
      iter_error[(timecounter)*grid_zonal+lon]=grid_h_iter[lat,lon]

    timecounter+=1

  print(str(lat))
  print(np.mean(np.absolute(orig_error)) )
  print(np.mean(np.absolute(iter_error)) )
  f = open("truth_diff_METRICnorm_iter_1_pseut_Lats.txt", "a")
  f.write(str(lat)+' '+str(np.mean(np.absolute(orig_error))) +' '+ str(np.mean(np.absolute(orig_error)))  +' '+str(np.mean(np.absolute(iter_error)) ) + '\n')
  f.close()


  
