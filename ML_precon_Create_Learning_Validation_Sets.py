import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

path='../data/data_ADI_NOTgcr_EXP3_Dp_1M10_res05'

grid_zonal=64
grid_merid=32

zonal_ext=2
zonal_ext_o=0
merid_ext=2

start_learning= 5000 # 5000
end_learning=  22000 #10000

start_valid= 22001
end_valid=  27000

learning_set=[]
for i in range(5040,10080,1):
  learning_set.append(i)
for i in range(12600,17640,1):
  learning_set.append(i)
for i in range(20160,25200,1):
  learning_set.append(i)
for i in range(27720,32760,1):
  learning_set.append(i)
for i in range(35280,40320,1):
  learning_set.append(i)

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

# 43200
def save_pair(x,y,fname):
    import h5py
    hf = h5py.File(fname, 'w')
    hf.create_dataset('x', data=x, dtype='float64')
    hf.create_dataset('y', data=y, dtype='float64')
    hf.close()
    return

#inputs=max(zonal_ext,merid_ext*2+1)*3+ ((merid_ext-1)*2+1)*2+ ((merid_ext-2)*2+1)*2  #zonal bands for 3 inputs+stencil
inputs=1+(zonal_ext*2+1)*7*(merid_ext*2+1)
outputs=zonal_ext_o*2+1 #1
print(inputs, outputs)

x=np.zeros(((len(learning_set))*grid_zonal,inputs))
y=np.zeros(((len(learning_set))*grid_zonal,outputs))

print ((len(learning_set))*grid_zonal)

x_val=np.zeros(((len(validation_set))*grid_zonal,inputs))
y_val=np.zeros(((len(validation_set))*grid_zonal,outputs))

print ((len(validation_set))*grid_zonal)
latitudes=np.zeros(32)

file_grid = path+'/Precon5_H_exp3_time0_codes_FT_bits52.txt'
xcoord, ycoord=np.loadtxt( file_grid, usecols=(0,1), unpack=True)

ncols, nrows = len(set(xcoord)), len(set(ycoord)) 
latitudes= sorted(set(ycoord))


file_H0_in = path+'/Topo/'+'Timestep'+str(0)+'.txt'
H0_in=np.loadtxt( file_H0_in, usecols=(0), unpack=True)
grid_H0_in = np.flipud(H0_in.reshape((nrows, ncols), order='F'))

for lat in range(0,nrows):


  indexing=0
  for time in learning_set: #range(start_learning,end_learning+1):      

    file_h_in = path+'/H_in/'+'Timestep'+str(time)+'.txt'
    file_h_out = path+'/H_in/'+'Timestep'+str(time+1)+'.txt'

    file_R_iter = path+'/R_iter/'+'Timestep'+str(time)+'_iter0.txt'
    file_a11 = path+'/a11/'+'Timestep'+str(time)+'.txt'
    file_a12 = path+'/a12/'+'Timestep'+str(time)+'.txt'
    file_a21 = path+'/a21/'+'Timestep'+str(time)+'.txt'
    file_a22 = path+'/a22/'+'Timestep'+str(time)+'.txt'
    file_b11 = path+'/b11/'+'Timestep'+str(time)+'.txt'
    file_b22 = path+'/b22/'+'Timestep'+str(time)+'.txt'

    h_in=np.loadtxt( file_h_in, usecols=(0), unpack=True)
    h_out=np.loadtxt( file_h_out, usecols=(0), unpack=True)

    R_iter_in=np.loadtxt( file_R_iter, usecols=(0), unpack=True)
    a11_in=np.loadtxt( file_a11, usecols=(0), unpack=True)
    a12_in=np.loadtxt( file_a12, usecols=(0), unpack=True)
    a21_in=np.loadtxt( file_a21, usecols=(0), unpack=True)
    a22_in=np.loadtxt( file_a22, usecols=(0), unpack=True)
    b11_in=np.loadtxt( file_b11, usecols=(0), unpack=True)
    b22_in=np.loadtxt( file_b22, usecols=(0), unpack=True)

    grid_h_in = np.flipud(h_in.reshape((nrows, ncols), order='F'))
    grid_h_out = np.flipud(h_out.reshape((nrows, ncols), order='F'))

    grid_R_iter_in = np.flipud(R_iter_in.reshape((nrows, ncols), order='F'))
    grid_a11_in = np.flipud(a11_in.reshape((nrows, ncols), order='F'))
    grid_a12_in = np.flipud(a12_in.reshape((nrows, ncols), order='F'))
    grid_a21_in = np.flipud(a21_in.reshape((nrows, ncols), order='F'))
    grid_a22_in = np.flipud(a22_in.reshape((nrows, ncols), order='F'))
    grid_b11_in = np.flipud(b11_in.reshape((nrows, ncols), order='F'))
    grid_b22_in = np.flipud(b22_in.reshape((nrows, ncols), order='F'))

    grid_h_in = grid_h_in /9.80616
    grid_h_out=grid_h_out /9.80616
    grid_R_iter_in = grid_R_iter_in /9.80616

    grid_h_out= grid_h_in-grid_h_out   # I want to predict the error

    #plt.imshow(grid_R_in, extent=(xcoord.min(), xcoord.max(), ycoord.min(), ycoord.max()),interpolation='nearest', cmap=cm.bwr)
    #plt.show()

    #for lat in range(0,nrows):
    for lon in range(0,ncols):
      
      x[indexing,0]=latitudes[lat]
      counter=1

      counter_out=0
      for lon_box in range(-zonal_ext_o,zonal_ext_o+1):
         pos_lon=(lon+lon_box)%ncols
         if (pos_lon<0):
           pos_lon=ncols-pos_lon+1
         y[indexing,counter_out]=grid_h_out[lat, pos_lon]
         #print(lat, pos_lon, counter_out)
         counter_out+=1
         
      for lat_box in range(-merid_ext,merid_ext+1):
        for lon_box in range(-zonal_ext,zonal_ext+1):
           pos_lat= lat+lat_box
           pos_lon=(lon+lon_box)%ncols
           if (pos_lat<0):
             pos_lat=abs(pos_lat)-1
             pos_lon=(lon+int(ncols/2)+lon_box)%ncols
           elif (pos_lat>=nrows):
             pos_lat=pos_lat-(pos_lat%nrows+1)
             pos_lon=(lon+int(ncols/2)+lon_box)%ncols
           if (pos_lon<0):
             pos_lon=ncols-pos_lon+1
           
           x[indexing,counter]=  grid_R_iter_in[pos_lat, pos_lon]
           x[indexing,counter+1]=grid_a11_in[pos_lat, pos_lon]
           x[indexing,counter+2]=grid_a12_in[pos_lat, pos_lon]
           x[indexing,counter+3]=grid_a21_in[pos_lat, pos_lon]
           x[indexing,counter+4]=grid_a22_in[pos_lat, pos_lon]
           x[indexing,counter+5]=grid_b11_in[pos_lat, pos_lon]
           x[indexing,counter+6]=grid_b22_in[pos_lat, pos_lon]

           counter=counter+7  # continue counting
      indexing+=1
           #print(counter)
  print(indexing)
  print ('writing file')
  save_pair(x,y,'../data/ML_data_ADI_NOTgcr_EXP3_Dp_1M10_res05/training_Lcoeff_R_'+str(zonal_ext)+'x'+str(merid_ext)+'_1x0'+str(lat)+'DP.h5')


  indexing=0
  for time in validation_set: #range(start_valid, end_valid+1):
    file_h_in = path+'/H_in/'+'Timestep'+str(time)+'.txt'
    file_h_out = path+'/H_in/'+'Timestep'+str(time+1)+'.txt'

    file_R_iter = path+'/R_iter/'+'Timestep'+str(time)+'_iter0.txt'
    file_a11 = path+'/a11/'+'Timestep'+str(time)+'.txt'
    file_a12 = path+'/a12/'+'Timestep'+str(time)+'.txt'
    file_a21 = path+'/a21/'+'Timestep'+str(time)+'.txt'
    file_a22 = path+'/a22/'+'Timestep'+str(time)+'.txt'
    file_b11 = path+'/b11/'+'Timestep'+str(time)+'.txt'
    file_b22 = path+'/b22/'+'Timestep'+str(time)+'.txt'

    h_in=np.loadtxt( file_h_in, usecols=(0), unpack=True)
    h_out=np.loadtxt( file_h_out, usecols=(0), unpack=True)

    R_iter_in=np.loadtxt( file_R_iter, usecols=(0), unpack=True)
    a11_in=np.loadtxt( file_a11, usecols=(0), unpack=True)
    a12_in=np.loadtxt( file_a12, usecols=(0), unpack=True)
    a21_in=np.loadtxt( file_a21, usecols=(0), unpack=True)
    a22_in=np.loadtxt( file_a22, usecols=(0), unpack=True)
    b11_in=np.loadtxt( file_b11, usecols=(0), unpack=True)
    b22_in=np.loadtxt( file_b22, usecols=(0), unpack=True)

    grid_h_in = np.flipud(h_in.reshape((nrows, ncols), order='F'))
    grid_h_out = np.flipud(h_out.reshape((nrows, ncols), order='F'))

    grid_R_iter_in = np.flipud(R_iter_in.reshape((nrows, ncols), order='F'))
    grid_a11_in = np.flipud(a11_in.reshape((nrows, ncols), order='F'))
    grid_a12_in = np.flipud(a12_in.reshape((nrows, ncols), order='F'))
    grid_a21_in = np.flipud(a21_in.reshape((nrows, ncols), order='F'))
    grid_a22_in = np.flipud(a22_in.reshape((nrows, ncols), order='F'))
    grid_b11_in = np.flipud(b11_in.reshape((nrows, ncols), order='F'))
    grid_b22_in = np.flipud(b22_in.reshape((nrows, ncols), order='F'))

    grid_h_in = grid_h_in /9.80616
    grid_h_out=grid_h_out /9.80616
    grid_R_iter_in = grid_R_iter_in /9.80616

    grid_h_out= grid_h_in-grid_h_out   # I want to predict the error

    #plt.imshow(grid_R_in, extent=(xcoord.min(), xcoord.max(), ycoord.min(), ycoord.max()),interpolation='nearest', cmap=cm.bwr)
    #plt.show()

    #for lat in range(0,nrows):
    for lon in range(0,ncols):
      
      x_val[indexing,0]=latitudes[lat]
      counter=1

      counter_out=0
      for lon_box in range(-zonal_ext_o,zonal_ext_o+1):
         pos_lon=(lon+lon_box)%ncols
         if (pos_lon<0):
           pos_lon=ncols-pos_lon+1
         y_val[indexing,counter_out]=grid_h_out[lat, pos_lon]
         #print(lat, pos_lon, counter_out)
         counter_out+=1
         
      for lat_box in range(-merid_ext,merid_ext+1):
        for lon_box in range(-zonal_ext,zonal_ext+1):
           pos_lat= lat+lat_box
           pos_lon=(lon+lon_box)%ncols
           if (pos_lat<0):
             pos_lat=abs(pos_lat)-1
             pos_lon=(lon+int(ncols/2)+lon_box)%ncols
           elif (pos_lat>=nrows):
             pos_lat=pos_lat-(pos_lat%nrows+1)
             pos_lon=(lon+int(ncols/2)+lon_box)%ncols
           if (pos_lon<0):
             pos_lon=ncols-pos_lon+1
           
           x_val[indexing,counter]=  grid_R_iter_in[pos_lat, pos_lon]
           x_val[indexing,counter+1]=grid_a11_in[pos_lat, pos_lon]
           x_val[indexing,counter+2]=grid_a12_in[pos_lat, pos_lon]
           x_val[indexing,counter+3]=grid_a21_in[pos_lat, pos_lon]
           x_val[indexing,counter+4]=grid_a22_in[pos_lat, pos_lon]
           x_val[indexing,counter+5]=grid_b11_in[pos_lat, pos_lon]
           x_val[indexing,counter+6]=grid_b22_in[pos_lat, pos_lon]

           counter=counter+7  # continue counting
      indexing+=1
           #print(counter)
  print(indexing)
  print ('writing file')
  save_pair(x_val,y_val,'../data/ML_data_ADI_NOTgcr_EXP3_Dp_1M10_res05/validation_Lcoeff_R_'+str(zonal_ext)+'x'+str(merid_ext)+'_1x0'+str(lat)+'DP.h5')


def load_pair(fname):
    import h5py
    hf = h5py.File(fname, 'r')
    x=hf['x'].value
    y=hf['y'].value
    hf.close()
    return x,y
