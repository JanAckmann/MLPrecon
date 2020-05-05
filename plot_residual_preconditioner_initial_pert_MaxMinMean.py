# -*- coding: utf-8 -*-
#module load python/2.7-ve3


import os as os
import os.path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import matplotlib
## data_NOADI_NOTgcr_EXP3_Dp_1M10_pert_res05
## data_ADI_NOTgcr_EXP3_Dp_1M10_pseut_res05
## data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_res05
## data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_pert_res05
## data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_Approx10_res05
## data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_Initial_res05

## data_ADIprecon_NOTgcr_EXP3_Dp_1M10_res05
## data_NOprecon_NOTgcr_EXP3_Dp_1M10_res05
## data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_Approx10_pert_res05

#path='../data/data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_FULL_res05/'
#path='../data/data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_INIT_res05/'
#path='../data/data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_REDU_res05/'
#path='../data/data_ADIprecon_NOTgcr_EXP3_Dp_1M10_res05/' 
#Precon='5'

path='../data/data_NOprecon_NOTgcr_EXP3_Dp_1M10_res05/'
Precon='0'

#data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_Approx10_res05/'#'../data/data_ADI_NOTgcr_EXITcond2M4_Dp_res2/'# '../data/data_ADI_NOTgcr_EXITcond5M4_Dp_res2/' #'../data/data_ADI_NOTgcr_D4_2Mplot_residual_preconditioner_initial_pert.py4_RP_res2/'
path_ref= 'data/'
plotpath= path

varlist=['R']
exp='3'
#Precon='0'
codesD='F'
codesQ='F'

solver='implicit_HP_implPD'

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 13}

matplotlib.rc('font', **font)

start_time= 360
end_time=  2856#20*240
jump=1 #24#*5

# initial plot
#start_time= 0
#end_time=   360#20*240
#jump=1



import numpy as np


def calc_average(in_file, start, end, jump):
  file_timeslice= in_file+str(start)
  counter=1
  x,y,average = np.loadtxt(file_timeslice, usecols=(0,1,2), unpack=True)

  for timestep in range(start,end+1, jump):
    print( timestep)
    file_timeslice= filename+str(timestep)
    counter+=1
    x,y,temp = np.loadtxt(file_timeslice, usecols=(0,1,2), unpack=True)
    average+= temp
 
  return x, y, average/float(counter)
 
def calc_std(in_file, average, start, end, jump):
  counter=0
  std=np.zeros(len(average), dtype=float)

  for timestep in range(start,end+1, jump):
    file_timeslice= filename+str(timestep)
    counter+=1
    x,y,temp = np.loadtxt(file_timeslice, usecols=(0,1,2), unpack=True)
    std+= (temp-average)**2.0
 
  return x, y, np.sqrt(std/float(counter-1))

 
 
############## MAIN PROGRAMM ######################  
number=(end_time-start_time)/jump
cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in np.linspace(1, 0, int(number)+1)]




for var in varlist:
  for bits in range(52,51,-2):
    plt.subplot(2, 1, 1)
    errors= np.zeros((3, 17))
    count=0
    max_iteration=0
    #print((end_time-start_time)/jump)
    residuum_infnorm= np.zeros((int((end_time-start_time)/jump)+1,100))
    for time in range(start_time, end_time+1,jump):
      plotted= False
      #print(time)
      residuum_norm= np.zeros((100))

      exit= np.zeros((100))
      for iteration in range(0,100,1):

        filenameT = path+'Precon'+Precon+'_' + var+'_exp'+exp+'_time'+str(time)+'_iter_'+str(iteration)+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'.txt'
        filenameF = path_ref + var+'_exp'+exp+'_time'+str(time)+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'.txt'
        #file_refT = path+ var+'_exp'+exp+'_time'+str(time)+'_codes_T'+'_bits52.txt'
        file_refF = path_ref+ var+'_exp'+exp+'_time'+str(time)+'_codes_'+str(codesD)+str(codesQ)+'_bits52.txt'
        #print 'calc_average'
        #print(filenameT)

        #xcoord, ycoord, mean_var= calc_average(filename, start_time, end_time, jump)
        #print 'done'
        #print 'savetxt mean'
        #np.savetxt( filename+'avg'+str(start_time)+'to'+str(end_time), np.transpose([xcoord,ycoord,mean_var]))

 

        #print 'calc_standard dev'


        #xcoord, ycoord, std_var= calc_std(filename, mean_var, start_time, end_time, 5000)
        #print 'done'
        #print 'savetxt std'
        #np.savetxt( filename+'std'+str(start_time)+'to'+str(end_time), np.transpose([xcoord,ycoord,std_var]))


        ### plot std and mean


        #print( 'load residual')
        #xcoord, ycoord, mean_var_refF =np.loadtxt( file_refF, usecols=(0,1,2), unpack=True)
        #xcoord, ycoord, mean_var_refT =np.loadtxt( file_refT, usecols=(0,1,2), unpack=True)
        if (os.path.exists(filenameT)):
          xcoord, ycoord, Residuum, exitcon   =np.loadtxt( filenameT, usecols=(0,1,2,3), unpack=True)
          #xcoord, ycoord, mean_varF     =np.loadtxt( filenameF, usecols=(0,1,2), unpack=True)
          ncols, nrows = len(set(xcoord)), len(set(ycoord))
          grid = np.flipud(Residuum.reshape((nrows, ncols), order='F'))
          residuum_norm[iteration]=np.linalg.norm (Residuum)
          #print((time-start_time)/jump)
          residuum_infnorm[int((time-start_time)/jump),iteration]=np.linalg.norm (Residuum , np.inf)
          #if iteration ==0:
            #ind=np.argmax(grid[0,:])
            #print(np.max(grid[0,:]), np.max(grid[1,:]), np.max(grid[2,:]), np.max(grid[3,:]), np.max(grid[4,:]), np.max(grid[5,:]))
          #print(residuum_infnorm[iteration])
          exit[iteration]=np.linalg.norm (exitcon , np.inf)
          #print(exit[iteration])
          max_iteration=max(iteration,max_iteration)
          #print(max_iteration)
        else:
          residuum_infnorm[int((time-start_time)/jump),iteration]=-1.0

    number_of_iters=0
    for iteration in range(0,100,1):
      if(max(residuum_infnorm[:,iteration]) != min(residuum_infnorm[:,iteration])):
        #print(residuum_infnorm[:,iteration])
        #print(max(residuum_infnorm[:,iteration]))
        #print(min(residuum_infnorm[:,iteration]))
        number_of_iters=iteration+1
    print(number_of_iters)

    for time in range(int((0)/jump), int((end_time-start_time)/jump)+1,int(jump/jump)):
      #print((residuum_infnorm[time,0])*(5.0*10**(-9)))
      residuum_infnorm[time,:]=residuum_infnorm[time,:]/(residuum_infnorm[time,0]*(1.0*10**(-10)))
      #print(time)
      #print(residuum_infnorm[time,:number_of_iters])

    #residuum_infnorm[:,:]=residuum_infnorm[:,:]*(5.0*10**(9))
    Minresiduum_infnorm= np.zeros((number_of_iters))
    Maxresiduum_infnorm= np.zeros((number_of_iters))
    Meanresiduum_infnorm= np.zeros((number_of_iters))

    for iteration in range(0,number_of_iters,1):
      Maxresiduum_infnorm[iteration]=max(residuum_infnorm[:,iteration])
    print('max')
    print(Maxresiduum_infnorm)
    for iteration in range(0,number_of_iters,1):
       Minresiduum_infnorm[iteration]=min(residuum_infnorm[:,iteration])
       if(Minresiduum_infnorm[iteration] >0.0):
         number_of_iters_min=iteration+1

    print('min')
    print(Minresiduum_infnorm)
    for iteration in range(0,number_of_iters,1):
       Meanresiduum_infnorm[iteration]=np.median(residuum_infnorm[:,iteration])
       if(Meanresiduum_infnorm[iteration] >0.0):
         number_of_iters_mean=iteration+1     
    print('mean')
    print(Meanresiduum_infnorm)       
       
    number_of_iters=np.arange(0,number_of_iters,1)#np.zeros((iteration))
    plt.plot(number_of_iters, Maxresiduum_infnorm[:], 'g--', label ='Max') #/exitcon[0]

    number_of_iters=np.arange(0,number_of_iters_mean,1)#np.zeros((iteration))
    plt.plot(number_of_iters, Meanresiduum_infnorm[:number_of_iters_mean],'r-', label ='Median') #/exitcon[0]

    number_of_iters=np.arange(0,number_of_iters_min,1)#np.zeros((iteration))
    plt.plot(number_of_iters, Minresiduum_infnorm[:number_of_iters_min],'b--', label ='Min') #/exitcon[0]

    plt.xlabel('Solver Iterations')
    plt.ylabel('Linf Residual')
    #number_of_iters=np.arange(0,10,1)
    plt.xticks(np.arange(0,14,1))
    #plt.plot((0.0,float(max_iteration) ), (1.0,1.0), 'k-')
    plt.plot((0.0,float(13) ), (1.0,1.0), 'k-')
    #plt.legend()
    #plt.legend(loc='upper left', prop={'size':8}, ncol=2, bbox_to_anchor=(0.78,1.37), framealpha=1)
    #plt.legend(loc='upper left', prop={'size':7}, ncol=2, bbox_to_anchor=(0.80,1), framealpha=1)
    #plt.legend(loc='upper left', prop={'size':7}, bbox_to_anchor=(1,1)) #(loc='upper right')  #['day0'], loc='upper right', 'day1','day2', 'day3','day4', 'day5','day6', 'day7', 'day8', 'day9','day10', 'day11','day12', 'day13','day14', 'day15', 'day16'
    plt.yscale('log')
    plt.ylim(ymin = 10**(-3))
    plt.axis
    #plt.show()
    #plt.savefig(plotpath +'L2_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_FF'+'_bits'+str(bits)+'.pdf', bbox_inches=0)
    #plt.close()
 
    #print( plotpath +'L2_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_FF'+'_bits'+str(bits)+'.pdf')

    plt.savefig(plotpath +'Norm_Errors_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'_'+str(end_time)+'unnorm_WoLegend.png')
    plt.close()
 
    print( plotpath +'Norm_Errors_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'_'+str(end_time)+'unnorm_WoLegend.png')


#compare_histograms_normal('h', 31, sigma_h) #, -0.0526478751975,0.0526478751975 )
#compare_histograms_talone('t', 31, sigma_t, -0.0280507825342,0.0280507825342 )
#compare_histograms('vn', 31, sigma_vn,-0.012620211929 , 0.012620211929)

#call_truncation_errors('h',n_cells, 1, num_timesteps)
#call_truncation_errors('t',n_cells, 32, num_timesteps)
#call_truncation_errors('vn',n_edges, 32, num_timesteps)

