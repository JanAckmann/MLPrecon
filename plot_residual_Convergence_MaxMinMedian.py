import os as os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import matplotlib


#path='../data/data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_FULL_res05/'  # full input convergence; perturbed run
#path='../data/data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_INIT_res05/'  # initial phase convergence; unperturbed run
#path='../data/data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_REDU_res05/'  # reduced input convergence; perturbed run
#path='../data/data_ADIprecon_NOTgcr_EXP3_Dp_1M10_res05/'                  # implicit Richardson convergence; unperturbed run
#Precon='5'                                                                # '5' means preconditioning is active -- needs to be set for all 4 possible paths above

path='../data/data_NOprecon_NOTgcr_EXP3_Dp_1M10_res05/'                    # no preconditioner; unperturbed run
Precon='0'                                                                 # '0' means no preconditioner is active -- use only in combination with "path='../data/data_NOprecon_NOTgcr_EXP3_Dp_1M10_res05/'"

plotpath= path                                                             # plot in data directory


##  other experiment identifiers, keep unchanged
varlist=['R']
exp='3'
codesD='F'
codesQ='F'


## fonts and font sizes
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 13}

matplotlib.rc('font', **font)


## in the figures, two timeperiods are plotted; comment out according to desired plot
# days 15-120
start_time= 360
end_time=  2856
jump=1 

# initial phase: days 0-15
#start_time= 0
#end_time=   360
#jump=1






 
 
############## MAIN PROGRAMM ######################  


for var in varlist:
  for bits in range(52,51,-2):
    plt.subplot(2, 1, 1)
    errors= np.zeros((3, 17))
    count=0
    max_iteration=0

    residuum_infnorm= np.zeros((int((end_time-start_time)/jump)+1,100))
    for time in range(start_time, end_time+1,jump):
      plotted= False

      residuum_norm= np.zeros((100))

      exit= np.zeros((100))
      for iteration in range(0,100,1):

        filenameT = path+'Precon'+Precon+'_' + var+'_exp'+exp+'_time'+str(time)+'_iter_'+str(iteration)+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'.txt'

        if (os.path.exists(filenameT)):
          xcoord, ycoord, Residuum, exitcon   =np.loadtxt( filenameT, usecols=(0,1,2,3), unpack=True)

          ncols, nrows = len(set(xcoord)), len(set(ycoord))
          grid = np.flipud(Residuum.reshape((nrows, ncols), order='F'))
          residuum_norm[iteration]=np.linalg.norm (Residuum)

          residuum_infnorm[int((time-start_time)/jump),iteration]=np.linalg.norm (Residuum , np.inf)

          exit[iteration]=np.linalg.norm (exitcon , np.inf)

          max_iteration=max(iteration,max_iteration)

        else:
          residuum_infnorm[int((time-start_time)/jump),iteration]=-1.0

    number_of_iters=0
    for iteration in range(0,100,1):
      if(max(residuum_infnorm[:,iteration]) != min(residuum_infnorm[:,iteration])):
        number_of_iters=iteration+1
    print(number_of_iters)

    for time in range(int((0)/jump), int((end_time-start_time)/jump)+1,int(jump/jump)):

      residuum_infnorm[time,:]=residuum_infnorm[time,:]/(residuum_infnorm[time,0]*(1.0*10**(-10)))

    Minresiduum_infnorm= np.zeros((number_of_iters))
    Maxresiduum_infnorm= np.zeros((number_of_iters))
    Meanresiduum_infnorm= np.zeros((number_of_iters))

    for iteration in range(0,number_of_iters,1):
      Maxresiduum_infnorm[iteration]=max(residuum_infnorm[:,iteration])
    #print('max')
    #print(Maxresiduum_infnorm)
    for iteration in range(0,number_of_iters,1):
       Minresiduum_infnorm[iteration]=min(residuum_infnorm[:,iteration])
       if(Minresiduum_infnorm[iteration] >0.0):
         number_of_iters_min=iteration+1
    #print('min')
    #print(Minresiduum_infnorm)

    for iteration in range(0,number_of_iters,1):
       Meanresiduum_infnorm[iteration]=np.median(residuum_infnorm[:,iteration])
       if(Meanresiduum_infnorm[iteration] >0.0):
         number_of_iters_mean=iteration+1     
    #print('mean')
    #print(Meanresiduum_infnorm)       
       
    number_of_iters=np.arange(0,number_of_iters,1)
    plt.plot(number_of_iters, Maxresiduum_infnorm[:], 'g--', label ='Max') 

    number_of_iters=np.arange(0,number_of_iters_mean,1)
    plt.plot(number_of_iters, Meanresiduum_infnorm[:number_of_iters_mean],'r-', label ='Median')

    number_of_iters=np.arange(0,number_of_iters_min,1)
    plt.plot(number_of_iters, Minresiduum_infnorm[:number_of_iters_min],'b--', label ='Min') 

    plt.xlabel('Solver Iterations')
    plt.ylabel('Linf Residual')
    plt.xticks(np.arange(0,14,1))

    plt.plot((0.0,float(13) ), (1.0,1.0), 'k-')        # 10**0 line means convergence to epsilon=1.0*10**(-10)

    plt.yscale('log')
    plt.ylim(ymin = 10**(-3))
    plt.axis
    plt.savefig(plotpath +'Norm_Errors_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'_'+str(end_time)+'unnorm_WoLegend.png')
    plt.close()
 
    print( plotpath +'Norm_Errors_Precon'+Precon+'_' + var+'_exp'+exp+'_codes_'+str(codesD)+str(codesQ)+'_bits'+str(bits)+'_'+str(end_time)+'unnorm_WoLegend.png')



