import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

font = {'family' : 'normal', 'weight' : 'normal', 'size' : 13} 
matplotlib.rc('font', **font)

data1='model2x2_1x0_l5w200_METRICnorm_Lcoeff_R_uselinear_Lats_ext_DP_noinputact.txt'
data2='model2x2_1x0_l1w5_METRICnorm_Lcoeff_R_uselinear_Lats_ext_DP_noinputact.txt'
data3='model2x2_1x0_l0w1_METRICnorm_Lcoeff_R_uselinear_Lats_ext_DP_noinputact.txt'

data4='model1x1_1x0_l5w200_METRICnorm_Lcoeff_R_uselinear_Lats_ext_DP_noinputact.txt'
data5='model1x1_1x0_l1w5_METRICnorm_Lcoeff_R_uselinear_Lats_ext_DP_noinputact.txt'
data6='model1x1_1x0_l0w1_METRICnorm_Lcoeff_R_uselinear_Lats_ext_DP_noinputact.txt'

data7='model2x2_1x0_l5w200_METRICnorm_Lcoeff_R_uselinear_Lats_ext_DP_noinputact_tanh.txt'

path='../data/data_ADI_NOTgcr_EXP3_Dp_1M10_res05'
file_grid = path+'/Precon5_H_exp3_time0_codes_FT_bits52.txt'
xcoord, ycoord=np.loadtxt( file_grid, usecols=(0,1), unpack=True)

latitudes=np.zeros(32)
latitudes=sorted(set(ycoord), reverse=True)
print(latitudes)

plt.ylim(3*10.0**(-5.0), 2*10.0**(-1.0))
Lat1, base1, baseval1, remainder1   =np.loadtxt( "truth_diff_METRICnorm_iter_1_pseut_Lats.txt", usecols=(0,1,2,3), unpack=True)
plt.plot(latitudes, remainder1[0:32]/base1[0:32], label ='Implicit Richardson, Iteration 1')

Lat2, base2, remainder2   =np.loadtxt( data1, usecols=(0,1,3), unpack=True)
plt.plot(latitudes, remainder2[0:32]/base2[0:32], label ='L5N200, 5x5, 250 Epochs',color='r', linestyle='dotted')


Lat2, base2, remainder2   =np.loadtxt( data2, usecols=(0,1,3), unpack=True)
plt.plot(latitudes, remainder2[0:32]/base2[0:32], label ='L1N5, 5x5',color='r', linestyle='--')
Lat2, base2, remainder2   =np.loadtxt( data3, usecols=(0,1,3), unpack=True)
plt.plot(latitudes, remainder2[0:32]/base2[0:32], label ='L0N0, 5x5',color='r', linestyle='-')

Lat2, base2, remainder2   =np.loadtxt( data6, usecols=(0,1,3), unpack=True)
plt.plot(latitudes, remainder2[0:32]/base2[0:32], label ='L0N0, 3x3',color='g', linestyle='-')

#Lat2, base2, remainder2   =np.loadtxt( data1, usecols=(0,1,3), unpack=True)
#plt.plot(latitudes[14], remainder2[32]/base2[32], 'r*',label ='L5N200, 5x5; 5xEpochs')

#Lat2, base2, remainder2   =np.loadtxt( data7, usecols=(0,1,3), unpack=True)
#plt.plot(latitudes, remainder2[0:32]/base2[0:32], label ='L5N200, 5x5; tanh',color='C1', linestyle='dotted')

plt.xlabel("Latitude")
plt.ylabel("Relative MAE Decrease")

plt.legend(loc='upper left', prop={'size':9}, bbox_to_anchor=(0.2,1.06), framealpha=1)

plt.yscale('log')
plt.savefig('lineplot_iter1_comparison.pdf', bbox_inches=0)


  
