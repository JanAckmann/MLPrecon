import numpy as np
import keras
from keras.models import Sequential,load_model
from keras.layers import Dense, Activation
from keras.callbacks import EarlyStopping,ModelCheckpoint
from keras.optimizers import rmsprop,adam
import keras.backend as K
import matplotlib.pyplot as plt

plotpath='Model_plots/'

def main(activation = 'relu',learning_rate=10.0**(-5),hidden_layers=5,width=200,Lat=9):
    K.set_floatx('float64')
    #n_train =4481792 #448064 #384064 #2689792 #4481792 # 4480000 #1612800 #1843200 #13312000
    #n_valid =1254400 #44800 # 539392 #897792 #1612800 #179200  #204800  #2455552
    inputs= 176 #136 #508 #976 #256 #76 #76  #976 #
    outputs= 1 #65 #


    print('load validation data')
    #x_val, y_val =  load_pair('../data/ML_data_ADI_NOTgcr_EXP3_Dp_res05/validation32x2_32x0.h5')
    x_val, y_val =  load_pair('../data/ML_data_ADI_NOTgcr_EXP3_Dp_1M10_res05/validation_Lcoeff_R_2x2_1x0'+str(Lat)+'DP.h5')
    print(np.sqrt(np.mean(np.square(y_val))) )
    print(np.mean(np.absolute(y_val)) )

    Maxima=np.zeros((inputs)) 
    Minima=np.zeros((inputs))
    Mean=np.zeros((inputs))
    MaximaO=np.zeros((outputs)) 
    MinimaO=np.zeros((outputs))
    MeanO=np.zeros((outputs))

    print(len(y_val))
    Relative_dec=np.zeros((len(y_val))) 
    print(len(y_val))
    ScalingO=np.zeros((len(y_val))) 


    for k in range(0, len(y_val)):
     for i in range(1, inputs,7):
      if( abs(x_val[k,i]) > ScalingO[k] ):
       ScalingO[k]    =abs(x_val[k,i])


    Maxima[0]=3.5/2.0 #max(x[:,0])
    Minima[0]=-3.5/2.0 #min(x[:,0])
    Mean[0]=0.0

    Maxima[1]=21.0 #max(x[:,i])
    Minima[1]=-21.0 #min(x[:,i])  
    Mean[1]=0.0

    Maxima[2]=8.54163233e+11 #21.0 #max(x[:,i])
    Minima[2]=2.48437139e+10
    Mean[2]=1.17612126e+11 #0.0

    Maxima[3]=7.32606531e+08 #21.0 #max(x[:,i])
    Minima[3]=-7.32606531e+08 #-21.0 #min(x[:,i])  
    Mean[3]=0.0

    Maxima[4]=4.34271611e+10 #1.0 #max(x[:,i])
    Minima[4]=1.28429346e+09  
    Mean[4]=2.57373289e+10

    Maxima[5]=7.32606531e+08
    Minima[5]=-7.32606531e+08
    Mean[5]=0.0

    Maxima[6]=5.26020440e+09 #21.0 #max(x[:,i])
    Minima[6]=-5.26020440e+09 #-21.0 #min(x[:,i])  
    Mean[6]=0.0

    Maxima[7]=2.63625191e+09 #21.0 #max(x[:,i])
    Minima[7]=-2.63625191e+09 #-21.0 #min(x[:,i])  
    Mean[7]=0.0


    print('Maxima')
    print (Maxima)
    print('Minima')
    print (Minima)
    print('Mean')
    print (Mean)    

    print ('before')

    print(min(y_val), max(y_val) )

    print(np.mean(y_val) )


    x_val[:,0]=(x_val[:,0] )/(Maxima[0])



    for k in range(0, len(y_val)):    
     for i in range(1, inputs,7):
      x_val[k,i]= x_val[k,i] /ScalingO[k]

    for i in range(2, inputs,7):

      x_val[:,i]=(x_val[:,i] -Mean[2])/(Maxima[2]-Minima[2])
    for i in range(3, inputs,7):

      x_val[:,i]=(x_val[:,i] -Mean[3])/(Maxima[3]-Minima[3])
    for i in range(4, inputs,7):

      x_val[:,i]=(x_val[:,i] -Mean[4])/(Maxima[4]-Minima[4])
    for i in range(5, inputs,7):

      x_val[:,i]=(x_val[:,i] -Mean[5])/(Maxima[5]-Minima[5])
    for i in range(6, inputs,7):

      x_val[:,i]=(x_val[:,i] -Mean[6])/(Maxima[6]-Minima[6])
    for i in range(7, inputs,7):

      x_val[:,i]=(x_val[:,i] -Mean[7])/(Maxima[7]-Minima[7])




    for k in range(0, len(y_val)):    
     y_val[k]= y_val[k] /ScalingO[k]





    print(min(y_val), max(y_val) )
    print('Mean of Testoutput and Validationoutput')

    print(np.mean(y_val) )

    print(np.mean(ScalingO) )

    print(min(ScalingO), max(ScalingO) )

    print ('normalised')

    print (x_val[0,:])
    print (y_val[0,:])
    print(ScalingO[0])

    assert x_val.shape[1] == inputs
    assert y_val.shape[1] == outputs



    #Load the best model
    model = load_model('ModelWeights/model2x2_1x0_l'+str(hidden_layers)+'w'+str(width)+'_METRICnorm_Lat_Lcoeff_R_uselinear_noinputact'+str(Lat)+'_ext_DP.h5')
    #Check it gets the same answer on the validation data
    error=model.evaluate(x_val,y_val)
    y_ML_model = model.predict(x_val)

    for k in range(0, len(y_val)): 
     y_ML_model[k]=y_ML_model[k]*ScalingO[k]
     y_val[k]=y_val[k]*ScalingO[k]
    plt.xlabel(r'$\Delta \Phi$', fontsize=16)
    plt.ylabel("Relative Decrease in Absolute Error", fontsize=16)
    plt.ylim(10.0**(-9.0), 10.0**(4.0))
    Relative_dec[:]= np.absolute((y_ML_model[:,0]-y_val[:,0])/y_val[:,0])
    plt.ylim(10.0**(-9.0), 10.0**(4.0))
    plt.yscale('log')
    plt.ylim(10.0**(-9.0), 10.0**(4.0))
    plt.axhline(y=10.0**(-1.0), color='r', linestyle='--')
    plt.scatter(y_val[:,0], Relative_dec, vmin=10.0**(-9.0), vmax=10.0**(4.0))
    plt.ylim(10.0**(-9.0), 10.0**(4.0))
    plt.savefig(plotpath +'model2x2_1x0_l'+str(hidden_layers)+'w'+str(width)+'_METRICnorm_Lat_Lcoeff_R_uselinear_noinputact'+str(Lat)+'_ext_DP.jpg', bbox_inches=0)
    plt.close()
    print('relative drecrease')
    print(np.mean(Relative_dec) )
    print(np.amax(Relative_dec) )
    print(np.amin(Relative_dec) )
    counter=0
    for i in range(len(Relative_dec)):
      if(Relative_dec[i]>=1.0):
        counter+=1
        print(y_val[i,0])
    print('how many larger than 1')
    print(float(counter))

    print(y_ML_model.shape[0])
    print(np.sqrt(np.mean(np.square(y_ML_model))) )
    print ('MAE of y_ML_model and y')
    print(np.mean(np.absolute(y_ML_model)) )

    print(np.sqrt(np.mean(np.square(y_val-y_ML_model))) )
    print(np.mean(np.absolute(y_val-y_ML_model)) )
    print(error)

    print(np.mean(np.absolute( (y_val[0:5]-y_ML_model[0:5]) )))

    print(np.mean(np.absolute( (y_val[0:1800*64]-y_ML_model[0:1800*64]) )))
    print(np.mean(np.absolute( (y_val[1*1800*64 :2*1800*64]-y_ML_model[1*1800*64 :2*1800*64])  )))
    print(np.mean(np.absolute( (y_val[2*1800*64 :3*1800*64]-y_ML_model[2*1800*64 :3*1800*64])  )))
    print(np.mean(np.absolute( (y_val[3*1800*64 :4*1800*64]-y_ML_model[3*1800*64 :4*1800*64])  )))
    print(np.mean(np.absolute( (y_val[4*1800*64 :]-y_ML_model[4*1800*64 :])  )))


    
    #for i in range(0,700):
      #print(i)
      #print(y_val[i]*(MaximaO[0]-MinimaO[0])+MeanO[0] , y_ML_model[i]*(MaximaO[0]-MinimaO[0])+MeanO[0], y_val[i]*(MaximaO[0]-MinimaO[0])+MeanO[0] - (y_ML_model[i]*(MaximaO[0]-MinimaO[0])+MeanO[0]))
      #print(y_val[i], y_ML_model[i], y_val[i] - y_ML_model[i])
      #raw_input("Press Enter to continue...")

def save_pair(x,y,fname):
    import h5py
    hf = h5py.File(fname, 'w')
    hf.create_dataset('x', data=x)
    hf.create_dataset('y', data=y)
    hf.close()
    return

def load_pair(fname):
    import h5py
    hf = h5py.File(fname, 'r')
    x=hf['x'].value
    y=hf['y'].value
    hf.close()
    return x,y

def relative_Loss(yTrue,yPred):
    return K.mean(K.abs((yPred-yTrue)/yTrue))

for w_width in [ 1]:
  for h_layers in range(0,1,2):
    for Latitude in range(0,32,1):
      w_width
      main(width=w_width,hidden_layers=h_layers,Lat=Latitude)
