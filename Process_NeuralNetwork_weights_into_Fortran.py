import numpy as np
import keras
from keras.models import Sequential,load_model
from keras.layers import Dense, Activation
from keras.callbacks import EarlyStopping,ModelCheckpoint
from keras.optimizers import rmsprop,adam
from keras.utils import plot_model
import keras.backend as K

def main(latitude=0):
    K.set_floatx('float64')
    print(latitude)
    inputs=176
    inputs1=5


    #Load the best model
    model = load_model('ModelWeights/model2x2_1x0_l0w1_METRICnorm_Lat_Lcoeff_R_uselinear_noinputact_reduced_r_only_r'+str(latitude)+'_ext_DP.h5')
    weights, biases = model.layers[0].get_weights()
    #print(weights)
    #print(biases)
    f = open('FortranWeights_model2x2_1x0_l0w1coeff_R_uselinear_noinputact_reduced_r_only_rLat'+str(latitude)+'_layer0.txt', "a")
    for inp in range(inputs):
      #print(weights[inputs,:])
      for i in range(0,1):
        f.write(str(weights[inp,i])+' ') 
      f.write('\n') 
    for i in range(0,1):
      f.write(str(biases[i])+' ')
    f.write('\n') 
    f.close()

    #weights, biases = model.layers[1].get_weights()
    #print(weights)
    #print(biases)
    #f = open('FortranWeights_model2x2_1x0_l1w5coeff_R_Lat'+str(latitude)+'_layer1.txt', "a")
    #print('FortranWeights_model2x2_1x0_l1w5coeff_R_Lat'+str(latitude)+'_layer1.txt')
    #for inp in range(inputs1):
    #  #print(weights[inputs,:])
    #  for i in range(1):
    #    f.write(str(weights[inp,i])+' ') 
    #  f.write('\n') 
    #for i in range(1):
    #  f.write(str(biases[i])+' ')
    #f.write('\n') 
    #f.close()
      

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

for latitude in range(0,32):
  main(latitude)
