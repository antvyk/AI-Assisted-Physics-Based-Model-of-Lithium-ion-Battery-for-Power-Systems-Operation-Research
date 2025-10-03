#NN model of LIBESS based on SPM with SEI

import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import pandas as pd 

DTYPE = 'float32'
tf.keras.backend.set_floatx(DTYPE)

def plot_loss(history):
  
  plt.plot(history.history['val_loss'], label='val_loss')
  plt.ylim([0, 1])
  plt.xlabel('Epoch')
  plt.ylabel('Error [Final SoC]')
  plt.legend()
  plt.grid(True)
  
  



dir_NN = 
dir_Raw_Data = 


column_names = ['Initial SoC', 'Power', 'Initial_SoC_Max','Final SoC','Degradation','Status']
raw_dataset_ch1 = pd.read_csv(dir_Raw_Data+"Charging_1.csv", names=column_names)
raw_dataset_ch2 = pd.read_csv(dir_Raw_Data+'Charging_2.csv', names=column_names)
raw_dataset_ch3 = pd.read_csv(dir_Raw_Data+"Charging_3.csv", names=column_names)
raw_dataset_ch4 = pd.read_csv(dir_Raw_Data+"Charging_4.csv", names=column_names)
raw_dataset_ch5 = pd.read_csv(dir_Raw_Data+"Charging_5.csv", names=column_names)
raw_dataset_ch6 = pd.read_csv(dir_Raw_Data+"Charging_6.csv", names=column_names)
raw_dataset_ch7 = pd.read_csv(dir_Raw_Data+'Charging_7.csv', names=column_names)
raw_dataset_ch8 = pd.read_csv(dir_Raw_Data+"Charging_8.csv", names=column_names)
raw_dataset_ch9 = pd.read_csv(dir_Raw_Data+"Charging_9.csv", names=column_names)
raw_dataset_ch10 = pd.read_csv(dir_Raw_Data+"Charging_10.csv", names=column_names)
raw_dataset_dis1 = pd.read_csv(dir_Raw_Data+"Discharging_1.csv", names=column_names)
raw_dataset_dis2 = pd.read_csv(dir_Raw_Data+"Discharging_2.csv", names=column_names)
raw_dataset_dis3 = pd.read_csv(dir_Raw_Data+"Discharging_3.csv", names=column_names)
raw_dataset_dis4 = pd.read_csv(dir_Raw_Data+"Discharging_4.csv", names=column_names)
raw_dataset_dis5 = pd.read_csv(dir_Raw_Data+"Discharging_5.csv", names=column_names)
raw_dataset_dis6 = pd.read_csv(dir_Raw_Data+"Discharging_6.csv", names=column_names)
raw_dataset_dis7 = pd.read_csv(dir_Raw_Data+"Discharging_7.csv", names=column_names)
raw_dataset_dis8 = pd.read_csv(dir_Raw_Data+"Discharging_8.csv", names=column_names)
raw_dataset_dis9 = pd.read_csv(dir_Raw_Data+"Discharging_9.csv", names=column_names)
raw_dataset_dis10 = pd.read_csv(dir_Raw_Data+"Discharging_10.csv", names=column_names)
raw_dataset = pd.concat([raw_dataset_ch1,raw_dataset_dis1,raw_dataset_ch2,raw_dataset_dis2,
                         raw_dataset_ch3,raw_dataset_dis3,
                         raw_dataset_ch4,raw_dataset_dis4,
                         raw_dataset_ch5,raw_dataset_dis5,
                         raw_dataset_ch6,raw_dataset_dis6,
                         raw_dataset_ch7,raw_dataset_dis7,
                         raw_dataset_ch8,raw_dataset_dis8,
                         raw_dataset_ch9,raw_dataset_dis9,
                         raw_dataset_ch10,raw_dataset_dis10],
                        
                        )

raw_dataset = raw_dataset.reset_index(drop=True)


raw_dataset.loc[raw_dataset['Degradation'] == 0, 'Degradation'] = 1.8 #Degradation output: change from 0 to 3
raw_dataset.loc[raw_dataset['Final SoC'] == 3, 'Final SoC'] = 1.8 #Degradation output: change from 0 to 3

dataset = raw_dataset.copy()
train_dataset = dataset
test_dataset = train_dataset



train_features = train_dataset.copy()
train_labels = train_features[['Final SoC','Degradation']].copy()
weights_train = train_features[['Status']].copy()
train_features = train_features.drop(['Final SoC','Degradation','Status'], axis=1) 


weights_train .loc[weights_train ['Status'] == 1, 'Status'] = 0.2
weights_train .loc[weights_train ['Status'] == 0, 'Status'] = 1

 



normalizer = tf.keras.layers.Normalization(axis=-1)

def build_and_compile_model(norm):
  model = tf.keras.Sequential([
      norm,
      tf.keras.layers.Dense(12, activation='relu'),
      tf.keras.layers.Dense(12, activation='relu'),
      #tf.keras.layers.Dense(8, activation='relu'),
      tf.keras.layers.Dense(2)
  ])

  model.compile(loss='mean_absolute_error', optimizer=tf.keras.optimizers.Adam(0.0005)) #before 0.001
  return model

model = build_and_compile_model(normalizer)


history = model.fit(
    train_features,
    train_labels,
    sample_weight=weights_train,
    validation_split=0,
    verbose=1, epochs=300)


#plot_loss(history)

# plt.plot(history.history['loss'], label='loss')
# plt.ylim([0, 0.2])
# plt.xlabel('Epoch')
# plt.ylabel('Error [Final SoC]')
# plt.legend()
# plt.grid(True)
  
weights_as_numpy = model.get_weights()


# np.savetxt(dir_NN+'w1.csv', weights_as_numpy[3], delimiter=",")
# np.savetxt(dir_NN+'b1.csv', weights_as_numpy[4], delimiter=",")
# np.savetxt(dir_NN+'w2.csv', weights_as_numpy[5], delimiter=",")
# np.savetxt(dir_NN+'b2.csv', weights_as_numpy[6], delimiter=",")
# np.savetxt(dir_NN+'w3.csv', weights_as_numpy[7], delimiter=",")
# np.savetxt(dir_NN+'b3.csv', weights_as_numpy[8], delimiter=",")


# np.savetxt(dir_NN+'w4.csv', weights_as_numpy[9], delimiter=",")
# np.savetxt(dir_NN+'b4.csv', weights_as_numpy[10], delimiter=",")

