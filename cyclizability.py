from keras.models import Sequential
from keras.models import Model
from keras.layers import Dropout, Flatten, BatchNormalization, TimeDistributed, Input, Add
from keras.layers import Dense, Conv2D, MaxPooling2D, LSTM, TimeDistributed, Reshape
import keras.backend as K

import pandas as pd
import numpy as np
from numpy import array
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

from Bio import SeqIO

## define DL model

def model_lstm():
    inputs = Input(shape=(50, 1, 4))
    x = Conv2D(64, (3, 1),
               padding='same',
               activation='relu')(inputs)
    x = MaxPooling2D(pool_size=(2, 1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.2)(x)

    x = Conv2D(16,
               kernel_size=(8,1),
               activation='relu',
               padding='same')(x)
    x = MaxPooling2D((2,1), padding='same')(x)
    x = BatchNormalization()(x)
    x = Dropout(0.2)(x)

    x = Conv2D(8,
               kernel_size=(13,1),
               activation='relu',
               padding='same')(x)
    x = MaxPooling2D((3, 1),padding='same')(x)
    x = BatchNormalization()(x)
    x = Dropout(0.2)(x)

    x = Reshape((K.int_shape(x)[1], K.int_shape(x)[3]))(x)
    x = LSTM(20, return_sequences=False)(x)
    outputs = Dense(1, activation='linear')(x)
    network = Model(inputs, outputs)
    network.compile(optimizer='rmsprop',
                    loss='mean_squared_error')
    return network

## one-hot function
    
def dnaOneHot(sequence):

    seq_array = array(list(sequence))
    
    code={"A":[0],"C":[1],"G":[2],"T":[3],"N":[4],
         "a":[0],"c":[1],"g":[2],"t":[3],"n":[4]}
    onehot_encoded_seq = []
    for char in seq_array:
        onehot_encoded = np.zeros(5)
        onehot_encoded[code[char]] = 1
        onehot_encoded_seq.append(onehot_encoded[0:4])
    
    return onehot_encoded_seq

## data to train 
# X6 is an array of sequences
# Y6 is numerical

data_chr5 = pd.read_csv("/Users/keren/Documents/HumanChecm/Cyclization/cycle6.txt",delimiter = ",")
X6 = []
for sequence_nt in data_chr5["Sequence"]:
    X6.append(dnaOneHot(sequence_nt)) 
X6 = array(X6)    
X6 = X6.reshape((X6.shape[0],50,1,4))
Y6 = data_chr5["C0"]    

## train model

network_lstm_chrv = model_lstm()
network_lstm_chrv.fit(X6,Y6,
                   epochs=80,
                   batch_size=128,
                   validation_split=0.2)    

## Mus musculus genome sequences

genome_file = SeqIO.parse(open('/Users/keren/Documents/HumanChecm/Cyclization/fa/Mus_musculus.GRCm39.dna.toplevel.fa'),'fasta')
cycle_path = "/Users/keren/Documents/HumanChecm/Cyclization/mus_cycle/"

for fasta in genome_file:
    chrom = fasta.id
    print(chrom)
    genome_sequence = str(fasta.seq)
    onehot_sequence = dnaOneHot(genome_sequence)
    onehot_sequence = array(onehot_sequence)
    onehot_sequence = onehot_sequence.reshape((onehot_sequence.shape[0],1,4))
    print(onehot_sequence.shape)
    fit = []
    for ind_local in np.array_split(range(25, onehot_sequence.shape[0]-24), 100):
        onehot_sequence_local = []
        for i in ind_local:
            s = onehot_sequence[(i-25):(i+25),]
            onehot_sequence_local.append(s)
        onehot_sequence_local = array(onehot_sequence_local)
        onehot_sequence_local = onehot_sequence_local.reshape((onehot_sequence_local.shape[0],50,1,4))
        fit_local = network_lstm_chrv.predict(onehot_sequence_local)
        fit_local = fit_local.reshape((fit_local.shape[0]))
        fit.append(fit_local)
    fit = [item for sublist in fit for item in sublist]
    fit = array(fit)
    n = fit.shape[0]
    fitall = np.vstack((range(25,25+n),fit))
    fitall = pd.DataFrame([*zip(*fitall)])
    fitall.columns = ["pos","C0"]
    fitall.to_csv(cycle_path+"chr_"+chrom+".txt", index = False)

