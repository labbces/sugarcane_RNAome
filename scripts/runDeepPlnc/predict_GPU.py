#!/usr/bin/env python

import os, sys
os.environ["CUDA_VISIBLE_DEVICES"] = "0" #put '-1' instead of '0' if u don't want to use GPU
import tensorflow as tf, numpy, keras, tensorflow, numpy as np
from keras_preprocessing.sequence import pad_sequences 
from keras.models import load_model
from multiprocessing import Pool #this package will help in faster loading of data

def mono_hot_encode1(seq):
	mapping = dict(zip(['.',')','('], range(3)))  
	text=[seq[i:i+1] for i in range(len(seq)-(1-1))]
	seq1 = [mapping[i] for i in text]
	return np.eye(3)[seq1]


def mono_hot_encode(seq):
	mapping = dict(zip(['A','T','G','C'], range(4)))  
	text=[seq[i:i+1] for i in range(len(seq)-(1-1))]
	seq1 = [mapping[i] for i in text]
	return np.eye(4)[seq1]


filename=str(sys.argv[1]) #change your file here
model = load_model("Model_A.h5") #model file
def s1(k):
	texts_mono_fold = []
	texts_mono_fold.append(mono_hot_encode1(test12[k]))
	padded_docs3 = pad_sequences(texts_mono_fold, maxlen=400, padding='post') 
	texts_mono = []
	lab = testlab[k]
	texts_mono.append(mono_hot_encode(testseq[k]))
	padded_docs4 = pad_sequences(texts_mono, maxlen=400, padding='post')
	return lab, padded_docs4, padded_docs3


print("Step 1: Library, Functions and Model loaded")

testseq=[]
testlab=[];test12=[]
for i in open(filename):
	z=i.split()
	testlab.append(str(z[0]))
	testseq.append(str(z[1]))
	test12.append(z[2])


name=[]
name = Pool().map(s1, [sent for sent in range(len(testseq))])
input_ids=[]
attention_masks=[]
attention_masks_fold=[]
for i, j in enumerate(name):
	input_ids.append(name[i][0])	
	attention_masks.append(name[i][1])
	attention_masks_fold.append(name[i][2])

train_inp = numpy.array(attention_masks).reshape(len(numpy.array(attention_masks)),400,4,1);attention_masks=[]
train_label = numpy.array(input_ids);input_ids=[]
train_fold = numpy.array(attention_masks_fold).reshape(len(numpy.array(attention_masks_fold)),400,3,1);attention_masks_fold=[]
print("Step 2: File loaded")
y_pred = model.predict([train_inp,train_fold], verbose=1)
f1=open(filename+"_prediction.txt",'w')
for j,i in enumerate(y_pred):
	f1.writelines(str(train_label[j])+"\t"+str(i)+"\n")


f1.close()
print("Step 3: Process Complete\n\nYours results is in "+filename+"_prediction.txt file")
