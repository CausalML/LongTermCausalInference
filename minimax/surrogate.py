#!/usr/bin/env python
# coding: utf-8

# In[6]:


import numpy as np
import random
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch
import torchvision
from scipy.stats import multivariate_normal
import sys


#
#Change the dimension of X, S as well as the number of 
#neurons in the first hidden layer numneuron according 
#to the setting specified in each cell of the table; 
#for example, if you would like to replicate the simulation 
#in the bottom-left cell, then you have to set dimX = 5,
#dimS = 5 and set numneuron = 30. The choice of number 
#of neurons in each result is described in the first 
#paragraph of Section F.3 in the Supplementary Material.
#
dimX = 10
dimS = 5
numneuron = 50

nsim = 20
sample_size = 2000
random_seed = 1

random.seed(random_seed)

momentum = 0.95
log_interval = 10
dimX = dimX + 1 #here "dimX" refers to the dimension of X plus the dimension of S2

#batch = int(sample_size/40)
#learning_rate = 0.0002
#n_epochs = 20
batch = int(sample_size/10)
learning_rate = 0.0002
n_epochs = 40

idx = sys.argv[1]

#the order of random variables in the csv file is: X, S2, S1, S3, Y, A

###get experimental and observational data
data = np.genfromtxt("tmp/obs_" + idx + ".csv", delimiter=',', dtype='f')
obsX = data[:,:dimX]##this includes both X and S2
obsS1 = data[:,dimX:(dimX + dimS)]
obsS3 = data[:,(dimX + dimS):(dimX + 2 * dimS)]
obsY = data[:,(dimX + 2 * dimS):(dimX + 2 * dimS + 1)]
obsA = data[:,(dimX + 2 * dimS + 1):(dimX + 2 * dimS + 2)]

data = np.genfromtxt("tmp/exp_" + idx + ".csv", delimiter=',', dtype='f')
expX = data[:, :dimX]##this includes both X and S2
expS1 = data[:, dimX:(dimX + dimS)]
expS3 = data[:,(dimX + dimS):(dimX + 2 * dimS)]
expY = data[:,(dimX + 2 * dimS):(dimX + 2 * dimS + 1)]
expA = data[:,(dimX + 2 * dimS + 1):(dimX + 2 * dimS + 2)]


# In[15]:



x_list_torch  = torch.from_numpy(obsX[obsA.flatten() ==1])
z_list_torch  = torch.from_numpy(obsS1[obsA.flatten() ==1] )
w_list_torch  = torch.from_numpy(obsS3[obsA.flatten() ==1] )
y_list_torch  = torch.from_numpy(obsY[obsA.flatten() ==1] )
        
        
torch.backends.cudnn.enabled = False
torch.manual_seed(random_seed)

class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
            
        #self.fc1 = nn.Linear(dimX + dimS, dimX + dimS)
        #self.fc2 = nn.Linear(dimX + dimS, 10)
        self.fc1 = nn.Linear(dimX + dimS, numneuron)
        self.fc2 = nn.Linear(numneuron, 10)
        self.fc3 = nn.Linear(10, 1)

    def forward(self, x):
        ####x = x.view(-1, 320)        
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = (self.fc3(x))
        ###x = self.fc3(x)
        return x

def kernel_output(data1,data2,data3,data4,bandwidth,bandwidth2):

    loss = torch.tensor(0)
    kkk = (data1-data2)
    gram = torch.zeros([batch, batch])
    #####print(kkk.size())
    for i in range(batch):
        for j in range(batch):
            gram[i,j]=torch.exp(-torch.sum((data3[i,:]-data3[j,:])*(data3[i,:]-data3[j,:]))/(2.0*bandwidth)-torch.sum((data4[i,:]-data4[j,:])*(data4[i,:]-data4[j,:]))/(2.0*bandwidth2) )
    loss = torch.matmul(torch.matmul(torch.transpose(kkk, 0, 1),gram),kkk)
    return loss/(batch*batch*1.0)

def kernel_output_stab(data1,data2,data3,data4,bandwidth,bandwidth2):

    loss = torch.tensor(0)
    kkk = (data1-data2)
    gram = torch.zeros([batch, batch])
    #####print(kkk.size())
    for i in range(batch):
        for j in range(batch):
            gram[i,j]=torch.exp(-torch.sum((data3[i,:]-data3[j,:])*(data3[i,:]-data3[j,:]))/(2.0*bandwidth)-torch.sum((data4[i,:]-data4[j,:])*(data4[i,:]-data4[j,:]))/(2.0*bandwidth2) )
    l, V = torch.linalg.eig(gram)
    #gram = V @ torch.diag(l / (5 + l)) @ torch.linalg.inv(V)
    gram = V @ torch.diag(l / 5) @ torch.linalg.inv(V)
    #gram = V @ torch.diag(l / (1 + 0.01 * l)) @ torch.linalg.inv(V)
    gram = gram.real
    loss = torch.matmul(torch.matmul(torch.transpose(kkk, 0, 1), gram),kkk)
    return loss/(batch*batch*1.0)
       
    
    
def train(epoch, train_losses,sample_size2,var1,var2):
    network.train()
    listlist = np.random.permutation(sample_size2)
    loss_track = 0.0 

    for i in range(int(sample_size2/batch)):
        optimizer.zero_grad()
        data1 = x_list_torch[listlist[i*batch:(i+1)*batch],:]
        data3 = w_list_torch[listlist[i*batch:(i+1)*batch],:]
        data4 = z_list_torch[listlist[i*batch:(i+1)*batch],:]
        data = torch.cat((data1,data3),1)
        output = network(data)
        target = y_list_torch[listlist[i*batch:(i+1)*batch],:]
        #loss =  kernel_output(output, target,data1,data4,var1,var2) 
        loss =  kernel_output_stab(output, target,data1,data4,var1,var2) 
        loss.backward()
        optimizer.step()
        loss_track  =   loss_track  + loss.item()
    train_losses.append(loss_track/(1.0*sample_size2))
            
    
        
network = Net()
optimizer = optim.RMSprop(network.parameters(), lr=learning_rate, momentum=momentum)
var1 = torch.var(x_list_torch).item() #change into obsX
var2 = torch.var(z_list_torch).item() #change into obsS1
train_losses = []

loss_list =[]
sample_size2 = y_list_torch.size()[0]
        
for epoch in range(1, n_epochs + 1):
    train(epoch,  train_losses,sample_size2,var1,var2)

network.eval()
with torch.no_grad():
    data_pred = torch.cat((torch.from_numpy(expX[expA.flatten() == 1]),torch.from_numpy(expS3[expA.flatten() == 1])),1)
    pred_exp = network(data_pred)
    data_pred = torch.cat((torch.from_numpy(obsX[obsA.flatten() == 1]),torch.from_numpy(obsS3[obsA.flatten() == 1])),1)
    pred_obs = network(data_pred)
    #pred_reg[iii,count] = torch.mean(pred).item()
            
#print(np.mean((pred_list_sample_dm[0:iii+1,count]-answer)*(pred_list_sample_dm[0:iii+1,count]-answer)))
        
        
########## IP
batch = int(sample_size/10)

learning_rate = 0.0002
n_epochs = 40
momentum = 0.95

#batch = int(sample_size/40)
#learning_rate = 0.0002
#n_epochs = 20
        
x_list_torch  = torch.from_numpy(np.concatenate((obsX[obsA.flatten() == 1], expX[expA.flatten() == 1]), axis=0))
z_list_torch  = torch.from_numpy(np.concatenate((obsS1[obsA.flatten() == 1], expS1[expA.flatten() == 1]), axis=0))
w_list_torch  = torch.from_numpy(np.concatenate((obsS3[obsA.flatten() == 1], expS3[expA.flatten() == 1]), axis=0))
#a_list_torch  = torch.from_numpy(np.concatenate(np.onesobsY[obsA.flatten() == 1], expX[expA.flatten() == 1], 2))
a_list_torch  = torch.from_numpy(np.concatenate((np.ones_like(obsY[obsA.flatten() == 1]), np.zeros_like(expY[expA.flatten() == 1])), axis=0))
y_list_torch  = torch.from_numpy(np.concatenate((obsY[obsA.flatten() == 1], expY[expA.flatten() == 1]), axis=0))
        
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
            
        self.fc1 = nn.Linear(dimX + dimS, numneuron)
        self.fc2 = nn.Linear(numneuron, 10)
        self.fc3 = nn.Linear(10, 1)

    def forward(self, x):
        ####x = x.view(-1, 320)        
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.softplus(self.fc3(x))
        ###x = self.fc3(x)
        return x
        
        

def kernel_output(data1,data2,data3,data4,bandwidth,bandwidth2):

    loss = torch.tensor(0)
    kkk = (data1-data2)
    gram = torch.zeros([batch, batch])
    for i in range(batch):
        for j in range(batch):
            gram[i,j]=torch.sum((data3[i,:]*data3[j,:])+ torch.sum((data4[i,:]*data4[j,:])))
    loss = torch.matmul(torch.matmul(torch.transpose(kkk, 0, 1),gram),kkk)
    return loss/(batch*batch*1.0)

def kernel_output_stab(data1,data2,data3,data4,bandwidth,bandwidth2):

    loss = torch.tensor(0)
    kkk = (data1-data2)
    gram = torch.zeros([batch, batch])
    for i in range(batch):
        for j in range(batch):
            gram[i,j]=torch.sum((data3[i,:]*data3[j,:])+ torch.sum((data4[i,:]*data4[j,:])))
    l, V = torch.linalg.eig(gram)
    #gram = V @ torch.diag(l / (1 + 0.05 * l)) @ torch.linalg.inv(V)
    gram = V @ torch.diag(l / 5) @ torch.linalg.inv(V)
    gram = gram.real
    loss = torch.matmul(torch.matmul(torch.transpose(kkk, 0, 1), gram),kkk)
    return loss/(batch*batch*1.0)

            
def train(epoch, train_losses,sample_size2,var1,var2):
    network2.train()
    listlist = np.random.permutation(sample_size2)
    loss_track = 0.0 

    for i in range(int(sample_size2/batch)):
        optimizer2.zero_grad()
        data1 = x_list_torch[listlist[i*batch:(i+1)*batch],:]
        data3 = w_list_torch[listlist[i*batch:(i+1)*batch],:]
        data6 = a_list_torch[listlist[i*batch:(i+1)*batch],:] #one means observational data
        data4 = z_list_torch[listlist[i*batch:(i+1)*batch],:]
        data = torch.cat((data1,data4),1)
        output = network2(data)
        target = y_list_torch[listlist[i*batch:(i+1)*batch],:]
        #loss =  kernel_output(data6 * output + data6, torch.ones_like(target), data1,data3,var1,var2) 
        loss =  kernel_output_stab(data6 * output + data6, torch.ones_like(target), data1,data3,var1,var2) 
        loss.backward()
        optimizer2.step()
        loss_track  =   loss_track  + loss.item()
    train_losses.append(loss_track/(1.0*sample_size2))
        
        
network2 = Net()
optimizer2 = optim.RMSprop(network2.parameters(), lr=learning_rate, momentum=momentum)
var1 = torch.var(x_list_torch).item()
var2 = torch.var(w_list_torch).item()
train_losses = []


loss_list =[]
sample_size2 = y_list_torch.size()[0]
for epoch in range(1, n_epochs + 1):
    train(epoch,  train_losses,sample_size2,var1,var2)
        
network2.eval()
with torch.no_grad():
    data_pred = torch.cat((torch.from_numpy(obsX[obsA.flatten() == 1]),torch.from_numpy(obsS1[obsA.flatten() == 1])),1)
    pred_prop = network2(data_pred)

y_list_torch  = torch.from_numpy(obsY[obsA.flatten() ==1] )
ate_dr = torch.mean(pred_exp).item() + torch.mean(pred_prop * (y_list_torch - pred_obs)).item()
ate_or = torch.mean(pred_exp).item()

sigma = torch.mean((pred_exp - ate_dr) * (pred_exp - ate_dr)).item() / np.sum(expA) + torch.mean(pred_prop * pred_prop * (y_list_torch - pred_obs) * (y_list_torch - pred_obs)).item() / np.sum(obsA)
sd = np.sqrt(sigma)

# In[18]:
np.savetxt("tmp/result_" + idx + ".csv", np.array([ate_dr, ate_or, sd]), delimiter=",")
