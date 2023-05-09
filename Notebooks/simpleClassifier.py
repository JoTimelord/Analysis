import numpy as np
import torch.nn as nn
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader
import torch
from sklearn.preprocessing import StandardScaler


mean_sig = np.array([2.5, 2.5, 2])
variance_sig = 1.5
covariance_sig = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
])*(variance_sig)**2


mean_bkg = np.array([0, 0, 0])
variance_bkg = 1.5
covariance_bkg = np.array([
    [1, 0.8, 0],
    [0.8, 1, 0],
    [0, 0, 1]
])*(variance_bkg**2)


size=1000
sig=np.random.multivariate_normal(mean_sig,covariance_sig,size)
bkg=np.random.multivariate_normal(mean_bkg,covariance_bkg,size)

X=np.concatenate((sig[:,1:],bkg[:,1:]))
Y=np.concatenate((np.ones(size),np.zeros(size)))
sc = StandardScaler()
X = sc.fit_transform(X)
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.15, random_state=42)

class Data(Dataset):
  def __init__(self, X_train, y_train):
    self.X = torch.from_numpy(X_train.astype(np.float32))
    self.y = torch.from_numpy(y_train.astype(np.float32))
    self.len = self.X.shape[0]
  
  def __getitem__(self, index):
    return self.X[index], self.y[index]
  def __len__(self):
    return self.len
  
traindata = Data(X_train, Y_train)
batch_size=10
trainloader = DataLoader(traindata, batch_size=batch_size, 
                         shuffle=True, num_workers=1)
input_size=2
hidden_sizes=128
output_size=1
model = nn.Sequential(nn.Linear(input_size, hidden_sizes), # 1st hidden
                      nn.ReLU(),
                      nn.Linear(hidden_sizes, hidden_sizes), # 2nd hidden
                      nn.ReLU(),
                      nn.Linear(hidden_sizes, hidden_sizes), # 3rd hidden
                      nn.ReLU(),
                      nn.Linear(hidden_sizes, output_size), # output layer
                      nn.Sigmoid()
                     )

learning_rate = 0.01
epochs = 100
optimizer = torch.optim.Adam(model.parameters(),lr=learning_rate)
loss_fn = nn.BCELoss()

model.train()
for i in range(1,epochs+1):
    epoch_loss = 0
    for j,(x_train,y_train) in enumerate(trainloader):
        optimizer.zero_grad()
        
        #calculate output
        output = model(x_train)

        #calculate loss
        loss = loss_fn(output,y_train.unsqueeze(1))

        #backprop
        loss.backward()
        optimizer.step()
        epoch_loss += loss.item()

    if i%5==0 or i==1: print(f'Epoch {i+0:03}: | Loss: {epoch_loss/len(trainloader):.5f}')