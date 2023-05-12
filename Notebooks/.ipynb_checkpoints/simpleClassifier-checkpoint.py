import numpy as np
import torch.nn as nn
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader
import torch
from sklearn.metrics import roc_curve,roc_auc_score
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Take in options')
parser.add_argument('--epoch', default=100, type=int,
                    help='number of epochs to train')
parser.add_argument('--batchsize', default=2000, type=int,
                    help='number of samples per batch')
parser.add_argument('--lr', default=0.001, type=float,
                   help='learning rate')

args = parser.parse_args()

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


size=2000000
sig=np.random.multivariate_normal(mean_sig,covariance_sig,size)
bkg=np.random.multivariate_normal(mean_bkg,covariance_bkg,size)

X=np.concatenate((sig[:,1:],bkg[:,1:]))
Y=np.concatenate((np.ones(size),np.zeros(size)))
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
batch_size= args.batchsize
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

learning_rate = args.lr
epochs = args.epoch
optimizer = torch.optim.Adam(model.parameters(),lr=learning_rate)
loss_fn = nn.BCELoss()

model.train()
loss_history=np.zeros(epochs)
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
    loss_history[i-1]=epoch_loss/len(trainloader)

plt.figure()
plt.plot(loss_history)
plt.title("loss history")
plt.xlabel("epoch")
plt.ylabel("loss")
plt.savefig('classifierHist.png')

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.9, random_state=42)
testdata = Data(X_test, Y_test)
batch_size=len(Y_test)
testloader = DataLoader(testdata, batch_size=batch_size, 
                         shuffle=True, num_workers=1)
y_pred_list = []
y_true_list = []
model.eval()
with torch.no_grad():
    for X_batch, Y_batch in testloader:
        y_test_pred = model(X_batch)
        y_pred_list.append(y_test_pred.cpu().numpy())
        y_true_list.append(Y_batch.cpu().numpy())

y_pred_list = [a.squeeze().tolist() for a in y_pred_list]
y_true_list = [a.squeeze().tolist() for a in y_true_list]

truth_arr=np.asarray(y_true_list[0])
pred_arr=np.asarray(y_pred_list[0])
fpr, tpr, thresholds = roc_curve(truth_arr, pred_arr)
auc = roc_auc_score(truth_arr, pred_arr)

plt.figure()
plt.plot(fpr,tpr,label="AUC="+str(auc))
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(loc=4)
plt.title("ROC curve")
plt.savefig('classifierROC.png')

bkg_y_pred=y_pred_list[y_true_list==0]


