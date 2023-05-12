import matplotlib.pyplot as plt
import torch
import pandas as pd

def checkInput(train_path,test_path):
    train=torch.load(train_path)
    test=torch.load(test_path)
    print("There are",train.shape[1]+test.shape[1], "events.")
    return None
 
    
def dCorr(var_1, var_2, normed_weight, power=1):
    """
    Stolen from https://github.com/gkasieczka/DisCo/blob/master/Disco.py
    var_1: First variable to decorrelate (eg mass)
    var_2: Second variable to decorrelate (eg classifier output)
    normed_weight: Per-example weight. Sum of weights should add up to N (where N is the number of examples)
    power: Exponent used in calculating the distance correlation

    var_1, var_2 and normed_weight should all be 1D torch tensors with the same number of entries
    """

    xx = var_1.view(-1, 1).repeat(1, len(var_1)).view(len(var_1),len(var_1))
    yy = var_1.repeat(len(var_1),1).view(len(var_1),len(var_1))
    amat = (xx-yy).abs()

    xx = var_2.view(-1, 1).repeat(1, len(var_2)).view(len(var_2),len(var_2))
    yy = var_2.repeat(len(var_2),1).view(len(var_2),len(var_2))
    bmat = (xx-yy).abs()

    amatavg = torch.mean(amat*normed_weight,dim=1)
    Amat=amat-amatavg.repeat(len(var_1),1).view(len(var_1),len(var_1))\
        -amatavg.view(-1, 1).repeat(1, len(var_1)).view(len(var_1),len(var_1))\
        +torch.mean(amatavg*normed_weight)

    bmatavg = torch.mean(bmat*normed_weight,dim=1)
    Bmat=bmat-bmatavg.repeat(len(var_2),1).view(len(var_2),len(var_2))\
        -bmatavg.view(-1, 1).repeat(1, len(var_2)).view(len(var_2),len(var_2))\
        +torch.mean(bmatavg*normed_weight)

    ABavg = torch.mean(Amat*Bmat*normed_weight,dim=1)
    AAavg = torch.mean(Amat*Amat*normed_weight,dim=1)
    BBavg = torch.mean(Bmat*Bmat*normed_weight,dim=1)

    if power == 1:
        dCorr=(torch.mean(ABavg*normed_weight))/torch.sqrt((torch.mean(AAavg*normed_weight)*torch.mean(BBavg*normed_weight)))
    elif power == 2:
        dCorr=(torch.mean(ABavg*normed_weight))**2/(torch.mean(AAavg*normed_weight)*torch.mean(BBavg*normed_weight))
    else:
        dCorr=((torch.mean(ABavg*normed_weight))/torch.sqrt((torch.mean(AAavg*normed_weight)*torch.mean(BBavg*normed_weight))))**power

    return dCorr


# plot training history stored as dictionary
def plotHistory(history):
    train_loss = history['train_loss']
    train_bce = history['train_bce']
    train_disco = history['train_disco']
    test_disco = history['test_disco']
    test_loss = history['test_loss']
    test_bce = history['test_bce']

    

    fig, axes=plt.subplots(figsize=(20,10))

    axes.plot(train_loss, alpha=0.8, linewidth=4, label='training loss')
    axes.plot(train_disco, alpha=0.8, linewidth=4, label='disco loss')
    axes.plot(train_bce, alpha=0.8, linewidth=4, label='bce loss')
    axes.plot()
    axes.legend()
    axes.set_title('Training Loss')

    fig, axes=plt.subplots(figsize=(20,10))

    axes.plot(test_loss, alpha=0.8, linewidth=4, label='testing loss')
    axes.plot(test_disco, alpha=0.8, linewidth=4, label='disco loss')
    axes.plot(test_bce, alpha=0.8, linewidth=4, label='bce loss')
    axes.plot()
    axes.legend()
    axes.set_title('Testing Loss')

    fig, axes=plt.subplots(figsize=(20,10))
    axes.plot(test_loss, linewidth=4, label='test loss')
    axes.plot(train_loss, linewidth=4, label='train loss')
    axes.legend()
    axes.set_title("Loss")