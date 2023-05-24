import matplotlib.pyplot as plt
import torch
import pandas as pd

def checkInput(train_path,test_path):
    train=torch.load(train_path)
    test=torch.load(test_path)
    print("There are",train.shape[1]+test.shape[1], "events.")
    return None
 
# plot training history stored as dictionary
def plotHistory(history):
    train_loss = history['train_loss']
    train_bce = history['train_bce']
    train_disco = history['train_disco']
    test_disco = history['test_disco']
    test_loss = history['test_loss']
    test_bce = history['test_bce']

    plt.style.use('seaborn-talk')

    fig, axes=plt.subplots(figsize=(20,10))
    axes.xticks(fontsize=15)
    axes.yticks(fontsize=15)
    axes.minorticks_on()

    axes.plot(train_loss, alpha=0.8, label='training loss')
    axes.plot(train_disco, alpha=0.8, label='disco loss')
    axes.plot(train_bce, alpha=0.8, label='bce loss')
    axes.plot()
    axes.legend()
    axes.set_title('Training Loss')

    fig, axes=plt.subplots(figsize=(20,10))

    axes.plot(test_loss, alpha=0.8, label='testing loss')
    axes.plot(test_disco, alpha=0.8, label='disco loss')
    axes.plot(test_bce, alpha=0.8, label='bce loss')
    axes.plot()
    axes.legend()
    axes.set_title('Testing Loss')

    fig, axes=plt.subplots(figsize=(20,10))
    axes.plot(test_loss, label='test loss')
    axes.plot(train_loss, label='train loss')
    axes.legend()
    axes.set_title("Loss")