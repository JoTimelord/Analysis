import matplotlib.pyplot as plt
import torch
import pandas as pd
import numpy as np
import yahist
from yahist import Hist1D
import mplhep as hep


# calculate a dataframe that stores the background and signal in A+B versus C+D region
# and a dataframe for closure in background only (for decorrelation purpose)
def ABCDClosure(df, truth_name, weight_name, ABCD_name, ABCD_range, ABCD_no, DisCo_name, DisCo_range, DisCo_no):
    result=pd.DataFrame(columns=['f cut', 'g cut', 'A', 'B', 'C', 'D', r'Predicted $A_{BKG}$', r'Actual $A_{BKG}$',
                                r'Actual $A_{SIG}$'])
    decorrelation=pd.DataFrame(columns=['f cut', 'g cut', 'predicted_A_bkg', 'actual_A_bkg', 'ratio'])
    NN_cuts=np.linspace(*ABCD_range, ABCD_no)
    Disco_cuts=np.linspace(*DisCo_range, DisCo_no)
    sig_df=df[df[truth_name]==1]
    bkg_df=df[df[truth_name]==0]
    total_sig=np.sum(sig_df[weight_name])
    total_bkg=round(np.sum(bkg_df[weight_name]))
    total_wgt=total_sig+total_bkg
    i=0
    for NNcut in NN_cuts:
        for Discocut in Disco_cuts:
            A_sig_wgt=np.sum(sig_df[(sig_df[ABCD_name]>=NNcut) & (sig_df[DisCo_name]>=Discocut)][weight_name])
            A_bkg_wgt=np.sum(bkg_df[(bkg_df[ABCD_name]>=NNcut) & (bkg_df[DisCo_name]>=Discocut)][weight_name])
            A=A_sig_wgt+A_bkg_wgt
            B_sig_wgt=np.sum(sig_df[(sig_df[ABCD_name]>=NNcut) & (sig_df[DisCo_name]<Discocut)][weight_name])
            B_bkg_wgt=np.sum(bkg_df[(bkg_df[ABCD_name]>=NNcut) & (bkg_df[DisCo_name]<Discocut)][weight_name])
            B=B_sig_wgt+B_bkg_wgt
            C_sig_wgt=np.sum(sig_df[(sig_df[ABCD_name]<NNcut) & (sig_df[DisCo_name]>=Discocut)][weight_name])
            C_bkg_wgt=np.sum(bkg_df[(bkg_df[ABCD_name]<NNcut) & (bkg_df[DisCo_name]>=Discocut)][weight_name])
            C=C_sig_wgt+C_bkg_wgt
            D_sig_wgt=np.sum(sig_df[(sig_df[ABCD_name]<NNcut) & (sig_df[DisCo_name]<Discocut)][weight_name])
            D_bkg_wgt=np.sum(bkg_df[(bkg_df[ABCD_name]<NNcut) & (bkg_df[DisCo_name]<Discocut)][weight_name])
            D=D_sig_wgt+D_bkg_wgt
            result.loc[i]=[NNcut, Discocut, A,  B, C, D, B*C/D, A_bkg_wgt, A_sig_wgt]
            decorrelation.loc[i]=[NNcut, Discocut, B_bkg_wgt*C_bkg_wgt/D_bkg_wgt, A_bkg_wgt, B_bkg_wgt*C_bkg_wgt/(D_bkg_wgt*A_bkg_wgt)]
            i+=1
    return (result, decorrelation)


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

# return a list of working points for brute Scan conditions
def bruteCond(sel_dict, bins=20):
    name_list=sel_dict['name']
    cuts_list=sel_dict['range']
    # make working points
    wk_cuts_list=[]
    for i, name in enumerate(name_list):
        cut=cuts_list[i]
        cut_list=np.linspace(*cut, bins).tolist()
        wk_cuts_list.append(cut_list)
    wk_pts=list(it.product(*wk_cuts_list))
    return wk_pts

# apply scans on various cuts to select signal region A
def bruteScan(data_df, sel_dict, bins=20, weight_name='xsec_sf', truth_name='is_signal'):
    sig_df=data_df[data_df[truth_name]==1]
    bkg_df=data_df[data_df[truth_name]==0]
    sig_tot=sig_df[weight_name].sum()
    scan_df=pd.DataFrame(columns=['Cuts', 'A', 'SIG(A)', 'SIG(A)/SQRT(BKG(A))', 'SIG(A)/SIG(TOT)'])
    name_list=sel_dict['name']
    cuts_list=sel_dict['range']
    # make working points
    wk_cuts_list=[]
    for i, name in enumerate(name_list):
        cut=cuts_list[i]
        cut_list=np.linspace(*cut, bins).tolist()
        wk_cuts_list.append(cut_list)
    wk_pts=list(it.product(*wk_cuts_list))
    # get selection string for each working point
    for i in tqdm(range(len(wk_pts))):
        wk=wk_pts[i]
        string=''
        for j, name in enumerate(name_list):
            string+=f"{name} > {wk[j]}"
            if j!=len(name_list)-1: string+=' and '
        # evaluate
        bkg_A=bkg_df[bkg_df.eval(string)][weight_name].sum()
        sig_A=sig_df[sig_df.eval(string)][weight_name].sum()
        scan_df.loc[i]=[string, bkg_A+sig_A, sig_A, sig_A/np.sqrt(bkg_A), sig_A/sig_tot]
    return scan_df.sort_values(by=['SIG(A)/SQRT(BKG(A))'], ascending=False,ignore_index=True)

# plot Histogram of inferred score of signal versus background
def getInferScoreHist(infer_name, data_frame, truth_name, weight_name, bins_no, x_range, signal_scale):
    bins_edges=np.linspace(*x_range, bins_no)
    bkg_df=data_frame[data_frame[truth_name]==0]
    sig_df=data_frame[data_frame[truth_name]==1]
    
    total_bkg_wgt=round(np.sum(bkg_df[weight_name]),4)
    total_sig_wgt=round(np.sum(sig_df[weight_name]),4)
    bkg_hist=Hist1D(bkg_df[infer_name], bins=bins_edges, weights=bkg_df[weight_name], label="bkg"+" ("+str(total_bkg_wgt)+")")
    sig_hist=Hist1D(sig_df[infer_name], bins=bins_edges, weights=sig_df[weight_name]*signal_scale, label=fr"sig$\times${signal_scale} ({total_sig_wgt})")
    
    fig, axes=plt.subplots()
    hep.cms.label("Preliminary", data=True, lumi=138, loc=0, ax=axes)
    axes.set_xlabel("DisCo NN Score")
    axes.set_ylabel("count")
    bkg_hist.plot(ax=axes,alpha=1)
    sig_hist.plot(ax=axes,alpha=1)
    

# plot training history stored as dictionary
def plotHistory(history):
    train_loss = history['train_loss']
    train_bce = history['train_bce']
    train_disco = history['train_disco']
    test_disco = history['test_disco']
    test_loss = history['test_loss']
    test_bce = history['test_bce']

    

    fig, axes=plt.subplots()

    axes.plot(train_loss, alpha=0.8, linewidth=3, label='training loss')
    axes.plot(train_disco, alpha=0.8, linewidth=3, label='disco loss')
    axes.plot(train_bce, alpha=0.8, linewidth=3, label='bce loss')
    axes.plot()
    axes.legend()
    axes.set_title('Training Loss')

    fig, axes=plt.subplots()

    axes.plot(test_loss, alpha=0.8, linewidth=3, label='testing loss')
    axes.plot(test_disco, alpha=0.8, linewidth=3, label='disco loss')
    axes.plot(test_bce, alpha=0.8, linewidth=3, label='bce loss')
    axes.plot()
    axes.legend()
    axes.set_title('Testing Loss')

    fig, axes=plt.subplots()
    axes.plot(test_loss, linewidth=3, label='test loss')
    axes.plot(train_loss, linewidth=3, label='train loss')
    axes.legend()
    axes.set_title("Loss")
    
# plot distribution in DisCo Target for a certain ABCD cut
def plotDisCoWithCuts(ABCD_cut, NN_score, DisCo_name, df, truth_name, weight_name, x_range, bins_no):
    bins_edges=np.linspace(*x_range, bins_no)
    
    # make cut on ABCD score
    data_frame=df[df[NN_score]>=ABCD_cut]
    
    bkg_df=data_frame[data_frame[truth_name]==0]
    sig_df=data_frame[data_frame[truth_name]==1]
    total_bkg_wgt=round(np.sum(bkg_df[weight_name]))
    total_sig_wgt=round(np.sum(sig_df[weight_name]))/1000
    
    df_list=[bkg_df[DisCo_name],sig_df[DisCo_name]]
    weight_list=[bkg_df[weight_name],sig_df[weight_name]/100]
    label_list=[fr'bkg({total_bkg_wgt})', fr'sig$\times$10({total_sig_wgt})']
    
    fig, axes=plt.subplots()
    hep.cms.label("Preliminary", data=True, lumi=138, loc=0, ax=axes)
    plt.hist(df_list, stacked=True, alpha=0.6, bins=bins_no, label=label_list, weights=weight_list)
    axes.legend()
    axes.set_xlabel(r"$\Delta \eta_{jj}$")
    axes.set_ylabel("count")