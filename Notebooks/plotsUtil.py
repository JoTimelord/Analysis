import matplotlib.pyplot as plt
import torch
import pandas as pd
import numpy as np
import yahist
from yahist import Hist1D
import mplhep as hep
import itertools as it
from tqdm import tqdm
import os
import warnings
import openpyxl


# from https://github.com/aminnj/yahist/blob/master/yahist/utils.py#L133 
def clopper_pearson_error(passed, total, level=0.6827):
    """
    matching TEfficiency::ClopperPearson(),
    >>> ROOT.TEfficiency.ClopperPearson(total, passed, level, is_upper)
    """
    import scipy.stats

    alpha = 0.5 * (1.0 - level)
    low = scipy.stats.beta.ppf(alpha, passed, total - passed + 1)
    high = scipy.stats.beta.ppf(1 - alpha, passed + 1, total - passed)
    return low, high

def simplifyError(passed,total,level=0.6827):
    low,high=clopper_pearson_error(passed, total, level)
    err=high-passed
    return err

# style a dataframe table
def makePretty(styler,color_code):
    styler.hide(axis='index')
    styler.format(precision=3)
    css_indexes=f'background-color: {color_code}; color: white;'
    styler.applymap_index(lambda _: css_indexes, axis=1)
    return styler
    
# calculate a dataframe that stores the background and signal in ABCD
# and a dataframe for closure in background only
# and return a dataframe for region A decomposition
def ABCDClosure(all_df, sel_dict, export_name, truth_name="is_signal",process_name="process",beautify=True):
    ABCD_df=pd.DataFrame(columns=['Region','Cut','BKG','SIG','TOT'])
    
    weight_name=sel_dict["weight"]
    total=np.sum(all_df[weight_name])
    sig=np.sum(all_df[all_df[truth_name]==True][weight_name])
    bkg=np.sum(all_df[all_df[truth_name]==False][weight_name])
    ABCD_df.loc[0]=["Total","/",bkg,sig,total]
    
    # selection strings
    presel_strings=sel_dict["presel"]
    region_A=sel_dict["A"]
    region_B=sel_dict["B"]
    region_C=sel_dict["C"]
    region_D=sel_dict["D"]
    
    
    df=all_df[all_df.eval(presel_strings)]
    sig_df=df[df[truth_name]==1]
    bkg_df=df[df[truth_name]==0]
    total_sig=np.sum(sig_df[weight_name])
    total_bkg=round(np.sum(bkg_df[weight_name]))
    total_wgt=total_sig+total_bkg
    ABCD_df.loc[1]=["After Preselection",presel_strings,total_bkg,total_sig,total_wgt]
    
    A_sig_wgt=np.sum(sig_df[sig_df.eval(region_A)][weight_name])
    A_bkg_wgt=np.sum(bkg_df[bkg_df.eval(region_A)][weight_name])
    A=A_sig_wgt+A_bkg_wgt
    ABCD_df.loc[2]=["A",region_A,A_bkg_wgt,A_sig_wgt,A]
    A_bkg_df=bkg_df[bkg_df.eval(region_A)].loc[:,[process_name,weight_name]]
    A_grouped=A_bkg_df.groupby([process_name]).sum()
    
    B_sig_wgt=np.sum(sig_df[sig_df.eval(region_B)][weight_name])
    B_bkg_wgt=np.sum(bkg_df[bkg_df.eval(region_B)][weight_name])
    B=B_sig_wgt+B_bkg_wgt
    ABCD_df.loc[3]=["B",region_B,B_bkg_wgt,B_sig_wgt,B]
    
    C_sig_wgt=np.sum(sig_df[sig_df.eval(region_C)][weight_name])
    C_bkg_wgt=np.sum(bkg_df[bkg_df.eval(region_C)][weight_name])
    C=C_sig_wgt+C_bkg_wgt
    ABCD_df.loc[4]=["C",region_C,C_bkg_wgt,C_sig_wgt,C]
    
    D_sig_wgt=np.sum(sig_df[sig_df.eval(region_D)][weight_name])
    D_bkg_wgt=np.sum(bkg_df[bkg_df.eval(region_D)][weight_name])
    D=D_sig_wgt+D_bkg_wgt
    ABCD_df.loc[5]=["D",region_D,D_bkg_wgt,D_sig_wgt,D]
    
    predict_A_bkg=B*C/(D)
    
    print(f"The predicted background in region A is {predict_A_bkg}.")
    
    ABCD_df.to_csv(f"{export_name}.csv")
    if beautify==True: 
        ABCD_styler=ABCD_df.style.pipe(makePretty,'#133760')
        css_alt_rows = 'background-color: yellow; color: black;'
        ABCD_styler=ABCD_styler.applymap(lambda x: css_alt_rows, subset=(slice(2,2), "BKG"))
        ABCD_styler.to_excel(f"{export_name}.xlsx", engine='openpyxl')
        ABCD_styler.to_latex(f"{export_name}.tex",position='ht',position_float={"centering"},caption="ABCD ")
        return (ABCD_styler,A_grouped)
    else: return (ABCD_df,A_grouped)


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
def bruteScan(data_df, sel_dict, save_name, bins=20, weight_name='xsec_sf', truth_name='is_signal'):
    warnings.filterwarnings("ignore")
    sig_df=data_df[data_df[truth_name]==1]
    bkg_df=data_df[data_df[truth_name]==0]
    sig_tot=sig_df[weight_name].sum()
    scan_df=pd.DataFrame(columns=['Cuts', 'A', 'SIG(A)', 'SIG(A)/SQRT(BKG(A))', 'SIG(A)/SIG(TOT)', 'BKG(A)'])
    name_list=sel_dict['name']
    cuts_list=sel_dict['range']
    # make working points
    wk_cuts_list=[]
    for i, name in enumerate(name_list):
        cut=cuts_list[i]
        cut_list=np.round(np.linspace(*cut, bins), decimals=2).tolist()
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
        scan_df.loc[i]=[string, bkg_A+sig_A, sig_A, sig_A/np.sqrt(bkg_A), sig_A/sig_tot, bkg_A]
    sorted_scan=scan_df.sort_values(by=['SIG(A)/SQRT(BKG(A))'], ascending=False,ignore_index=True)
    sorted_scan.to_csv(f'{save_name}.csv')
    return sorted_scan

# plot Histogram of inferred score of signal versus background
def getInferScoreHist(infer_name, data_frame, truth_name, weight_name, bins_no, x_range, signal_scale, save_fig, dir_name,separate=False):
    bins_edges=np.linspace(*x_range, bins_no)
    bkg_df=data_frame[data_frame[truth_name]==0]
    sig_df=data_frame[data_frame[truth_name]==1]
    
    total_bkg_wgt=round(np.sum(bkg_df[weight_name]),4)
    total_sig_wgt=round(np.sum(sig_df[weight_name]),4)
    bkg_hist=Hist1D(bkg_df[infer_name], bins=bins_edges, weights=bkg_df[weight_name], label="bkg"+" ("+str(total_bkg_wgt)+")")
    sig_hist=Hist1D(sig_df[infer_name], bins=bins_edges, weights=sig_df[weight_name]*signal_scale, label=fr"sig$\times${signal_scale} ({total_sig_wgt})")
    
    fig, axes=plt.subplots(figsize=(10,10))
    hep.cms.label("Preliminary", data=True, lumi=138, loc=0, ax=axes)
    axes.set_xlabel("DisCo NN Score")
    axes.set_ylabel("count")
    bkg_hist.plot(ax=axes,alpha=1)
    sig_hist.plot(ax=axes,alpha=1)
    if save_fig==True: plt.savefig(f'{dir_name}/NN_hist.png')
    
# plot training history stored as dictionary
def plotHistory(history, save_fig, save_dir, epoch_step=5):
    train_loss = history['train_loss']
    train_bce = history['train_bce']
    train_disco = history['train_disco']
    test_disco = history['test_disco']
    test_loss = history['test_loss']
    test_bce = history['test_bce']
    
    xaxis=np.arange(0,len(train_loss))
    
    fig, axes=plt.subplots(figsize=(15,10))

    axes.plot(xaxis[::epoch_step], train_loss[::epoch_step],alpha=0.8, color='green', linewidth=2,label=r'$L_{training}$')
    axes.plot(xaxis[::epoch_step],train_disco[::epoch_step], alpha=0.8, color='green',ls='dashed',linewidth=3,
              label=r'$L_{Disco}$')
    axes.plot(xaxis[::epoch_step], train_bce[::epoch_step], alpha=0.8,color='green',ls='dotted',linewidth=3, label=r'$L_{BCE}$')
    
    axes.plot(xaxis[::epoch_step], test_loss[::epoch_step], alpha=0.8,color='darkorange',linewidth=2, label=r'$L_{testing}$')
    axes.plot(xaxis[::epoch_step], test_disco[::epoch_step], alpha=0.8,
              linewidth=3,color='darkorange',ls='dashed',label=r'$L_{Disco}$')
    axes.plot(xaxis[::epoch_step], test_bce[::epoch_step], alpha=0.8,
              linewidth=3,color='darkorange',linestyle='dotted',label=r'$L_{BCE}$')
    axes.plot()
    axes.legend(prop={"family":"URW Bookman","size":16})
    axes.tick_params(axis='both',labelsize=12)
    axes.set_title('Training Loss',fontname='Nimbus Roman',fontsize=24)
    if save_fig==True: plt.savefig(f"{save_dir}/history_loss.png")
    
# plot distribution in DisCo Target for background for a range of ABCD cuts
def plotDisCoWithCuts(ABCD_cuts, NN_score, DisCo_name, df, weight_name, dir_name, epoch_no, truth_name='is_signal', x_range=(0,9), bins_no=20, save_fig=True, title=None):
    plt.style.use(hep.style.CMS)
    plt.rcParams.update({"figure.facecolor":  (1,1,1,0)})
    
    bins_edges=np.linspace(*x_range, bins_no)
    
    fig, axes=plt.subplots()
    hep.cms.label("Preliminary", data=True, lumi=138, loc=0, ax=axes)
    bkg_df=df[df[truth_name]==0]
    
    no_cut_hist=Hist1D(bkg_df[DisCo_name],bins=bins_edges,weights=bkg_df[weight_name],label='no cut').normalize()
    no_cut_hist.plot(ax=axes,alpha=0.5,histtype='stepfilled')
    
    j=1
    # make cut on ABCD score
    for i, cut in enumerate(ABCD_cuts):
        mask=bkg_df[NN_score] >=cut
        data_frame=bkg_df[mask]
        bkg_hist=Hist1D(data_frame[DisCo_name], bins=bins_edges, weights=data_frame[weight_name], label=f"{cut}<=NN").normalize()
        bkg_hist.plot(ax=axes,alpha=1,errors=True,color=f'C{j}')
        bkg_hist.plot(ax=axes,alpha=1,color=f'C{j}',label=None)
        j+=1
        mask=bkg_df[NN_score] <cut
        data_frame=bkg_df[mask]
        bkg_hist=Hist1D(data_frame[DisCo_name], bins=bins_edges, weights=data_frame[weight_name], label=f"{cut}>NN").normalize()
        bkg_hist.plot(ax=axes,alpha=1,errors=True,color=f'C{j}')
        bkg_hist.plot(ax=axes,alpha=1,color=f'C{j}',label=None)
        j+=1

    axes.set_xlabel(r"$\Delta \eta_{jj}$")
    axes.set_ylabel("count")
    axes.set_xlim(left=0)
    axes.set_ylim(bottom=0)
    # axes.set_title(f"epoch={epoch_no}")
    
    if save_fig==True: plt.savefig(f'{dir_name}/NNDisCo_hist_epoch{epoch_no}.png')
    
