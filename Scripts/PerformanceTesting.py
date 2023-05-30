#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
from sklearn import svm
from sklearn import metrics
from sklearn.svm import SVC
import pandas as pd
import numpy as np
from sklearn.calibration import CalibratedClassifierCV
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from sys import argv
import os


# In[2]:


resultslist = {}
withoutnames = ['Without Lake', 'Without Liao', 'Without Menon', 'Without Wu', 'Without Young']


# In[3]:


def rename(sobj):
    sobj.obs['Original_Cell_Types'] = sobj.obs['Original_Cell_Types'].astype(str)
    
    celltypeslist = []
    bettername = []
    celltypeslist.append(['Glomerular_endothelium_Young','EC-AEA_Menon','Endothelial Cells - AEA & DVR _Lake','GC-EC_Menon','Descending_vasa_recta_Young', 'Endothelial Cells - glomerular capillaries_Lake','EC-PT_Menon','Ascending_vasa_recta_Young','Endothelial Cells - AVR_Lake','EC_Wu'])
    bettername.append('Endothelium')
    celltypeslist.append(['cSMC/MC_Menon','Mesangial_Young','Vascular Smooth Muscle Cells and pericytes_Lake','Pericyte_Wu','Mesangial Cells_Lake'])
    bettername.append('Perivascular and Mesangium')
    celltypeslist.append(['Fibroblast_ureter_Young','FIB_Menon','Myofibroblast_Wu','Interstitium_Lake','Fibroblast_Wu'])
    bettername.append('Fibrobasts')
    celltypeslist.append(['PODO_Menon','Glomerular_epithelium_and_podocytes_Young','Podocytes_Lake'])
    bettername.append('Podocytes')
    celltypeslist.append(['Proximal Tubule Epithelial Cells (S2)_Lake','DTL_Menon','Unknown - Novel PT CFH+ Subpopulation (S2)_Lake','LOH (DL)_Wu','Decending thin Limb_Lake','Proximal Tubule Epithelial Cells - Fibrinogen+ (S3)_Lake','PEC_Menon','Proximal Tubule Epithelial Cells (S3)_Lake', 'Glomerular parietal eptihelial_Liao', 'PT_Wu'])
    bettername.append('Parietal Epithelium, Late Proximal Tubule, Descending Thin Limb')
    celltypeslist.append(['Mono1_Wu','MYL-2_Menon','Immune Cells - Macrophages_Lake','Monocytes_Liao','MNP_Young','MYL-1_Menon','Mono2_Wu','B Cells_Wu'])
    bettername.append('Monocyte')
    celltypeslist.append(['NK_Young','NK_Menon','NK_proliferating_Young','NKT_Young','NK-T_Liao','T cells_Wu','CD8T_Young','Th_Young','T_Menon'])
    bettername.append('Natural Killer, T')
    celltypeslist.append(['B_Menon','B_Liao','B_Young','Plasmacytoid_Young','Plasma2_Wu','Plasma1_Wu'])
    bettername.append('B, Plasma, Plasmacytoid')
    celltypeslist.append(['Intercalated_Liao','Collecting_duct_type_B_Young','IC_Menon','Collecting Duct - Intercalated Cells Type A (medulla)_Lake','Collecting Duct - Intercalated Cells Type A (cortex)_Lake','IC-B_Menon','Collecting_duct_type_A_Young','Collecting Duct - Intercalated Cells Type B_Lake'])
    bettername.append('Intercalated')
    celltypeslist.append(['PC_Menon','CD_Wu','Distal_tubule_and_collecting_duct_principle_cells_Young','Principal_Liao','Collecting Duct - Principal Cells (medulla)_Lake','Collecting system - Principal Cells (cortex)_Lake'])
    bettername.append('Principal')
    celltypeslist.append(['Urothelium_superficial_Young','Urothelium_Young','Urothelium_intermediate_and_basal_Young'])
    bettername.append('Urothelium')
    celltypeslist.append(['Thick Ascending Limb_Lake','TAL_Menon','Distal tubule_Liao','Loop_of_Henle_Young','LOH (AL)_Wu','ATL_Menon','Thin ascending limb_Lake'])
    bettername.append('Ascending Loop of Henle')
    celltypeslist.append(['DCT_Menon','Distal Convoluted Tubule_Lake','Connecting Tubule_Lake','Connecting Tubule_Lake','CNT_Menon'])
    bettername.append('Distal Convoluted and Connecting Tubules')
    celltypeslist.append(['Mast_Young','Mast cells_Wu'])
    bettername.append('Mast')
    celltypeslist.append(['Neutrophil_Young'])
    bettername.append('Neutrophil')
    celltypeslist.append(['Proximal_tubule_convoluted_and_straight_Young','Proximal Tubule Epithelial Cells (S1/S2)_Lake','Proximal_tubule_convoluted_tubular_Young','PT_Menon','Proximal convoluted tubule_Liao','Proximal tubule_Liao','Proximal straight tubule_Liao'])
    bettername.append('Proximal Tubule')
    
    sobj.obs['Leveled_Names'] = sobj.obs['Leveled_Names'].astype(str)
    
    count = 0

    for celltype in celltypeslist:
        for names in celltype:
            sobj.obs['Leveled_Names'][sobj.obs['Original_Cell_Types'] == names] = bettername[count]
        count += 1
    
    sobj.obs['Original_Cell_Types'] = sobj.obs['Original_Cell_Types'].astype('category')
    sobj.obs['Leveled_Names'] = sobj.obs['Leveled_Names'].astype('category')
    
    return(sobj)


# In[4]:


def run(anndata, classifier):
    trainanndata = anndata[anndata.obs['Study'] != 'Single']
    testanndata = anndata[anndata.obs['Study'] == 'Single']

    trainX = trainanndata.obsm['X_pca'] 
    testX = testanndata.obsm['X_pca']

    trainY = trainanndata.obs['Leveled_Names']
    testY = testanndata.obs['Leveled_Names']

    classifier.fit(trainX, trainY)
    predicted = classifier.predict(testX)
    prob = np.max(classifier.predict_proba(testX), axis = 1)
    unlabeled = np.where(prob < 0.6)
    predicted[unlabeled] = 'Unknown'
    
    outx = pd.crosstab(testY, predicted, rownames=['True'], colnames=['Predicted'])
    outy = len(predicted[predicted == 'Unknown'])/len(predicted)
    
    return(outx)


# In[5]:


resultslist = [[],[],[],[],[]]
for i in range(5):
    anndata = sc.read(argv[i+1])
    anndata = rename(anndata)
    
    l = SVC()
    l = CalibratedClassifierCV(l)
    
    resultslist[i].append(run(anndata, l))
    
    m = RandomForestClassifier()
    m = CalibratedClassifierCV(m)

    resultslist[i].append(run(anndata, m))
    
    n = MLPClassifier()
    n = CalibratedClassifierCV(n)

    resultslist[i].append(run(anndata, n))

    o = KNeighborsClassifier()
    o = CalibratedClassifierCV(o)

    resultslist[i].append(run(anndata, o))

    p = XGBClassifier()
    p = CalibratedClassifierCV(p)

    resultslist[i].append(run(anndata, p))


if not os.path.exists('scores'):
    os.mkdir('scores')
if not os.path.exists('confusion/RF'):
    os.mkdir('confusion/RF')
if not os.path.exists('confusion/SVC'):
    os.mkdir('confusion/SVC')
if not os.path.exists('confusion/KN'):
    os.mkdir('confusion/KN')
if not os.path.exists('confusion/XGB'):
    os.mkdir('confusion/XGB')
if not os.path.exists('confusion/MLP'):
    os.mkdir('confusion/MLP')

def medianf1(confusionmatrix):
    data = confusionmatrix
    
    data0 = data / data.sum(0)
    
    data1 = data / data.sum(1)
    
    dataf1 = (2 * data0 * data1)/(data0 + data1)
    dataf1 = dataf1.fillna(0)
    
    f1slist = []
    for j in dataf1.columns:
        if j in dataf1[j].keys():
            f1slist.append(dataf1[j][dataf1[j].keys() == j])
    
    
    
    f1slist = np.array(f1slist)
    
    f1slist = f1slist[f1slist > 0]
    return(pd.Series(f1slist).median())


# In[7]:


resultslist1 = [[],[],[],[],[]]
count = 0
for q in resultslist:
    for u in q:
        resultslist1[count].append(medianf1(u))
    count += 1


# In[8]:


tobemapped = pd.DataFrame(resultslist1)


# In[24]:


def formatting(df):
    df.columns = ['Lake','Liao','Menon','Wu','Young']
    df.index = ['SVC','RF','MLP','KN','XGB']
    return(df)


# In[25]:


tobemapped = formatting(tobemapped)


# In[9]:


tobemapped.to_csv('scores/F1_scores.csv')


# In[12]:


def uk(confusionmatrix):
    data = confusionmatrix
    percents = data['Unknown'] / data.sum(1) 
    
    return(percents.median())


# In[13]:


unknownlist1 = [[],[],[],[],[]]
count = 0
for q in resultslist:
    for u in q:
        unknownlist1[count].append(uk(u))
    count += 1


# In[28]:


ukdf = pd.DataFrame(unknownlist1)
ukdf = formatting(ukdf)
ukdf.to_csv('confusion/ukdf.csv')


# In[17]:


# In[47]:


ukdf1 = ukdf * 100
cmap = sns.color_palette('coolwarm_r', as_cmap = True)
sns.set_theme(font_scale = 3)
plt.figure(figsize=(20,15))
percentunknownheatmap = sns.heatmap(ukdf1, annot=True, cmap = cmap, robust = True)
percentunknownheatmap = percentunknownheatmap.get_figure()
percentunknownheatmap.savefig("confusion/Percent_Unknown_heatmap.pdf")

# In[45]:


cmap = sns.color_palette('coolwarm', as_cmap = True)
sns.set_theme(font_scale = 3)
plt.figure(figsize=(20,15))
medianf1heatmap = sns.heatmap(tobemapped, annot = True, cmap=cmap, robust = True)
medianf1heatmap = medianf1heatmap.get_figure()
medianf1heatmap.savefig("scores/Median_F1_heatmap.pdf")

# In[21]:

datasetlist = ['Lake','Liao','Menon','Wu','Young']
algorithmlist = ['SVC','RF','MLP','KN','XGB']


def medianf1maker(confusionmatrix):
    data = confusionmatrix
    
    data0 = data / data.sum(0)
    
    data1 = data / data.sum(1)
    
    dataf1 = (2 * data0 * data1)/(data0 + data1)
    dataf1 = dataf1.fillna(0)
    
    return(dataf1)
    
avglist=[]
for i in range(5):
    mat = resultslist[i][4]
    avglist.append(medianf1maker(mat))

columns = []
lst = []
for mat in avglist:
    for col in mat.columns:
        if not (col in lst):
            columns.append(col)
        lst.append(col)
        
rows = []
lst2 = []
for mat in avglist:
    for row in mat.index:
        if not (row in lst2):
            rows.append(row)
        lst2.append(row)

template=pd.DataFrame(columns=columns, index = rows)

for i in template.columns:
  #  print(i)
    for j in template.index:
     #   print(j)
        values = []
        for mat in avglist:
            if (i in mat.columns) and (j in mat.index):
             #   print(mat.loc[j,i])
                values.append(mat.loc[j,i])
        
        if len(values) > 0:
            template.loc[j,i] = sum(values) / len(values)
        #    print(len(values))

            
template = template.fillna(0)
template = template.sort_index()
template = template.sort_index(axis=1)
cmap = sns.color_palette('coolwarm', as_cmap = True)
sns.set(font='Arial', font_scale=0.85)
celltypef1heatmap = sns.heatmap(template, cmap=cmap, robust = True)
celltypef1heatmap = celltypef1heatmap.get_figure()
celltypef1heatmap.savefig("scores/CellType_F1_heatmap.pdf")
r = [] 
for m in range(5):
    mat = resultslist[m][4]
    r.append(mat['Unknown']/mat.sum(1))
k = []
for i in r:
    for j in i.index:
        if j not in k:
            k.append(j)
l = {}
for i in k:
    ilist = []
    for j in r:
        if i in j.index:
              ilist.append(j[i])
    l[i]= sum(ilist) / len(ilist)
    
l = pd.DataFrame.from_dict(l, orient='index')
l = 100*l
l=l.sort_index()
cmap = sns.color_palette('coolwarm', as_cmap = True)
sns.set(font='Arial', font_scale=0.85)
ukbct = sns.heatmap(l, cmap=cmap, robust = True)
ukbct = ukbct.get_figure()
ukbct.savefig("scores/UnknownByCellType_F1_heatmap2.pdf")


