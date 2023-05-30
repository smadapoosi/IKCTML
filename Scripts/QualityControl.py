#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scanpy as sc
from sklearn import svm
from sklearn import metrics
from sklearn.svm import SVC
import pandas as pd
import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sys import argv
import os



# In[ ]:


anndata = sc.read(argv[1])



def rename(sobj):
    sobj.obs['Original_Cell_Types'] = sobj.obs['Original_Cell_Types'].astype(str)
    
    celltypeslist = []
    bettername = []
    celltypeslist.append(['Glomerular_endothelium_Young','EC-AEA_Menon','Endothelial Cells - AEA & DVR _Lake','GC-EC_Menon','Descending_vasa_recta_Young', 'Endothelial Cells - glomerular capillaries_Lake','EC-PT_Menon','Ascending_vasa_recta_Young','Endothelial Cells - AVR_Lake','EC_Wu','Endothelial Cells (unassigned)_Lake'])
    bettername.append('Endothelium')
    celltypeslist.append(['cSMC/MC_Menon','Mesangial_Young','Vascular Smooth Muscle Cells and pericytes_Lake','Pericyte_Wu','Mesangial Cells_Lake'])
    bettername.append('Perivascular and Mesangium')
    celltypeslist.append(['Fibroblast_ureter_Young','FIB_Menon','Myofibroblast_Wu','Interstitium_Lake','Fibroblast_Wu'])
    bettername.append('Fibrobasts')
    celltypeslist.append(['PODO_Menon','Glomerular_epithelium_and_podocytes_Young','Podocytes_Lake'])
    bettername.append('Podocytes')
    celltypeslist.append(['Proximal Tubule Epithelial Cells (S2)_Lake','DTL_Menon','Unknown - Novel PT CFH+ Subpopulation (S2)_Lake','LOH (DL)_Wu','Decending thin Limb_Lake','PEC_Menon','Proximal Tubule Epithelial Cells (S3)_Lake', 'Glomerular parietal eptihelial_Liao'])
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
    celltypeslist.append(['Proximal_tubule_convoluted_and_straight_Young','Proximal Tubule Epithelial Cells (S1/S2)_Lake','PT_Wu','Proximal_tubule_convoluted_tubular_Young','PT_Menon','Proximal convoluted tubule_Liao','Proximal tubule_Liao','Proximal straight tubule_Liao','Proximal Tubule Epithelial Cells - Fibrinogen+ (S3)_Lake'])
    bettername.append('Proximal Tubule')
    
    anndata.obs['Leveled_Names']='blank'
    sobj.obs['Leveled_Names'] = sobj.obs['Leveled_Names'].astype(str)
    
    count = 0

    for celltype in celltypeslist:
        for names in celltype:
            sobj.obs['Leveled_Names'][sobj.obs['Original_Cell_Types'] == names] = bettername[count]
        count += 1
    
    sobj.obs['Original_Cell_Types'] = sobj.obs['Original_Cell_Types'].astype('category')
    sobj.obs['Leveled_Names'] = sobj.obs['Leveled_Names'].astype('category')
    
    return(sobj)

anndata = rename(anndata)
# In[ ]:


resultslist = {}
o = SVC()
o = CalibratedClassifierCV(o)

trainanndata = anndata
testanndata = anndata

trainX = trainanndata.obsm['X_pca'] 
testX = testanndata.obsm['X_pca']

trainY = trainanndata.obs['Leveled_Names']
testY = testanndata.obs['Leveled_Names']

o.fit(trainX, trainY)
predicted = o.predict(testX)
prob = np.max(o.predict_proba(testX), axis = 1)
unlabeled = np.where(prob < 0.6)
predicted[unlabeled] = 'Unknown'


# In[ ]:


resultslist['confusion'] = pd.crosstab(testY, predicted, rownames=['True'], colnames=['Predicted'])
resultslist['rejection'] = len(predicted[predicted == 'Unknown'])/len(predicted)


# In[ ]:


resultslist['confusion'].to_csv('data/QCSVMConfusion.csv')
percentunknown = resultslist['confusion']['Unknown']/resultslist['confusion'].sum(axis=0)
percentunknown.to_csv('confusion/UnknownpercentsvmQCconfusion.csv')


# In[ ]:


results = anndata.obs['Leveled_Names'] == predicted
results.to_csv('data/cellstouse.csv')

