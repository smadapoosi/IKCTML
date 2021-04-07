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


resultslist['confusion'].to_csv('QCSVMConfusion.csv')
percentunknown = resultslist['confusion']['Unknown']/resultslist['confusion'].sum(axis=0)
percentunknown.to_csv('confusion/UnknownpercentsvmQCconfusion.csv')


# In[ ]:


results = anndata.obs['Leveled_Names'] == predicted
results.to_csv('data/cellstouse.csv')

