
# Load Packages: --------------------------------------------------------------

import shap
import pickle

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math
import statsmodels.api as sm
from statsmodels.discrete import discrete_model

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold

from sklearn.metrics import accuracy_score
from sklearn.metrics import log_loss, r2_score
from sklearn.metrics import plot_roc_curve
from sklearn import datasets, metrics, model_selection, svm

from skopt import BayesSearchCV
from skopt.space import Real, Integer
from skopt.utils import use_named_args

from statsmodels.stats.outliers_influence import variance_inflation_factor    



# Load Data: ------------------------------------------------------------------

df = pd.read_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/COG/expression_table/cog_table.csv')
meta = pd.read_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/COG/expression_table/cog_metadata.csv')

X = df.drop(columns = ['name'])
y = meta.drop(columns = ['name'])

conds = list(y.columns)
conds.remove('endobiotic')
conds.remove('epibiotic')
conds.remove('group attack')
conds.remove('non predator')

log_loss_score = []
imp_feats = []

ROC_X = []
ROC_Y = []
ROC_MOD = []


# Train Models: ---------------------------------------------------------------

for cond in conds:

    CONDITION = cond
    
    X_train, X_test, y_train, y_test = train_test_split(X, y[CONDITION], test_size=0.2, random_state=42, stratify=y[CONDITION])
    	
    # Use SKOPT for Bayesian Parameter Optimisation
    model = LogisticRegression(penalty = 'l1', solver='liblinear')    
    params = {'C':[0.000001, 0.000003, 0.00001, 0.00003, 0.0001, 0.0003,
                   0.001, 0.003, 0.01, 0.03, 0.1, 0.15, 0.2]}
    
    opt = BayesSearchCV(
        model,
        params,
        cv = StratifiedKFold(
            n_splits = 5,
            shuffle = True,
            #random_state = 42
        ),
        n_iter = 20,
        verbose = 0,
        refit = True,
        #random_state = 42,
        scoring = 'f1'
        )
        
    opt.fit(X_train, y_train)
    
    best = dict(opt.best_params_)
    
    # Build Best Model
    model = LogisticRegression(penalty = 'l1', solver = 'liblinear',
                               C = best['C'])
    
    model.fit(X_train, y_train)
    
    # Analyse Model
    y_pred = model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    
    # Calc Logloss
    lg_probas = model.predict_proba(X_test)
    log_loss_score.append(log_loss(y_test, lg_probas))
    lg_probas = model.predict_proba(X_train)
    log_loss_score.append(log_loss(y_train, lg_probas))
    
    # Find Feature Importances
    
    coef_df = pd.DataFrame(model.coef_)   
    coef_df = pd.DataFrame(coef_df)
    coef_df.columns = X.columns
    coef_df = pd.DataFrame.transpose(coef_df)
    coef_df = coef_df.sort_values(by=[0])
    coef_df = coef_df[~(coef_df == 0).any(axis=1)]
    
    imp_feats.append(coef_df)
    
    # Save info for ROC curve plots
    ROC_MOD.append(model)
    ROC_X.append(X_test)
    ROC_Y.append(y_test)


# Multicolinarity of Important Features:

def calculate_vif_(X, thresh=5.0):
    variables = list(range(X.shape[1]))
    dropped = True
    while dropped:
        dropped = False
        vif = [variance_inflation_factor(X.iloc[:, variables].values, ix)
               for ix in range(X.iloc[:, variables].shape[1])]

        maxloc = vif.index(max(vif))
        if max(vif) > thresh:
            print('dropping \'' + X.iloc[:, variables].columns[maxloc] +
                  '\' at index: ' + str(maxloc))
            del variables[maxloc]
            dropped = True

    print('Remaining variables:')
    print(X.columns[variables])
    return X.iloc[:, variables]

feats_cog = []    
for df in imp_feats:
    ind_df = list(df.index.values)
    cogs_x = list(X.columns)
    for q in ind_df:
        cogs_x.remove(q)
    new_df = X.drop(columns = cogs_x)
    feats_cog.append(new_df)
 
vif_outputs = []
for df in feats_cog:
    vifscore = calculate_vif_(df)
    vif_outputs.append(vifscore)

# Round 2:

log_loss_score = []
imp_feats = []

ROC_X = []
ROC_Y = []
ROC_MOD = []


# Train Models: ---------------------------------------------------------------
n=0
for cond in conds:
    X = vif_outputs[n]
    n+=1
    CONDITION = cond
    
    X_train, X_test, y_train, y_test = train_test_split(X, y[CONDITION], test_size=0.2, random_state=42, stratify=y[CONDITION])
    	
    # Use SKOPT for Bayesian Parameter Optimisation
    model = LogisticRegression(penalty = 'l1', solver='liblinear')    
    params = {'C':[0.000001, 0.000003, 0.00001, 0.00003, 0.0001, 0.0003,
                   0.001, 0.003, 0.01, 0.03, 0.1, 0.15, 0.2]}
    
    opt = BayesSearchCV(
        model,
        params,
        cv = StratifiedKFold(
            n_splits = 5,
            shuffle = True,
            #random_state = 42
        ),
        n_iter = 20,
        verbose = 0,
        refit = True,
        #random_state = 42,
        scoring = 'f1'
        )
        
    opt.fit(X_train, y_train)
    
    best = dict(opt.best_params_)
    
    # Build Best Model
    model = LogisticRegression(penalty = 'l1', solver = 'liblinear',
                               C = best['C'])
    
    model.fit(X_train, y_train)
    
    # Analyse Model
    y_pred = model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    
    # Calc Logloss
    lg_probas = model.predict_proba(X_test)
    log_loss_score.append(log_loss(y_test, lg_probas))
    lg_probas = model.predict_proba(X_train)
    log_loss_score.append(log_loss(y_train, lg_probas))
    
    # Find Feature Importances
    
    coef_df = pd.DataFrame(model.coef_)   
    coef_df = pd.DataFrame(coef_df)
    coef_df.columns = X.columns
    coef_df = pd.DataFrame.transpose(coef_df)
    coef_df = coef_df.sort_values(by=[0])
    coef_df = coef_df[~(coef_df == 0).any(axis=1)]
    
    imp_feats.append(coef_df)
    
    # Save info for ROC curve plots
    ROC_MOD.append(model)
    ROC_X.append(X_test)
    ROC_Y.append(y_test)

# Important Feature Tables: ---------------------------------------------------

n=0
for x in imp_feats:
    x.to_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/#Results/Analysis 2_ Machine Learning/LR/Count/Features/'+conds[n]+'.csv')
    n+=1

# Plot ROC Curves: ------------------------------------------------------------
    
x=0
for cond in conds:
    CONDITION = cond
    title = CONDITION.capitalize()
    
    ns_probs = [0 for _ in range(len(ROC_Y[x]))]
    # predict probabilities
    lr_probs = ROC_MOD[x].predict_proba(ROC_X[x])
    # keep probabilities for the positive outcome only
    lr_probs = lr_probs[:, 1]
    # calculate scores
    ns_auc = metrics.roc_auc_score(ROC_Y[x], ns_probs)
    lr_auc = metrics.roc_auc_score(ROC_Y[x], lr_probs)
    # calculate roc curves
    ns_fpr, ns_tpr, _ = metrics.roc_curve(ROC_Y[x], ns_probs)
    lr_fpr, lr_tpr, _ = metrics.roc_curve(ROC_Y[x], lr_probs)
    # plot the roc curve for the model
    plt.plot(ns_fpr, ns_tpr, linestyle='--', color='black')
    plt.plot(lr_fpr, lr_tpr, color='red', marker='.', label='AUC = %.3f' % (lr_auc))
    # axis labels
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title+ ' ROC Logistic Regression')
    # show the legend
    plt.legend()
    # show the plot
    plt.savefig('C:/Users/Stefan/Desktop/Mod_7_Predators/#Results/Analysis 2_ Machine Learning/LR/Count/ROC/'+CONDITION+'.png')
    plt.clf()
    x+=1





'''
# Plot the effect of C on logloss: --------------------------------------------
n=0
dfs = []
for cond in conds:
    X = vif_outputs[n]
    n+=1
    CONDITION = cond
    
    X_train, X_test, y_train, y_test = train_test_split(X, y[CONDITION], test_size=0.2, random_state=42, stratify=y[CONDITION])
    	
    # Use SKOPT for Bayesian Parameter Optimisation  
    params = [0.000001, 0.000003, 0.00001, 0.00003, 0.0001, 0.0003,
                   0.001, 0.003, 0.01, 0.03, 0.1, 0.15, 0.2, 0.25, 0.3,
                   0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8,
                   0.85, 0.9, 0.95, 1]
    
    log_loss_test = []
    log_loss_train = []
    
    for param in params:
            
        model = LogisticRegression(penalty = 'l1', solver = 'liblinear',
                                   C = param)
        
        model.fit(X_train, y_train)
        
        # Analyse Model
        y_pred = model.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)
        
        # Calc Logloss
        lg_probas = model.predict_proba(X_test)
        log_loss_test.append(log_loss(y_test, lg_probas))
        lg_probas = model.predict_proba(X_train)
        log_loss_train.append(log_loss(y_train, lg_probas))
        
        print(param)
    
    
    daf = pd.DataFrame({'Penalty':params,'Test':log_loss_test,'Train':log_loss_train})

    fig, ax = plt.subplots()
    ax.plot(daf['Penalty'], daf['Train'], label='Train')
    ax.plot(daf['Penalty'], daf['Test'], label='Test')
    ax.legend()
    plt.xlabel('Penalty')
    plt.ylabel('Logloss')
    plt.title('SKL-LR Logloss: '+CONDITION)
    plt.savefig('C:/Users/Stefan/Desktop/Mod_7_Predators/#Results/Analysis 2_ Machine Learning/LR/Count/Logloss/'+CONDITION+'.png')
'''  


