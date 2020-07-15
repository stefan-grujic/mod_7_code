
# Load Packages: --------------------------------------------------------------
import pickle
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math
import statsmodels.api as sm
from statsmodels.discrete import discrete_model
from numpy.polynomial.polynomial import polyfit
from scipy.ndimage.filters import gaussian_filter1d
from scipy.ndimage.filters import median_filter


from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold

from sklearn.metrics import accuracy_score
from sklearn.metrics import log_loss, r2_score
from sklearn.metrics import roc_auc_score
import statistics

from skopt import BayesSearchCV
from skopt.space import Real, Integer
from skopt.utils import use_named_args

from statsmodels.stats.outliers_influence import variance_inflation_factor    




df = pd.read_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/COG/expression_table/cog_table.csv')
meta = pd.read_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/COG/expression_table/cog_metadata.csv')
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
proportions = [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]

# Load Data: ------------------------------------------------------------------

X = df.drop(columns = ['name'])
y = meta.drop(columns = ['name'])

conds = list(y.columns)
conds.remove('endobiotic')
conds.remove('epibiotic')
conds.remove('group attack')
conds.remove('non predator')

for CONDITION in conds:
    auc_list = []
    avg_list = []
    for a in proportions:
        df = pd.read_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/COG/expression_table/cog_table.csv')
        meta = pd.read_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/COG/expression_table/cog_metadata.csv')
    
        df[CONDITION] = meta[CONDITION]
        
        ABCD = df.groupby(df[CONDITION]).apply(pd.DataFrame.sample, frac=a).reset_index(drop=True)
        
        X = ABCD.drop(columns = ['name', CONDITION])
        y = ABCD[CONDITION]
        
        imp_feats = []
        
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y)
        	
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
        
        # Find Feature Importances
        
        coef_df = pd.DataFrame(model.coef_)   
        coef_df = pd.DataFrame(coef_df)
        coef_df.columns = X.columns
        coef_df = pd.DataFrame.transpose(coef_df)
        coef_df = coef_df.sort_values(by=[0])
        coef_df = coef_df[~(coef_df == 0).any(axis=1)]
        
        imp_feats.append(coef_df)
        
        if len(imp_feats[0])>1:
            
            feats_cog = []    
            for daf in imp_feats:
                ind_df = list(daf.index.values)
                cogs_x = list(X.columns)
                for q in ind_df:
                    cogs_x.remove(q)
                new_df = X.drop(columns = cogs_x)
                feats_cog.append(new_df)
             
            vif_outputs = []
            for daf in feats_cog:
                vifscore = calculate_vif_(daf)
                vif_outputs.append(vifscore)
                
            X = vif_outputs[0]
            
            
            
            # Sensitivity Analysis: -------------------------------------------------------
        ROC_avg =[]
        for n in range(1,31):
    
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y)
            	
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
            
            ns_probs = [0 for _ in range(len(y_test))]
            lr_probs = model.predict_proba(X_test)
            lr_probs = lr_probs[:, 1]
            ROC_AUC = roc_auc_score(y_test, lr_probs)
            ROC_avg.append(ROC_AUC)
        
        auc_list.append(statistics.mean(ROC_avg))
        avg_list.append(ROC_avg)
        
    # PICKLE AVG LIST
'''     
    with open("C:/Users/Stefan/Desktop/Mod_7_Predators/#Results/Analysis 2_ Machine Learning/LR/Sens/avglist.txt", "wb") as fp:   #Pickling
        pickle.dump(avg_list, fp)
    
    with open("C:/Users/Stefan/Desktop/Mod_7_Predators/#Results/Analysis 2_ Machine Learning/LR/Sens/avglist.txt", "rb") as fp:   # Unpickling
        avg_list = pickle.load(fp)
'''
    
    # BOX PLOTS
    
    from scipy.interpolate import UnivariateSpline
    
    daf = pd.DataFrame({'Subsample':proportions,'ROC AUC':auc_list})
    
    medlist = []
    for x in avg_list:
        medlist.append(statistics.median(x))
    
    x_ax = daf['Subsample']
    s = UnivariateSpline(x_ax, medlist, s=0.1)
    xs = np.linspace(0.5,1,11)
    y_smooth = s(xs)
    
    box_p = pd.DataFrame(avg_list)
    box_p.index = proportions
    
    fig, ax = plt.subplots(figsize = (15,8))
    plt.ylim(0.5,1)
    plt.xlim(0.5,1.0)
    ax.plot(range(1,12), y_smooth, color = 'red')
    plt.boxplot(box_p, patch_artist=True, 
                medianprops=dict(color='orange'),
                boxprops=dict(facecolor='lightyellow', color='black', alpha = 0.5))
    ax.plot(range(1,12), y_smooth, color = 'red')
    plt.xticks(list(range(1,len(proportions)+1)), proportions)
    plt.xlabel('Subsample Size', fontdict = {'fontsize':14})
    plt.ylabel('ROC AUC', fontdict = {'fontsize':14})
    plt.title('Sensitivity Analysis: '+ CONDITION[0].upper()+CONDITION[1:], fontdict = {'fontsize':20})
    plt.savefig('C:/Users/Stefan/Desktop/Mod_7_Predators/#Results/Analysis 2_ Machine Learning/LR/Sens/'+CONDITION+'box.png')


