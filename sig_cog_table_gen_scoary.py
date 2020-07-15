import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import csv

scoary = 'C:/Users/Stefan/Desktop/Mod_7_Predators/SCOARY'

for file in os.listdir(scoary+'/OUTPUT/'):
    if '.csv' in file:
        df = pd.read_csv(scoary+'/OUTPUT/'+file)
        keep = []
        pval = []
        assoc = []
        for x in df['Gene']:
            ind = df.index[df['Gene'] == x].tolist()
            p = df['Worst_pairwise_comp_p'][ind]
            #p = df['Empirical_p'][ind]
            if float(p) <= 0.05:
                keep.append(x)
                pval.append(float(p))
                pos = df['Number_pos_present_in'][ind]
                neg = df['Number_neg_present_in'][ind]
                if float(pos) >= float(neg):
                    assoc.append('+')
                elif float(neg) >= float(pos):
                    assoc.append('-')
                else:
                    print('ERROR in '+file)
                
        
        if len(keep) == 0:
            continue
        
        else:
            newfilename = file[:file.find('_')]+"_cogs_scoary"
            data = open(scoary+'/ASSOC/'+"{0}.csv".format(newfilename), "w", newline = '')
            writer = csv.writer(data)
            writer.writerow(keep)
            writer.writerow(pval)
            writer.writerow(assoc)
            data.close()


func_table = pd.read_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/COG_FUNCTIONS.csv' )
func_table.COG = func_table.COG.astype(str)

dfs = []
names = []
dir_loc = 'C:/Users/Stefan/Desktop/Mod_7_Predators/SCOARY/ASSOC'

for file in os.listdir(dir_loc):
    names.append(file[:file.find('_')])
    df = pd.read_csv(dir_loc+'/'+file, index_col = False, header = None)
    df = pd.DataFrame.transpose(df)
    df.columns = ['COG', 'P', 'ASSOCIATION']
    df = df.sort_values(by='COG')
    df.COG = df.COG.astype(str)
    df = df.merge(func_table, on='COG')
    dfs.append(df)

all_cogs = pd.read_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/COG/expression_table/cog_table.csv')
all_cogs = all_cogs.columns
all_cogs = all_cogs[1:]

removelist = []
for x in func_table['COG']:
    if x not in all_cogs:
        removelist.append(x)

all_cogs = func_table[~func_table['COG'].isin(removelist)]
all_cogs = all_cogs.reset_index(drop=True)


COG_functions = {'A':'RNA processing and modification',
'B':'Chromatin Structure and dynamics',
'C':'Energy production and conversion',
'D':'Cell cycle control and mitosis',
'E':'Amino Acid metabolism and transport',
'F':'Nucleotide metabolism and transport',
'G':'Carbohydrate metabolism and transport',
'H':'Coenzyme metabolism',
'I':'Lipid metabolism',
'J':'Translation',
'K':'Transcription',
'L':'Replication and repair',
'M':'Cell wall/membrane/envelope biogenesis',
'N':'Cell motility',
'O':'Post-translational modification,\nprotein turnover, chaperone functions',
'P':'Inorganic ion transport and metabolism',
'Q':'Secondary Structure',
'T':'Signal Transduction',
'U':'Intracellular trafficing and secretion',
'V':'Defense mechanisms',
'W':'Extracellular Stuctures',
'Y':'Nuclear structure',
'Z':'Cytoskeleton',
'R':'General Functional Prediction only',
'S':'Function Unknown',
'X':'Mobilome: prophages, transposons'}

def function_annotator(data):
    n=0
    for fun in data['func']:
        if len(fun)==1:
            data['func'][n] = COG_functions[fun]
        if len(fun)>=2:
            x=len(fun)
            new = []
            for y in range(0,x):
                if y == x-1:
                    new.append(COG_functions[fun[y]])
                else:
                    new.append(COG_functions[fun[y]]+'+ ')
            data['func'][n] = ''.join(new)
        n+=1

function_annotator(all_cogs)

for x in dfs:
    function_annotator(x)

'''
# Save COGs
n=0
for x in dfs:  
    f = names[n]
    n+=1
    x = x.sort_values(by='P', ascending = True)
    x.to_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/#Results/Analysis 1_ Lineage Effect/Significant COGs + Function/CSV/'+f+'.csv', index=False)
'''

# Compile into dataframe of proportionate counts

names_pos_neg = []
for name in names:
    names_pos_neg.append(name+'_pos')
    names_pos_neg.append(name+'_neg')

count_list = []
for df in dfs:
    
    count_pos = []
    count_neg = []
    m=0
    for x in df['func']:
        if '+' in df['ASSOCIATION'][m]:
            if '+' in x:
                l = x.split('+ ')
                for n in l:
                    count_pos.append(n)
            else:
                count_pos.append(x)
        elif '-' in df['ASSOCIATION'][m]:
            if '+' in x:
                l = x.split('+ ')
                for n in l:
                    count_neg.append(n)
            else:
                count_neg.append(x)        
        m+=1
    
    c_sum = len(count_pos) + len(count_neg)
    
    count_pos = np.unique(count_pos, return_counts=True)
    count_pos = pd.DataFrame(count_pos)
    count_pos = pd.DataFrame.transpose(count_pos)
    count_pos[1] = count_pos[1]*100/c_sum
    
    count_neg = np.unique(count_neg, return_counts=True)
    count_neg = pd.DataFrame(count_neg)
    count_neg = pd.DataFrame.transpose(count_neg)
    count_neg[1] = count_neg[1]*100/c_sum
    
    n=0
    for x in count_pos[1]:
        count_pos[1][n] = round(count_pos[1][n], 2)
        n+=1
    n=0
    for x in count_neg[1]:
        count_neg[1][n] = round(count_neg[1][n], 2)
        n+=1
    
    count_list.append(count_pos)
    count_list.append(count_neg)

all_count = []
for x in all_cogs['func']:
    if '+' in x:
        l = x.split('+ ')
        for n in l:
            all_count.append(n)
    else:
        all_count.append(x)
    
all_count = np.unique(all_count, return_counts=True)
all_count = pd.DataFrame(all_count)
all_count = pd.DataFrame.transpose(all_count)
all_sum = sum(all_count[1])
all_count[1] = all_count[1]*100/all_sum

n=0
for x in all_count[1]:
    all_count[1][n] = round(all_count[1][n], 2)
    n+=1

expected = all_count
expected = expected.iloc[::-1]
expected = expected.reset_index(drop=True)

sizelist = []
for x in expected[0]:
    sizelist.append(float(0))
    
for x in names_pos_neg:
    expected[x] = sizelist


for x in expected[0]:
    n=0
    for y in count_list:       
        if x in list(y[0]):                  
            ind_val = y[0][y[0]==x].index[0]
            ind_val_2 = expected[0][expected[0]==x].index[0]
            prop = y[1][ind_val]
            expected[names_pos_neg[n]][ind_val_2] = float(prop)
        n+=1

expected = expected.sort_values(by=1, ascending = False)

# BAR CHARTS

# OXYGEN TOLERANCE

barWidth = 0.2           
r1 = np.arange(len(expected[0]))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
r4 = [x + barWidth for x in r3]


figure(figsize = (30,10))
plt.bar(r1, expected[1], color='black', width=barWidth, edgecolor='white', label='Expected')
plt.bar(r2, expected['aerobic_neg'], color='#6bd3ff', width=barWidth, edgecolor='black', label='Aerobic Neg', hatch='//')
plt.bar(r2, expected['aerobic_pos'], color='#6bd3ff', width=barWidth, edgecolor='black', label='Aerobic Pos', bottom=  expected['aerobic_neg'])
plt.bar(r3, expected['facultative_neg'], color='#b187ff', width=barWidth, edgecolor='black', label='Facultative Neg', hatch='//')
plt.bar(r3, expected['facultative_pos'], color='#b187ff', width=barWidth, edgecolor='black', label='Facultative Pos', bottom = expected['facultative_neg'])
plt.bar(r4, expected['anaerobic_neg'], color='#3af0b1', width=barWidth, edgecolor='black', label='Anaerobic_Neg', hatch='//')
plt.bar(r4, expected['anaerobic_pos'], color='#3af0b1', width=barWidth, edgecolor='black', label='Anaerobic_Pos', bottom = expected['anaerobic_neg'])

# Add xticks on the middle of the group bars
plt.title('Lineage Independent COG Functions Associated with Oxygen Tolerance', fontdict = {'fontsize':30})
plt.xlabel('Function', fontdict = {'fontweight':'bold', 'fontsize':20})
plt.xticks([r + barWidth for r in range(len(expected[0]))], expected[0], rotation=70, fontsize=14)
plt.ylabel('Proportion of Significant COGs (%)', fontdict = {'fontweight':'bold', 'fontsize':20}) 
# Create legend & Show graphic
plt.legend(prop = {'size':20})

plt.savefig('C:/Users/Stefan/Desktop/Mod_7_Predators/#Results/Analysis 1_ Lineage Effect/Bar Charts/O2_Tol.png',  bbox_inches='tight')


# PREDATORY STRATEGY

figure(figsize = (30,10))
plt.bar(r1, expected[1], color='black', width=barWidth, edgecolor='white', label='Expected')
plt.bar(r2, expected['group_neg'], color='#4de350', width=barWidth, edgecolor='black', label='Group Attack Neg', hatch='//')
plt.bar(r2, expected['group_pos'], color='#4de350', width=barWidth, edgecolor='black', label='Group Attack Pos', bottom=  expected['group_neg'])
plt.bar(r3, expected['epibiotic_neg'], color='#45d4ff', width=barWidth, edgecolor='black', label='Epibiotic Neg', hatch='//')
plt.bar(r3, expected['epibiotic_pos'], color='#45d4ff', width=barWidth, edgecolor='black', label='Epibiotic Pos', bottom=  expected['epibiotic_neg'])
plt.bar(r4, expected['non_neg'], color='#bdbdbd', width=barWidth, edgecolor='black', label='Non-Predator Neg', hatch='//')
plt.bar(r4, expected['non_pos'], color='#bdbdbd', width=barWidth, edgecolor='black', label='Non-Predator Pos', bottom=  expected['non_neg'])
# Add xticks on the middle of the group bars
plt.title('Lineage Independent COG Functions Associated with Predatory Behaviour', fontdict = {'fontsize':30})
plt.xlabel('Function', fontdict = {'fontweight':'bold', 'fontsize':20})
plt.xticks([r + barWidth for r in range(len(expected[0]))], expected[0], rotation=70, fontsize=14)
plt.ylabel('Proportion of Significant COGs (%)', fontdict = {'fontweight':'bold', 'fontsize':20})
# Create legend & Show graphic
plt.legend(prop = {'size':20})

plt.savefig('C:/Users/Stefan/Desktop/Mod_7_Predators/#Results/Analysis 1_ Lineage Effect/Bar Charts/Pred_1.png',  bbox_inches='tight')


# ML: -------------------------------------------------------------------------
     
ml = 'C:/Users/Stefan/Desktop/Mod_7_Predators/#Results/Analysis 2_ Machine Learning/LR/Count/Features'

keep = []
pval = []
assoc = []

for file in os.listdir(ml):
    if '.csv' in file:
        df = pd.read_csv(ml+'/'+file, index_col=0)
        keep = []
        pval = []
        assoc = []
        n=0
        for x in df['0']:
            keep.append(df.index[n])
            pval.append(abs(float(x)))
            n+=1
            if float(x) > 0:
                assoc.append('+')
            elif float(x) < 0:
                assoc.append('-')
            else:
                print('ERROR in '+file+str(x))

        else:
            newfilename = file[:file.find('.')]+"_cogs"
            data = open(ml+'/COGS/'+"{0}.csv".format(newfilename), "w", newline = '')
            writer = csv.writer(data)
            writer.writerow(keep)
            writer.writerow(pval)
            writer.writerow(assoc)
            data.close()
            
func_table = pd.read_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/COG_FUNCTIONS.csv' )
func_table.COG = func_table.COG.astype(str)

dfs = []
names = []

for file in os.listdir(ml+'/COGS'):
    names.append(file[:file.find('_')])
    df = pd.read_csv(ml+'/COGS'+'/'+file, index_col = False, header = None)
    df = pd.DataFrame.transpose(df)
    df.columns = ['COG', 'COEFF', 'ASSOCIATION']
    df = df.sort_values(by='COG')
    df.COG = df.COG.astype(str)
    df = df.merge(func_table, on='COG')
    dfs.append(df)

all_cogs = pd.read_csv('C:/Users/Stefan/Desktop/Mod_7_Predators/COG/expression_table/cog_table.csv')
all_cogs = all_cogs.columns
all_cogs = all_cogs[1:]

removelist = []
for x in func_table['COG']:
    if x not in all_cogs:
        removelist.append(x)

all_cogs = func_table[~func_table['COG'].isin(removelist)]
all_cogs = all_cogs.reset_index(drop=True)

for x in dfs:
    function_annotator(x)

n=0
for x in dfs:  
    f = names[n]
    n+=1
    x = x.sort_values(by='COEFF', ascending = False)
    x.to_csv(ml+'/Functions/'+f+'.csv', index=False)

