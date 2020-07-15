#========================================================================#
# IMPORT PACKAGES:                                                       #
#========================================================================#

import csv
import os
from shutil import copyfile
import pandas as pd
import numpy as np

#========================================================================#
# WELCOME INTERFACE:                                                     #
#========================================================================#

'''

Prompts user for locations of necessary files as well as output location.

'''

print('\n\nThis script converts Prokka output folders into filestructure required for COG based ML analysis.')
print('\nPlease ensure that your Prokka output files are named using your required sample identifier.\n Example: Staphylococcus_aureus.gff\n\n')

prokka_folder = input('Please enter the location of the directory in which your Prokka output folders are stored.\n\n Example: /home/user/project_folder/prokka_output\n>>> ', )
root = input('\n\nPlease enter the location of the directory you wish the output of this script to be sotred.\n\n Example: /home/user/project_folder\n>>> ', )

optional = input('\nDo have the a compatible metadata table prepared? (Y/N)\n (NOTE: A folder of fasta files will also be needed.)\n>>> ', )
metadata_loc = None

while optional.upper() != 'Y' and optional.upper() != 'N':
    optional = input('\n Invalid response, please enter Y or N (not case sensitive).\n>>> ', )

if optional.upper() == 'Y':
    metadata_loc = input('\n\nPlease enter the location of your Metadata file.\n\n Example: /home/user/project_folder/metadata_table.csv\n>>> ')
    if metadata_loc[-4:] != '.csv':
        print('\n\nError: Metadata must be provided in the form of a .csv file')
        exit()
    fasta_loc = input('\n\nPlease enter the location of the directory in which the Fasta files you used for Prokka are stored.\n\n Example: /home/user/project_folder/fasta_folder \n>>>', )
elif optional.upper() == 'N':
    print('\n\nWARNING: Directory will only be partially complete.')

#========================================================================#
# GFF MOVE:                                                              #
#========================================================================#

'''

Takes output of Prokka and copies GFFs to a new directory.

'''
print('\n\nMoving GFF files to new folder...')

gff_folder = root+'/GFF_FASTA'
os.mkdir(gff_folder)

for foldername in os.listdir(prokka_folder):
    for filename in os.listdir(prokka_folder+'/'+foldername):  
        if '.gff' in filename:
            copyfile(prokka_folder+'/'+foldername+'/'+filename, gff_folder+'/'+filename)

#========================================================================#
# GFF FASTA STRIP:                                                       #
#========================================================================#

'''

Takes GFFs and trims the fasta seq from the end.

'''
print('\n\nTrimming FASTA seqs from GFF files...')


output_loc = root+'/GFF'
os.mkdir(output_loc)

for filename in os.listdir(gff_folder):
	gff = open(gff_folder+'/'+filename, 'r')
	gff = gff.read()
	value = gff.find('##FASTA')
	gff_trimmed = gff[0:value]
	text_file = open(output_loc+'/'+filename, "w")
	text_file.write(gff_trimmed)
	text_file.close()

#========================================================================#
# GFF TO COG:                                                            #
#========================================================================#

'''
Recursively strips all COGs from each gff file in a folder,
then saves the list in a .csv format in target folder 
then looks through the destination folder and creates
a .csv file of all of the unique COGs found.

Takes these COG lists and creates a binary expression table.

Also drops columns with zero variance.

'''

print('\n\nConverting GFF to COG binary expression matrix...')

foldername = output_loc
destination = root+'/COGS'

os.mkdir(destination)

# Create a list of COGs for each .gff file

for filename in os.listdir(foldername):
    
    gff = open(foldername+'/'+filename, "r")
    gff= gff.read()
    
    coglist = []
    indexidentifier = 0
    endidentifier = 0
    
    while indexidentifier != -1:
        indexidentifier = gff.find("COG")
        gff = gff[indexidentifier:]
        endidentifier = gff.find(";")
        coglist.append(gff[:endidentifier])
        gff = gff[endidentifier:]
    
    strippedlist = []
    for x in coglist:
        strippedlist.append(x[4:])
        
    strippedlist = set(strippedlist)
    strippedlist = list(strippedlist)
    strippedlist.sort()
    
    for x in strippedlist:
        if "COG" not in x:
            strippedlist.remove(x)
    
    fileindex = filename.find('.')
    newfilename = filename[:fileindex]

    file = open(destination+'/'+"{0}.csv".format(newfilename), "w", newline = '')
    writer = csv.writer(file)
    writer.writerow(strippedlist)
    file.close()

# Create a list of COGs for all .gff files

coglistcomplete = []

for filename in os.listdir(destination+'/'):
    
   cog = open(destination+'/'+filename, "r")
   cog = cog.read()
   cog = cog.rstrip()
    
   coglist = list(cog.split(","))
    
   for x in coglist:
       coglistcomplete.append(x)
      
strippedlist = set(coglistcomplete)
strippedlist = list(strippedlist)
strippedlist.sort()

file = open(destination+'/all_cogs.csv', "w", newline = '')
writer = csv.writer(file)
writer.writerow(strippedlist)
file.close()

# Write a Binary Expression table for all COGs

cognames = open(destination+'/all_cogs.csv', "r")
cognames = cognames.read()
cognames = cognames.rstrip()
cognames = "name,"+ cognames

cognames = list(cognames.split(","))

os.mkdir(destination+'/expression_table')

file = open(destination+"/expression_table/cog_table.csv", "w", newline = '')
writer = csv.writer(file)
writer.writerow(cognames)

for filename in os.listdir(destination+'/'):
    
    if 'all_cogs' in filename:
        continue
    
    elif 'expression_table' in filename:
        continue
    
    else:
        cog = open(destination+'/'+filename, "r")
        cog = cog.read()
        cog = cog.rstrip()
        cog = list(cog.split(","))
        
        fileindex = filename.find('.')
        orgname = filename[:fileindex]
        
        cogbinary = []
        
        for x in cognames[1:]:
            if x in cog:
                cogbinary.append(1)
            else:
                cogbinary.append(0)

        cogbinary.insert(0, orgname)
        writer.writerow(cogbinary)

file.close()

# Sort csv by Name

table = pd.read_csv(destination+'/expression_table/cog_table.csv')
sorted_table = table.sort_values('name')

# Remove columns with zero variance

trimlist = []
for column in sorted_table:
    if column == 'name':
        continue
    elif sum(sorted_table[column]) == len(sorted_table.index):    
        trimlist.append(column)

sorted_table = sorted_table.drop(columns = trimlist)

sorted_table.to_csv(destination+'/expression_table/cog_table.csv', index = False)

# Exits if no metadata has been provided.

if optional.upper() == 'N':
    print('\n\n Processing complete, further processes require provision of a compatible metadata table.')
    exit()


#========================================================================#
# ONEHOTENCODE METADATA:                                                 #
#========================================================================#

'''
Takes compatible metadata format and produces a onehotencoded metadata table.
'''

print('\n\nOnehotencoding metadata...')

metatable = pd.read_csv(metadata_loc)

metadata = metatable[['Name', 'Oxygen Tolerance', 'Predation Strategy']]
metadata = metadata.sort_values('Name')

aerobe = []
anaerobe = []
facultative = []

for x in metadata['Oxygen Tolerance']:
    if x.lower() == 'aerobic':
        aerobe.append(1)
        anaerobe.append(0)
        facultative.append(0)
    elif x.lower() == 'anaerobic':
        aerobe.append(0)
        anaerobe.append(1)
        facultative.append(0)
    elif x.lower() == 'facultative':
        aerobe.append(0)
        anaerobe.append(0)
        facultative.append(1)
    else:
        print(x)

endobiotic = []
epibiotic = []
group = []
nonpred = []
pred = []

epi = ['epibiotic', 'lysis (diffusable)', 'undefined epibiotic', 
       'picket fence', 'lytic enzymes (close)', 'undefined epiobiotic']
end = ['periplasmic', 'endobiotic', 'cytoplasmic']
grp = ['group-attack', 'ixotrophy', 'myxobacteria-like', 'wolf-pack', 'secondary metabolites (far)']

for x in metadata['Predation Strategy']:
    if x.lower() in end:
        endobiotic.append(1)
        epibiotic.append(0)
        group.append(0)
        nonpred.append(0)
        pred.append(1)
    elif x.lower() in epi:
        endobiotic.append(0)
        epibiotic.append(1)
        group.append(0)
        nonpred.append(0)
        pred.append(1)
    elif x.lower() in grp:
        endobiotic.append(0)
        epibiotic.append(0)
        group.append(1)
        nonpred.append(0)
        pred.append(1)
    elif x.lower() == 'non-predator':
        endobiotic.append(0)
        epibiotic.append(0)
        group.append(0)
        nonpred.append(1)
        pred.append(0)
    else:
        print(x)

expression = pd.read_csv(destination+'/expression_table/cog_table.csv')

d = {'name': list(expression['name']), 'aerobic':aerobe, 'anaerobic':anaerobe, 
     'facultative':facultative, 'predator':pred, 'endobiotic':endobiotic, 
     'epibiotic':epibiotic, 'group attack':group, 'non predator':nonpred}

onehot = pd.DataFrame(data=d)

onehot.to_csv(destination+'/expression_table/cog_metadata.csv', index = False)

o_t = []
p_s =[]

for x in metadata['Oxygen Tolerance']:
    if x.lower() == 'anaerobic':
        o_t.append(0)
    elif x.lower() == 'facultative':
        o_t.append(1)
    elif x.lower() == 'aerobic':
        o_t.append(2)
    else:
        print(x)
        
for x in metadata['Predation Strategy']:
    if x.lower() == 'non-predator':
        p_s.append(0)
    elif x.lower() in grp:
        p_s.append(1)
    elif x.lower() in epi:
        p_s.append(2)
    elif x.lower() in end:
        p_s.append(3)
    else:
        print(x)
        
md = {'name': list(expression['name']), 'oxygen_tolerance':o_t,
      'predation_strategy':p_s}


multinom = pd.DataFrame(data = md)

multinom.to_csv(destination+'/expression_table/cog_multinomial_metadata.csv', index=False)

multitxt = 'Oxygen Tolerance:\n0 = Anaerobic\n1 = Facultative\n2 = Aerobic\n\nPredation Strategy:\n0 = Non Predator\n1 = Group Attack\n2 = Epibiotic\n3 = Endobiotic'
file = open(destination+'/expression_table/multinomial_info.txt', 'w')
file.write(multitxt)
file.close()


#========================================================================#
# LOCATE 16S rRNA AND GENERATE XSTRING:                                  #
#========================================================================#

'''
Searches for 16S rRNA in each gff and produces an Xstring compatible with
Xstring_to_newick.R 

'''

print('\n\nLocating 16S rRNA and generating Xstring...')

foldername = fasta_loc
destination = root+'/FASTA_TXT'
os.mkdir(destination)

for filename in os.listdir(foldername):
    filenametxt= filename[:-6]+'.txt'
    with open(foldername+'/'+filename, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(destination+'/'+filenametxt, 'w') as fout:
        fout.writelines(data[1:])
        
foldername = root+'/GFF'
destination = root+'/GFF_TSV'
os.mkdir(destination)

for filename in os.listdir(foldername):
    filenametsv= filename[:-4]+'.tsv'
    with open(foldername+'/'+filename, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(destination+'/'+filenametsv, 'w') as fout:
        fout.writelines(data[2:])
        
tsvname = root+'/GFF_TSV'
fastaname = root+'/FASTA_TXT'

RNA_dict = {}

for filename in os.listdir(tsvname):
    data = pd.read_csv(tsvname+'/'+filename, delimiter = '\t', header=None)
    f = [] 
    for x in data[8]:
        if '16S ribosomal RNA' in x:
            a = data.index[data[8]==x]
            b = int(data.iloc[a][3]-1)
            c = int(data.iloc[a][4]-1)
            f.append([b,c])
        RNA_dict[filename[:-4]] = f


RNA_code = {}

for filename in os.listdir(fastaname):
    fastafile = open(fastaname+'/'+filename,'r')
    fastafile = fastafile.read()
    fastafile = fastafile.rstrip()
    
    codeonly = []
    for x in fastafile:
        if x not in ['A','T','G','C']:
            continue
        else:
            codeonly.append(x)
        
    codelist = []
    for x in RNA_dict[filename[:-4]]:
        codelist.append(''.join(codeonly[x[0]:x[1]]))
    
    RNA_code[filename[:-4]] = codelist

mystring=[]

for x in RNA_code:
    mystring.append('>'+x.replace(' ','_')+'\n'+RNA_code[x][0])
    mystring.append('\n')
    
del mystring[-1]
mystring = ''.join(mystring)


os.mkdir(root+'/TREE')
os.mkdir(root+'/SCOARY')
os.mkdir(root+'/SCOARY/INPUT')
os.mkdir(root+'/SCOARY/OUTPUT')

text_file = open(root+'/TREE/16S_XString.fasta', "w")
text_file.write(mystring)
text_file.close()

text_file = open(root+'/SCOARY/INPUT16S_XString.fasta', "w")
text_file.write(mystring)
text_file.close()

#========================================================================#
# COG DATA TO SCOARY:                                                            #
#========================================================================#

'''
Produces a binary expression table in a scoary compatible format.

'''
print('\n\nProducting Scoary compatible binary expression matrix...')


meta = pd.read_csv(root+'/COGS/expression_table/cog_metadata.csv', index_col=None)
exp = pd.read_csv(root+'/COGS/expression_table/cog_table.csv', index_col='name')

exp_t = pd.DataFrame.transpose(exp)
exp_t = pd.DataFrame.reset_index(exp_t)

for x in exp_t.columns:
    exp_t = exp_t.rename(columns={x:x.replace(' ','_')})

for x in meta['name']:
    meta["name"].replace({x: x.replace(' ','_')}, inplace=True)

genes_col = list(exp_t['index'])

exp_t = exp_t.drop(columns=['index'])
# Change zero valuse into blanks
exp_t[exp_t.eq(0)] = ""

cols = ['Gene', 'Non-unique Gene name','Annotation','No. isolates','No. sequences','Avg sequences per isolate','Genome Fragment','Order within Fragment','Accessory Fragment', 'Accessory Order with Fragment','QC','Min group size nuc','Max group size nuc','Avg group size nuc']

blanks = []

for x in genes_col:
    blanks.append("")

blank_list = []
blank_list.append(genes_col)

for col in cols[1:]:
    blank_list.append(blanks)
    
df = pd.DataFrame(blank_list)
df = pd.DataFrame.transpose(df)
df.columns = cols
df['Gene'] = df['Gene'].astype(str)

joined = pd.concat([df,exp_t], axis=1)

meta = meta.rename(columns={'name':np.nan})

joined.to_csv(root+'/SCOARY/INPUT/cog_table.csv', index=False)
meta.to_csv(root+'/SCOARY/INPUT/cog_metadata.csv', index = False)

print('\n\nAll processes complete, please run Xstring_to_newick.R to complete next step.')
