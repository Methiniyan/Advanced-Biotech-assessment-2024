#!/usr/bin/env python
# coding: utf-8

# In[48]:


import os


# In[49]:


import pandas as pd


# In[50]:


dt_list = "Downloads/biotech/"   #import and merge the data files.
dt_file_list = []

for files in os.listdir(dt_list):
    if files.endswith(".txt"):
        fp = os.path.join(dt_list, files)
        df = pd.read_csv(fp, sep='\t')
        dt_file_list.append(df)

    
dt_file_list = pd.concat(dt_file_list, ignore_index=True)


# In[51]:


print(dt_file_list)


# In[52]:


# mutation = mut, this will help to identify the mutation.
def mut(row):
    if row['WildType.Sequence'] == row['Mutant.Sequence']:
        return "no mutation"
    elif len(row['WildType.Sequence']) == len(row['Mutant.Sequence']):
        return "substitution mutation"
    elif len(row['WildType.Sequence']) < len(row['Mutant.Sequence']):
        return "insertion mutation"
    else:
        return "deletion mutation"


# In[53]:


dt_file_list['mt_type'] = df.apply(mut, axis=1)


# In[54]:


#find the differences between cellviability, mRNA expression and proteirn expression.
dt_file_list['CellViability.Mut.Rep_mean'] = dt_file_list[['CellViability.Mut.Rep1', 'CellViability.Mut.Rep2', 'CellViability.Mut.Rep3']].mean(axis=1)
dt_file_list['CellViability.WT.Rep_mean'] = dt_file_list[['CellViability.WT.Rep1', 'CellViability.WT.Rep2', 'CellViability.WT.Rep3']].mean(axis=1)
dt_file_list['CelViability.Mut_difference']= dt_file_list['CellViability.Mut.Rep_mean'] -  dt_file_list['CellViability.WT.Rep_mean'] 


# In[55]:


dt_file_list['mRNA.Expression.Mut.Rep_mean'] = dt_file_list[['mRNA.Expression.Mut.Rep1', 'mRNA.Expression.Mut.Rep2', 'mRNA.Expression.Mut.Rep3']].mean(axis=1)
dt_file_list['mRNA.Expression.WT.Rep_mean'] = dt_file_list[['mRNA.Expression.WT.Rep1', 'mRNA.Expression.WT.Rep2', 'mRNA.Expression.WT.Rep3']].mean(axis=1)
dt_file_list['mRNA.Expression_difference']= dt_file_list['mRNA.Expression.Mut.Rep_mean'] -  dt_file_list['mRNA.Expression.WT.Rep_mean']


# In[56]:


dt_file_list['Protein.Expression.Mut.Rep_mean'] = dt_file_list[['Protein.Expression.Mut.Rep1', 'Protein.Expression.Mut.Rep2', 'Protein.Expression.Mut.Rep3']].mean(axis=1)
dt_file_list['Protein.Expression.WT.Rep_mean'] = dt_file_list[['Protein.Expression.WT.Rep1', 'Protein.Expression.WT.Rep2', 'Protein.Expression.WT.Rep3']].mean(axis=1)
dt_file_list['Protein.Expression_difference']= dt_file_list['Protein.Expression.Mut.Rep_mean'] -  dt_file_list['Protein.Expression.WT.Rep_mean']


# In[57]:


print(dt_file_list)


# In[82]:


#Find the location of mutation. location=loc
def loc(row):
    wt_seq = row["WildType.Sequence"]
    mut_seq = row["Mutant.Sequence"]
    
    mut_loc = [i for i in range(min(len(wt_seq), len(mut_seq))) if wt_seq[i] != mut_seq[i]]
    
    if any(position < 1000 for position in mut_loc):
        return "promoter"
    elif any(position >= 1000 for position in mut_loc):
        return "CDS"
    else:
        return "no mutation"

dt_file_list["mut_position"] = dt_file_list.apply(loc, axis=1)


# In[83]:


# Find the top 5 genes

def CellViavilityeffect(row):
    if row['mRNA.Expression.Mut.Rep_mean'] > 0:
        return 'Increase'
    elif row['mRNA.Expression.Mut.Rep_mean'] < 0:
        return 'Decrease'
    else:
        return 'No Effect'
    
dt_file_list['CellViavilityeffect'] = dt_file_list.apply(CellViavilityeffect, axis=1)
sort_genes = dt_file_list.sort_values(by='CellViavilityeffect', ascending=False)
top5_genes = sort_genes.head(5)
print(top5_genes[['Gene', 'mt_type', 'CellViavilityeffect', 'mut_position', 'CellViability.Mut.Rep_mean']]) 


# In[84]:


import matplotlib.pyplot as plt
import seaborn as sbn

fig, axes = plt.subplots(3, 1, figsize=(5, 10))

sbn.boxplot(x='Gene', y='CellViability.Mut.Rep_mean', hue='mt_type', data=top5_genes, palette='viridis', ax=axes[0])
axes[0].set_ylabel('CellViability.Mut.Rep_mean')

sbn.scatterplot(x='Gene', y='CellViability.Mut.Rep_mean', hue='mt_type', data=top5_genes, palette='viridis', ax=axes[1])
axes[1].set_xlabel('Gene')
axes[1].set_ylabel('CellViability.Mut.Rep_mean')

sbn.kdeplot(data=top5_genes, x='CellViability.Mut.Rep_mean', hue= 'mt_type', ax=axes[2])
axes[2].set_xlabel('CellViability.Mut.Rep_mean')
axes[2].set_ylabel('Density')

plt.tight_layout()
plt.show()


# In[ ]:




