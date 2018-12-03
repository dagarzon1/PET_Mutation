#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt


# In[3]:


AA3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W','ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
AA1to3 = {}
for i in range(len(list(AA3to1.keys()))):
    AA1to3[list(AA3to1.values())[i]] = list(AA3to1.keys())[i]


# In[4]:


AA = 'SSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ'
Am = list(AA)
Am3 = []
for i in range(len(Am)):
    Am3.append(AA1to3[Am[i]])


# In[5]:


CodonDict={'ATT':'I',   'ATC':'I',  'ATA':'I',  'CTT':'L',  'CTC':'L',  
'CTA':'L',  'CTG':'L',  'TTA':'L',  'TTG':'L',  'GTT':'V',  'GTC':'V',  
'GTA':'V',  'GTG':'V',  'TTT':'F',  'TTC':'F',  'ATG':'M',  'TGT':'C',  
'TGC':'C',  'GCT':'A',  'GCC':'A',  'GCA':'A',  'GCG':'A',  'GGT':'G',  
'GGC':'G',  'GGA':'G',  'GGG':'G',  'CCT':'P',  'CCC':'P',  'CCA':'P',  
'CCG':'P',  'ACT':'T',  'ACC':'T',  'ACA':'T',  'ACG':'T',  'TCT':'S',  
'TCC':'S',  'TCA':'S',  'TCG':'S',  'AGT':'S',  'AGC':'S',  'TAT':'Y',  
'TAC':'Y',  'TGG':'W',  'CAA':'Q',  'CAG':'Q',  'AAT':'N',  'AAC':'N',  
'CAT':'H',  'CAC':'H',  'GAA':'E',  'GAG':'E',  'GAT':'D',  'GAC':'D',  
'AAA':'K',  'AAG':'K',  'CGT':'R',  'CGC':'R',  'CGA':'R',  'CGG':'R',  
'AGA':'R',  'AGG':'R',  'TAA':'X',  'TAG':'X',  'TGA':'X'}


# In[6]:


AAtoGEN = {}
for i in range(len(list(CodonDict.keys()))):
    AAtoGEN[list(CodonDict.values())[i]] = list(CodonDict.keys())[i]


# In[7]:


Gen = []
for i in range(len(Am)):
    Gen.append(AAtoGEN[Am[i]])


# In[8]:


Gen = np.array(Gen)


# In[9]:


def mutation(gen):
    g = ''.join(gen)
    idx = np.random.randint(len(g))
    a = ['A','T','G','C']
    i = np.random.randint(4)
    d = list(g)
    d[idx] = a[i]
    g = ''.join(d)
    p = []
    for i in range(0,len(g),3):
        p.append(CodonDict[g[i:i+3]])
    return ''.join(p), g


# Gillespi

# In[ ]:


for k  in range(1):
    n_0 = 50
    ge = Gen
    n_max = 10000
    t = [0]
    pob = [0 for i in range(n_0)]
    n_mu = [0]
    n_mu2 = [0]
    n_n = [0]
    t_mu = [0]
    j=0
    f=0
    fit= [f]
    i=1
    b = 0
    stop=False
    while(True):
        g = 2.0*len(pob)+0.01*f*pob.count(1)
        ga = len(pob)+0.01*b*pob.count(2)
        mu = g * 0.01 
        s = g + ga + mu
        t_new = -1/s * np.log(np.random.random())
        t.append(t[i-1]+t_new)

        alpha = np.random.random()
        if(g > alpha * s and len(pob)<=n_max):
            r = np.array([0,1,2])
            prob = np.array([0, 0.33 + 0.1 * f, 0.33 + 0.1 * b])
            prob[0] = 1 - prob[1] - prob[2]
            pob.append(np.random.choice(r,p=prob))
        if(ga > alpha * s):
            idx = np.random.randint(len(pob))
            pob.pop(idx)
        if(mu > alpha * s):
            p, ge1 = mutation(ge)
            if (p[64] == 'Y'or p[138] == 'M' or p[137] == 'S' ):
                ge = ge1
                pob.append(1)
                if (p[67] == 'A' or p[67] == 'V' or p[94] == 'V' or p[215] == 'F' or p[185] == 'F' or p[217] == 'P'):
                    f = 0.34 * np.random.random()  
                if (p[67] == 'A' and p[94] == 'V' and p[215] == 'F' and p[185] == 'F' and p[217] == 'P'):
                    f = 0.80 * np.random.random()
                    print(t[i])
                    stop = True
            else:
                b = -0.45 * np.random.random()
                pob.append(2)
        n_mu.append(pob.count(1))
        n_n.append(pob.count(0))
        n_mu2.append(pob.count(2))
        fit.append(f)
        if(t[i]>10 or stop):
            break
        if(i%1000 == 0):
            print(i,t[i],f)
        i+=1
        
    plt.plot(t,n_mu)
    
    print(AA[67], AA[94], AA[215], AA[185], AA[217])
    print(p[67], p[94], p[215], p[185], p[217])

plt.xlabel('time')
plt.ylabel('number of mutants')
plt.show()
plt.savefig("mutants23.pdf") 


# In[ ]:


p


# In[11]:


AA


# In[12]:


plt.plot(t,n_mu) 
plt.show()


# In[ ]:




