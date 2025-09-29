import matplotlib.pyplot as plt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import numpy as np
import scipy.stats as st
import seaborn as sns
import pandas as pd
import pickle

data = pd.read_csv('all_features_predict_lineage_year.csv', sep=',', header=None)
array = data.to_numpy()

X = array[:,2:42]
y = array[:,1].astype('int')
lables = array[:,0].astype('int')
lineages = array[:,43]
years = array[:,42].astype('float')

lda = LinearDiscriminantAnalysis(n_components=2)
filename = 'lda_model.sav'
#pickle.dump(ldaFit, open(filename, 'wb'))
##load the model from disk
ldaFit = pickle.load(open(filename, 'rb'))

X_r2 = ldaFit.transform(X)

#print all avian and sporadic, Fig.4c (a huge image)
#sporadic = X_r2[lables == -1, 0:2]
#avian = X_r2[lables == -2, 0:2]
#persistent = X_r2[lables == 1, 0:2]
#plt.scatter(avian[:, 0], avian[:, 1], alpha=0.8, color="#DCDCDC", label="avian", s=5)
#plt.scatter(sporadic[:, 0], sporadic[:, 1], alpha=0.8, color="#778899", label="sporadic", s=5)
#plt.scatter(persistent[:, 0], persistent[:, 1], alpha=0.6, color="#990000", label="persistent", s=5)
#plt.legend(loc = "best")
#plt.xlim(-5, 6)
#plt.ylim(-5, 5)
#plt.savefig("allfeature_LDA_avian_sporadic_persistent.pdf")

#change the used index from 2:42 by every 5 postions to plot the LDA for each CDS
#this can generate Extended Data Fig. 4

sporadic_avian = X_r2[y == -1, 0:2]
sporadic_avian_labels = lables[y == -1];
sporadic_avian = sporadic_avian[sporadic_avian_labels < 0, 0:2];
x_sporadic_avian = sporadic_avian[:,0]
y_sporadic_avian = sporadic_avian[:,1]
persistent = X_r2[lables == 1, 0:2]
x_persistent = persistent[:,0]
y_persistent = persistent[:,1]

plt.figure(figsize=(10,6))

sns.kdeplot(x=x_sporadic_avian, y=y_sporadic_avian, fill=True, cmap='gist_yarg')
sns.kdeplot(x=x_persistent, y=y_persistent, fill=True, thresh=0.025, cmap='Reds')

#Fig. 5c
#G1-G14
data = pd.read_csv('H5_202204_202309_features_predict_seqname_label_genotype.csv', sep=',', header=None)
array8 = data.to_numpy()
X8 = array8[:,2:42]
X8_r2 = ldaFit.transform(X8)
lables8 = array8[:,-1]

color_top14=["#FDB462", "#1B9E77", "#D95F02", "#7570B3","#66A61E" , "#E7298A", "#D9D9D9","#BEBADA","#FB8072", "#80B1D3","#8DD3C7" , "#B3DE69" ,"#FCCDE5","#B15928"]
genotype_top14 = ["G1", "G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12","G13","G14"]
for index, genoytppe in enumerate(genotype_top14):
    if index == 0:
        plt.plot(X8_r2[lables8==genoytppe, 0], X8_r2[lables8==genoytppe, 1], alpha=0.8, color="black", marker='o', markerfacecolor=color_top14[index], markersize=4, markeredgewidth=0.25, linewidth=0, label=genoytppe)
    else:
        plt.scatter(X8_r2[lables8==genoytppe, 0], X8_r2[lables8==genoytppe, 1], alpha=0.8, color=color_top14[index], label=genoytppe, s=15)

#cattle 
data = pd.read_csv('feature_exact_predict_host_NorthAmerica_2334b.csv', sep=',', header=None)
array3 = data.to_numpy()
X3 = array3[:,2:42]
X3_r2 = ldaFit.transform(X3)
lables3 = array3[:,-1]
seqname3 = array3[:,1]

#plt.scatter( X3_r2[lables3=='B3.6', 0], X3_r2[lables3=='B3.6', 1], alpha=0.8, color="#FF9999", label='B3.6',  s=15,)
#plt.scatter(X3_r2[lables3=='B3.12', 0], X3_r2[lables3=='B3.12', 1], alpha=0.8,color="#F08080", label='B3.12',  s=15,)
#plt.scatter( X3_r2[lables3=='C2.1', 0], X3_r2[lables3=='C2.1', 1], alpha=0.8, color="#FFCC99", label='C2.1',  s=15,)
plt.scatter(X3_r2[lables3=='B3.13', 0], X3_r2[lables3=='B3.13', 1], alpha=0.8, color="#FF9933", label='B3.13',  s=15)
#plt.scatter(X3_r2[lables3=='B3.7', 0], X3_r2[lables3=='B3.7', 1], alpha=0.8,color="#FF6666", label='B3.7',  s=15)

#2.3.2.1c
clade_2321c = "2.3.2.1c"
plt.scatter(X8_r2[lables8==clade_2321c, 0], X8_r2[lables8==clade_2321c, 1], alpha=0.8, color="#EEF88A", label=clade_2321c, s=15)

#HuH3N8
data = pd.read_csv('HuH3N8.feature_predict.csv', sep=',', header=None)
array7 = data.to_numpy()
X7 = array7[:,2:42]
X7_r2 = ldaFit.transform(X7)
plt.scatter(X7_r2[:, 0], X7_r2[:, 1], color="#3CB371", label='HuH3N8', s=15)

#earliest
data = pd.read_csv('earliest_features_seqname_predict.csv', sep=',', header=None)
array1 = data.to_numpy()
X1 = array1[:,2:42]
y1 = array1[:,1].astype('int')
lables1 = array1[:,0].astype('int')
X1_r2 = ldaFit.transform(X1)
target_names=['Hu1','Sw2','Eq','Ca2','Hu2-4,Sw1,Sw3-5,Ca1','Sw2(1985)']
#plot all earliest first
plt.plot(X1_r2[lables1 != 18, 0], X1_r2[lables1 != 18, 1], color="black", marker='o', markerfacecolor="#00FFFF", markersize=5, markeredgewidth=0.25, linewidth=0)
#Hu1
plt.plot(X1_r2[lables1 == 10, 0], X1_r2[lables1 == 10, 1], color="black", marker='o', markerfacecolor="#FF0000", markersize=5, markeredgewidth=0.25, linewidth=0)
#Sw2
plt.scatter(X1_r2[lables1 == 18, 0], X1_r2[lables1 == 18, 1], alpha=1, color="#0000FF", s=20)
#Eq
plt.plot(X1_r2[lables1 == 16, 0], X1_r2[lables1 == 16, 1], color="black", marker='o', markerfacecolor="#00FF00", markersize=5, markeredgewidth=0.25, linewidth=0)
#Ca2
plt.plot(X1_r2[lables1 == 15, 0], X1_r2[lables1 == 15, 1], color="black", marker='o', markerfacecolor="#FFFF00", markersize=5, markeredgewidth=0.25, linewidth=0)
#Sw21985
plt.plot(X1_r2[lables1 == 19, 0], X1_r2[lables1 == 19, 1], color="black", marker='o', markerfacecolor="navy", markersize=5, markeredgewidth=0.25, linewidth=0)

#Fig. 5a
#Sw2 between 1979 to 1985
#data = pd.read_csv('earliest_Sw2_features_seqname_predict.csv', sep=',', header=None)
#array2 = data.to_numpy()
#X2 = array2[:,2:42]
#X2_r2 = ldaFit.transform(X2)
#corn flower blue 
#plt.plot(
#    X2_r2[:, 0], X2_r2[:, 1], 'o-', alpha=0.5, color="#6495ED", label='Sw2', linewidth=0.15, markersize=2,
#)
#data = pd.read_csv('EA_ancestral_feature.csv', sep=',', header=None)
#array3 = data.to_numpy()
#X3 = array3[:,2:42]
#lables3 = array3[:,0].astype('int')
#X3_r2 = ldaFit.transform(X3)
#X3_r2_real = X3_r2[lables3 == 1, :]
#X3_r2_synthetic = X3_r2[lables3 == 2, :]
#plt.plot(X3_r2_real[:, 0], X3_r2_real[:, 1], 'o-', color="#004C99", label='Sw2_real', linewidth=0.15, markersize=2)
#plt.plot(X3_r2_synthetic[:, 0], X3_r2_synthetic[:, 1], 'o-', color="#0066CC", label='Sw2_ancestral', linewidth=0.15, markersize=2)
#plt.savefig("allfeature_LDA_avianPDFBySVM_persistentPDFBySVM_Sw2_ancestral.pdf")

plt.legend(ncol=8, loc = "best")
plt.xlim(-5, 6)
plt.ylim(-5, 5)
plt.savefig("allfeature_LDA_avian_sporadic_persistent_genotype_cattle_2321c_huH3N8.pdf")