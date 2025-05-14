import matplotlib.pyplot as plt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import numpy as np
import scipy.stats as st
import seaborn as sns
import pandas as pd
import pickle

data = pd.read_csv('/media/data/ytye/influenza/all_features_predict_lineage_year.csv', sep=',', header=None)
array = data.to_numpy()

#earliest
data = pd.read_csv('/media/data/ytye/influenza/earliest_features_seqname_predict.csv', sep=',', header=None)
array1 = data.to_numpy()
X1 = array1[:,2:42]
y1 = array1[:,1].astype('int')
lables1 = array1[:,0].astype('int')

X_train = array[:,2:42]
y_train = array[:,1].astype('int')
lables = array[:,0].astype('int')
lineages = array[:,42]
years = array[:,43].astype('float')

#cds
plt.close()
beginIdx = 30
endIdx = beginIdx + 5
X_train1 = X_train[:,beginIdx:endIdx]
lda = LinearDiscriminantAnalysis(n_components=2)
ldaFit = lda.fit(X_train1, y_train)
X_1  = X[:,beginIdx:endIdx]
X_r2 = ldaFit.transform(X_1)
ldaFit.explained_variance_ratio_

sporadic_avian = X_r2[y == -1, 0:2]
sporadic_avian_labels = lables[y == -1];
sporadic_avian = sporadic_avian[sporadic_avian_labels < 0, 0:2];

persistent = X_r2[lables == 1, 0:2]

x_sporadic_avian = sporadic_avian[:,0]
y_sporadic_avian = sporadic_avian[:,1]
x_persistent = persistent[:,0]
y_persistent = persistent[:,1]

sns.kdeplot(x=x_sporadic_avian, y=y_sporadic_avian, fill=True, cmap='gist_yarg')
sns.kdeplot(x=x_persistent, y=y_persistent, fill=True, cmap='Reds')

X1_1 = X1[:,beginIdx:endIdx]
X1_r2 = ldaFit.transform(X1_1)
target_names=['Hu1','Sw2','Eq','Ca2','Hu2-4,Sw1,Sw3-5,Ca1','Sw2(1985)']
#plot all earliest first
#plt.scatter(X1_r2[lables1 != 19, 0], X1_r2[lables1 != 19, 1], alpha=1, color="#00FFFF", label=target_names[4], s=30)
plt.plot(X1_r2[lables1 != 18, 0], X1_r2[lables1 != 18, 1], color="black", marker='o', markerfacecolor="#00FFFF", markersize=5, markeredgewidth=0.25, linewidth=0)

#Hu1
#plt.scatter(X1_r2[lables1 == 10, 0], X1_r2[lables1 == 10, 1], alpha=1, color="#FF0000", label=target_names[0], s=30)
plt.plot(X1_r2[lables1 == 10, 0], X1_r2[lables1 == 10, 1], color="black", marker='o', markerfacecolor="#FF0000", markersize=5, markeredgewidth=0.25, linewidth=0)

#Sw2
#plt.scatter(X1_r2[lables1 == 18, 0], X1_r2[lables1 == 18, 1], alpha=1, color="#0000FF", label=target_names[1], s=30)
plt.scatter(X1_r2[lables1 == 18, 0], X1_r2[lables1 == 18, 1], alpha=1, color="#0000FF", label=target_names[1], s=20)

#Eq
#plt.scatter(X1_r2[lables1 == 16, 0], X1_r2[lables1 == 16, 1], alpha=1, color="#00FF00", label=target_names[2], s=30)
plt.plot(X1_r2[lables1 == 16, 0], X1_r2[lables1 == 16, 1], color="black", marker='o', markerfacecolor="#00FF00", markersize=5, markeredgewidth=0.25, linewidth=0)

#Ca2
#plt.scatter(X1_r2[lables1 == 15, 0], X1_r2[lables1 == 15, 1], alpha=1, color="#FFFF00", label=target_names[3], s=30)
plt.plot(X1_r2[lables1 == 15, 0], X1_r2[lables1 == 15, 1], color="black", marker='o', markerfacecolor="#FFFF00", markersize=5, markeredgewidth=0.25, linewidth=0)

#Sw21985
#plt.scatter(X1_r2[lables1 == 19, 0], X1_r2[lables1 == 19, 1], alpha=1, color="navy", label=target_names[5], s=30)
plt.plot(X1_r2[lables1 == 19, 0], X1_r2[lables1 == 19, 1], color="black", marker='o', markerfacecolor="navy", markersize=5, markeredgewidth=0.25, linewidth=0)

plt.xlim(-4.5, 4.5)
plt.ylim(-4, 4)
plt.savefig("/media/data/ytye/influenza/16weight/allfeature_LDA_avianPDF_persistentPDF_earliest_PB1.pdf")

#Sw2 between 1979 to 1985
data = pd.read_csv('/container_data/ytye/influenza/data/work_N/before20220322_all/contents_exact_nofilterN/svm/earliest_Sw2_features_seqname_predict.csv', sep=',', header=None)
array2 = data.to_numpy()
X2 = array2[:,2:42]
X2_r2 = ldaFit.transform(X2)
#corn flower blue 
plt.plot(
    X2_r2[:, 0], X2_r2[:, 1], 'o-', alpha=0.5, color="#6495ED", label='Sw2', linewidth=0.15, markersize=2,
)

HA
>>> X1_r2
array([[-0.05560784, -0.16434227],
       [-1.29656687, -0.23603841],
       [-0.89476115, -0.06324991],
       [ 1.16540477,  0.48699716],
       [ 2.15158523, -0.92682431],
       [ 2.50809053, -1.58452055],
       [ 0.67580811,  0.16352936],
       [-0.40243148, -0.72117788],
       [ 0.85512618, -1.06297634],
       [ 0.81631231, -1.16054832],
       [ 2.04659261, -0.35964894],
       [ 0.93705815,  0.40966844],
       [-0.56910057, -0.91558323]])
>>> ldaFit.explained_variance_ratio_
array([9.99796831e-01, 2.03168853e-04]) 
      
NA
       X1_r2
array([[ 1.55606825, -0.49745469],
       [-0.46960034,  1.62686881],
       [ 0.70504278,  1.83574161],
       [ 1.82938584, -0.22888314],
       [ 3.37675256, -1.73450693],
       [ 1.1616247 ,  0.89169018],
       [ 1.14624206,  2.10827308],
       [ 2.26142592, -0.34733426],
       [-0.60948926,  1.18690331],
       [ 0.53740861,  1.20701873],
       [ 0.33837602, -0.15749057],
       [ 2.30374693, -2.04211276],
       [ 1.06187563, -0.42073032]])
>>> ldaFit.explained_variance_ratio_
array([9.99775656e-01, 2.2434134e-04])

NP
>>> X1_r2
array([[-1.05339964, -1.21829567],
       [ 0.53700039,  1.06793769],
       [ 1.43235098,  0.64237193],
       [ 2.69226542, -0.92059234],
       [ 2.60238975, -0.7644263 ],
       [-0.63396529,  0.79604333],
       [-0.44521827,  0.47488287],
       [-0.13354472, -0.00666481],
       [-1.32374852,  1.03950911],
       [-0.87020793,  0.37232889],
       [ 3.30050998, -0.11456395],
       [ 2.55921432, -1.09134983],
       [ 3.87307939, -0.41769299]])
  ldaFit.explained_variance_ratio_
array([9.99991711e-01, 8.28872006e-06])
     
 PA
 >>> X1_r2
array([[ 1.05593921e+00, -6.73567842e-01],
       [ 1.75421877e+00,  1.01532415e+00],
       [ 1.04242450e+00,  4.43497898e-01],
       [ 5.54383975e-01, -6.78807325e-01],
       [ 2.88637378e+00,  1.82817743e+00],
       [-1.87550098e-01, -2.06427566e+00],
       [ 2.35994397e-01, -1.58880325e-01],
       [ 1.77093752e+00, -7.88590521e-01],
       [ 5.61392707e-01, -4.28192506e-01],
       [ 3.87070482e-01, -2.12990983e+00],
       [ 5.01401862e-01,  3.04243938e-01],
       [ 5.54241416e-01, -8.72507921e-02],
       [ 1.83744713e-03, -1.89423905e-01]])
  ldaFit.explained_variance_ratio_
array([9.99633202e-01, 3.66798380e-04])

     
PB1
 X1_r2
array([[-0.51376178,  1.52603786], Hu1
       [-1.95604917, -0.23353899], Hu2
       [-0.27081033,  0.20355217], Hu3
       [ 2.55716243, -0.64089918], Hu4
       [ 1.03360564,  0.68960664], Ca1
       [ 1.40687947,  2.06046405], Ca2
       [-4.22092491,  2.20973057], Eq
       [-0.5180163 ,  0.72335146], Sw1
       [-0.14886093,  0.30224125], Sw2
       [-0.10283263,  0.08807885], Sw21985
       [ 2.75831218,  1.01682921], Sw3
       [ 2.63353352,  0.00918827], Sw4
       [ 2.63234808,  1.03100878]]) Sw5
>>> ldaFit.explained_variance_ratio_
array([9.99969144e-01, 3.08560159e-05])

PB2
>>> X1_r2
array([[-1.11954158,  1.81445471],
       [-0.00525513,  0.83558814],
       [ 1.55669336, -0.51033187],
       [ 0.87883406,  0.30094765],
       [ 4.20044811, -1.17046098],
       [-1.25624862,  0.23861847],
       [ 1.05543605, -2.3620802 ],
       [-1.15311656,  1.70028277],
       [ 0.46166936,  0.54622384],
       [ 0.6259723 ,  1.51307906],
       [-0.98003627,  0.9511546 ],
       [-1.05606379,  0.53759624],
       [ 0.028295  ,  1.01836236]])
    ldaFit.explained_variance_ratio_
array([9.99975981e-01, 2.40189707e-05])
    
M1
>>> X1_r2
array([[-0.72511604,  0.01576785],
       [ 1.0423941 , -1.33415069],
       [ 1.66612758, -1.69891872],
       [ 1.26074648,  0.80656834],
       [ 2.11371046,  1.77939756],
       [-0.88559972,  0.81788898],
       [ 0.32306189,  1.19050083],
       [ 1.09608809,  2.55379146],
       [-0.74112587,  0.12759869],
       [-0.0893923 , -0.11781294],
       [ 0.79930315,  0.03307942],
       [ 1.95509548,  1.01263349],
       [ 1.30502962,  0.74246847]])
ldaFit.explained_variance_ratio_
array([9.99975759e-01, 2.42406768e-05])

NS1
 X1_r2
array([[-0.7959197 ,  0.77631142],
       [ 1.72376495, -0.12000536],
       [ 1.44325478,  0.13209316],
       [ 1.9637905 ,  0.90362868],
       [ 0.85218867, -0.01878745],
       [ 0.41604065, -0.10673916],
       [-1.29569618,  0.63111697],
       [ 1.40060365, -1.25627526],
       [-1.20340025,  0.58972001],
       [ 0.51242039, -0.17853719],
       [ 2.11839353,  1.11709089],
       [ 2.11895545,  1.29168493],
       [ 1.954314  ,  1.35947944]])
array([9.99900953e-01, 9.90466797e-05])