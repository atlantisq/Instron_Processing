# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 22:44:42 2020

@author: atlan
"""

import os
import pandas
import matplotlib.pyplot as plt
import re
from scipy import integrate
import numpy as np

directory = os.getcwd()
print(directory)


pandas.set_option('display.max_rows', None)
pandas.set_option('display.max_columns', None)
pandas.set_option('display.width', None)
pandas.set_option('display.max_colwidth', -1)

#Tan_delta = pandas.DataFrame()
plt.figure(dpi=100)
color_map = {0:'dodgerblue', 5:'green', 15:'orange', 25:'red',26:'palegoldenrod', None:'dodgerblue', 10:'purple'}
counter = 0
data = pandas.DataFrame(columns=['filename', 'TFK', 'Toughness','Elong_at_break','UTensileStr','E'])
avgdata = pandas.DataFrame(columns=['TFK', 'Toughness', 'Stdev', 'Elong_at_break', 'Elong_std', 'UTensileStr', 'UStr_std','E','E_std'])



for filename in os.listdir(directory):
    tfk = None
    file = (os.path.join(directory, filename))
    print(file)
    if os.path.splitext(file)[-1] != '.csv':
        continue
    counter += 1
    
    match = re.search('TFK(\d+)', filename, re.IGNORECASE)
    if match:
        tfk = (int(match.group(1)))
        
    df = pandas.read_csv(file, names = ['x','stress'], usecols = [5,6])
    df = df.drop(0, axis = 0)
    df = df.dropna(axis=0, how = 'any')
    df = df.astype(float)
    #print(df)
    
    i = 1
    while df['stress'][i] < 0:
        i+=1
    
    default_x = df['x'][i]
    #print(default_x)
    
    df['x'] = df['x'] - default_x
    
    df['strain'] = df['x']/17
    df['stress'] = df['stress']/3
    #print (df)
    
    for i in range(1,len(df)):
        if df['stress'][i] - df['stress'][i+10] > 4:
            elong_at_break = df['strain'][i+10]
            df_truncated = df[0:i+10]
            #print (df_truncated)
            break
    
    if True:
        plt.plot(df['strain'],df['stress'], c=color_map[tfk])
    #break
    s = df_truncated.apply(lambda x: integrate.trapz(x,df_truncated['x']))
    #print(s)
    
    df_elastic = df.loc[(df['strain'] < 0.01)]
    df_elastic = df_elastic.loc[(df_elastic['strain'] > 0.005)]
    df_preyield = df.loc[(df['strain'] < 0.05)]
    
    #print(df_elastic)
    #plt.plot(df_elastic['strain'],df_elastic['stress'], c=color_map[tfk])
    slopes = df_elastic.apply(lambda x: np.polyfit(df_elastic['strain'], x, 1))
    #print(slopes['stress'][0])

    data = data.append({'filename' : filename, 'TFK': tfk, 'Toughness': s['stress'], 'Elong_at_break' : elong_at_break, 'UTensileStr': df_preyield['stress'].max(), 'E': slopes['stress'][0]}, ignore_index=True)
    
    
    
plt.show()
data = data.sort_values(by = 'TFK')
print(data)

tfk_list = data.TFK.unique()

avg = {}
stdev = {}
for tfk in tfk_list:
    data_subset = data.loc[data['TFK'] == tfk]
    mean = data_subset.mean()['Toughness']
    std = data_subset.std()['Toughness']
    
    elong_at_br = data_subset.mean()['Elong_at_break']
    elong_std = data_subset.std()['Elong_at_break']
    
    ulti_str = data_subset.mean()['UTensileStr']
    ultistr_std = data_subset.std()['UTensileStr']
    
    E = data_subset.mean()['E']
    E_std = data_subset.std()['E']
    
    #print(mean, std)
    avgdata = avgdata.append({'TFK': int(tfk), 'Toughness': mean,'Stdev': std, 'Elong_at_break': elong_at_br, 'Elong_std' : elong_std, 'UTensileStr' : ulti_str, 'UStr_std' : ultistr_std, 'E' : E, 'E_std': E_std}, ignore_index=True)
    avg[tfk] = mean
    stdev[tfk] = std
    
print(avgdata)


def plot_data(average, std):
    plt.figure(dpi=100)
    x_pos = [i for i, _ in enumerate(avgdata['TFK'])]
#print(x_pos)

    for i in x_pos:
        plt.bar(i, avgdata[average][i], color=color_map[avgdata['TFK'][i]], yerr=avgdata[std][i])
        #plt.xlabel("Energy Source")
        #plt.ylabel("Energy Output (GJ)")
        plt.title(average)

    plt.xticks(x_pos, avgdata['TFK'])
    
    ann = plt.annotate('', xy=(0.5, 1), xycoords='data',
                  xytext=(1, 1), textcoords='data',
                  arrowprops=dict(arrowstyle="-",
                                  connectionstyle="bar",
                                  ec="k",
                                  shrinkA=5, shrinkB=5))

    plt.show()
    
plot_data('Toughness', 'Stdev')
plot_data('Elong_at_break', 'Elong_std')
plot_data('UTensileStr', 'UStr_std')
plot_data('E', 'E_std')



def t_test(tfks): # example input: [[0,5], [0,15], [0,25]]
    import scipy.stats
    for i in tfks:
        tfk1 = i[0]
        tfk2 = i[1]
        data_subset1 = data.loc[data['TFK'] == tfk1, ['Toughness', 'Elong_at_break', 'UTensileStr', 'E']]
        data_subset2 = data.loc[data['TFK'] == tfk2, ['Toughness', 'Elong_at_break', 'UTensileStr', 'E']]
        #print(data_subset1)
        p_value = scipy.stats.ttest_ind(data_subset1, data_subset2)
        print(p_value)
        print ('\nt_test between', tfk1, ',' ,tfk2)
        print('Toughness:', p_value[1][0])
        print('Elongation:', p_value[1][1])
        print('Ultimate_Tensile:', p_value[1][2])
        print('E:', p_value[1][3])
        
t_test([[0,5],[0,15],[0,25],[5,15],[15,25],[5,25]])
    
