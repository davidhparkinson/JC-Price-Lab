
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#**IMPORTANT**: SYS is the package for getting input from the termal
import sys
import csv
import plotly.express as px
import glob
import math
import scipy.stats

def get_file_name_from_path(file):
    import os
    return file.split(os.path.sep)[-1]

def Scatter(file1,file2,comparison,toCheck,PoI,showGraphs,tissue,style,detail):
    #fig, ax = plt.subplots()
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}
    NoSepCov = 'NoSep Coverage (%)'
    SepCov = 'Sep Coverage (%)'
    NoUni = 'NoSepUnique'
    SepUni = 'SepUnique'
    acc = "Accession"
    accP = "Protein Accession"
    files = [file1,file2]
    lengths = [len(file1),len(file2)]
    names = [None]*len(files)
    datacomb = [None]*len(files)
    datatot = [None]*len(files)
    datacombSD = [None]*len(files)
    countList = [None]*len(files)
    data1comb = pd.DataFrame(data = None)
    data2comb = pd.DataFrame(data = None)
    PTMofInt = PoI
    more1 = 0
    more2 = 0
    equal = 0
    more1all = 0
    more2all = 0
    equalall = 0
    if(style=='c'):
        style = 'inner'
    elif(style=='a'):
        style = 'outer'

    if(toCheck=='c' or toCheck=='q'):
        if(toCheck =='c'):
            if(comparison=='s'):
                names[0] = ['IPS Coverage (%)','ProteinUnique'," Protein ",'Coverage (%)','Coverage','sdPro','#Unique']
                names[1] = ['PS Coverage (%)','PeptideUnique'," Peptide ",'Coverage (%)','Coverage','sdPep','#Unique']
            elif(comparison=='c'):
                names[0] = ['NoSep Coverage (%)','NoSepUnique'," No ",'Coverage (%)','Coverage','sdNo','#Unique']
                names[1] = ['Sep Coverage (%)','SepUnique'," ",'Coverage (%)','Coverage','sdSep','#Unique']
        elif(toCheck=='q'):
            if(comparison=='s'):
                names[0] = ['IPS Peptide Count','ProteinUnique'," Protein ",'#Peptides','Peptides','sdPro','#Unique']
                names[1] = ['PS Peptide Count','PeptideUnique'," Peptide ",'#Peptides','Peptides','sdPep','#Unique']
            elif(comparison=='c'):
                names = ['NoSep Peptide Count','Sep Peptide Count','NoSepUnique','SepUnique'," No "," ",'#Peptides','Peptides']
                names[0] = ['NoSep Peptide Count','NoSepUnique'," No ",'#Peptides','Peptides','sdNo','#Unique']
                names[1] = ['Sep Peptide Count','SepUnique'," ",'#Peptides','Peptides','sdSep','#Unique']
        
        for k in range(len(files)):
            i=0
            data = [None]*len(files[k])
            datacomb[k] = pd.DataFrame(data = None)
            for filename in files[k]:
                data[i] = pd.read_csv(filename,index_col=acc,header=0,usecols=[acc,names[k][3],names[k][6]])
                data[i] = data[i].rename({names[k][3]: names[k][0], names[k][6]: names[k][1]}, axis='columns')
                data[i]["count"]=1
                datacomb[k] = datacomb[k].append(data[i])
                i+=1
            countList[k] = datacomb[k].groupby([acc]).sum()
            for ind in countList[k].index:
                df1 = pd.DataFrame([[0,0,0]],columns=[names[k][0],names[k][1],"count"], index=[ind])
                df1.index.name = acc
                if(countList[k].loc[ind,"count"]==(lengths[k]-2)):
                    datacomb[k]=datacomb[k].append(df1)
                    datacomb[k]=datacomb[k].append(df1)
                if(countList[k].loc[ind,"count"]==(lengths[k]-1)):
                    datacomb[k]=datacomb[k].append(df1)
            
            datacombSD[k] = datacomb[k].groupby(acc)[names[k][0]].std()
            datacomb[k] = datacomb[k].groupby(acc).sum()/lengths[k]
            datacomb[k][names[k][5]]=datacombSD[k]

        #choose inner our outer merge. Inner is only common proteins, outer is union of all.
        dataJoin = datacomb[0].merge(datacomb[1], how=style, on=acc)
        dataJoin = dataJoin.fillna(0)
        for j in dataJoin.index:
            if(pd.isna(dataJoin[names[0][0]][j])):
                dataJoin[names[0][0]][j] = 0
                dataJoin[names[0][1]][j] = 0
                dataJoin[names[0][5]][j] = 0
            if(pd.isna(dataJoin[names[1][0]][j])):
                dataJoin[names[1][0]][j] = 0
                dataJoin[names[1][1]][j] = 0
                dataJoin[names[1][5]][j] = 0
        dataJoin['Abundance'] = dataJoin[names[0][1]] + dataJoin[names[1][1]]
        dataJoin['tStat']=((dataJoin[names[0][0]]-dataJoin[names[1][0]])/(np.sqrt((dataJoin[names[0][5]].pow(2)/float(len(file1)))+(dataJoin[names[1][5]].pow(2))/float(len(file2)))))
                                                                                        
        alpha = 0.05
        t=scipy.stats.t.ppf(1 - alpha / 2, len(file1)+len(file2)-2)                                                                                    
        dataJoin['Keep?']=False

        for i in dataJoin.index:
            if(dataJoin['tStat'][i] >= t):
                more1+=1
                #dataJoin['Keep?'][i] = True
                dataJoin.loc[i,"Keep?"] = True
            elif(dataJoin['tStat'][i] <= -t):
                more2+=1
                #dataJoin['Keep?'][i] = True
                dataJoin.loc[i,"Keep?"] = True
            else:
                equal+=1
                #dataJoin['Keep?'][i] = False
                dataJoin.loc[i,"Keep?"] = False
        for i in dataJoin.index:
            if(dataJoin[names[0][0]][i] > dataJoin[names[1][0]][i]):
                more1all+=1
            elif(dataJoin[names[0][0]][i] < dataJoin[names[1][0]][i]):
                more2all+=1
            else:
                equalall+=1
        dataJoin2 = dataJoin
        dataJoin = dataJoin[dataJoin['Keep?']==True]
        dataJoin2 = dataJoin2[dataJoin2['Keep?']==False]
        dataJoin['rank'] = dataJoin['Abundance'].rank()*100/len(dataJoin['Abundance'])

        print('\033[1m' + tissue+" "+detail+":"+'\033[0m')
        print("More "+names[0][4]+" with"+names[0][2]+"separation: " + str(more1all) + " proteins, ("+str(more1)+" significant).")
        print("More "+names[1][4]+" with"+names[1][2]+"separation: " + str(more2all) + " proteins, ("+str(more2)+" significant).")
        print("Equal "+names[0][4]+": " + str(equalall) + " proteins, ("+str(equal)+" are insignificant).")
        print("")
        print('p-value: '+str(alpha))
        print("")
        plt.scatter(dataJoin2[names[0][0]], dataJoin2[names[1][0]], c='gray',alpha=0.2)
        plt.scatter(dataJoin[names[0][0]], dataJoin[names[1][0]],c=dataJoin['rank'], cmap='viridis')
        plt.title(names[0][4] + " in " + tissue+" "+detail)
        plt.xlabel(names[0][0])
        plt.ylabel(names[1][0])
        cbar = plt.colorbar()
        cbar.set_label("Abundance Percentile")
        #ax1 = dataJoin.plot.scatter(x=names[0][0], y=names[1][0], c='rank',colormap='viridis')
        x = np.linspace(0, 100, 1000)
        plt.plot(x, x + 0, linestyle='solid',color='r')
        if showGraphs:
            plt.show()
    elif(toCheck=='p'):
        #create lists of strings needed for each set of parameters
        if(comparison=='s'):
            names[0] = ['IPS PTMs'," Protein ",'sdPro','PTM count']
            names[1] = ['PS PTMs'," Peptide ",'sdPep','PTM count']
        elif(comparison=='c'):
            names[0] = ['NoSep PTMs'," No ",'sdNo','PTM count']
            names[1] = ['WithSep PTMs'," ",'sdSep','PTM count']
        
        #loop to wrangle data for each of the 2 groups
        for k in range(len(files)):
            i=0
            data = [None]*len(files[k])
            datacomb[k] = pd.DataFrame(data = None)
            datatot[k] = pd.DataFrame(data = None)
            #goes through each of the replicate files in the group and append them all together, after filtering out PTMs
            for filename in files[k]:
                data[i] = pd.read_csv(filename,usecols=[accP,'PTM'])
                data[i] = pd.concat((data[i],data[i].PTM.str.get_dummies(sep='; ')),axis=1)
                data[i] = data[i].groupby(accP).sum()
                #these 2 PTMs are undesired
                data[i] = data[i].drop(columns = ["Pyro-glu from Q"])
                data[i] = data[i][data[i].columns.drop(list(data[i].filter(regex='Carbamidomethylation')))]
                #filters for the specific PTM, if specified
                if(PTMofInt != 'all'):
                    data[i] = data[i].filter(regex=PTMofInt)
                datacomb[k] = datacomb[k].append(data[i])
                i+=1
            
            #lets us know how many times each protein was seen though all of the replicates
            datatot[k] = datacomb[k].sum(axis=1)
            datatot[k] = datatot[k].to_frame()
            datatot[k] = datatot[k].rename(columns={0:names[k][0]})
            datatot[k] = datatot[k].reset_index()
            datatot[k]["count"]=[1]*len(datatot[k])
            #print(datatot[k])

            #adds additional rows of 0's for each instance that the protein is not seen in the replicates
            countList[k] = datatot[k].groupby([accP]).sum()
            for ind in countList[k].index:
                df1 = pd.DataFrame([[ind,0,0]],columns=[accP,names[k][0],"count"])
                if(countList[k].loc[ind,"count"]==(lengths[k]-2)):
                    datatot[k]=datatot[k].append(df1)
                    datatot[k]=datatot[k].append(df1)
                if(countList[k].loc[ind,"count"]==(lengths[k]-1)):
                    datatot[k]=datatot[k].append(df1)
            #gets the SDs and averages for each Protein
            datacombSD[k] = datatot[k].groupby(accP)[names[k][0]].std()
            datatot[k] = datatot[k].groupby(accP).sum()/lengths[k]
            datatot[k][names[k][2]]=datacombSD[k]
            #datatot[k] = datatot[k][datatot[k][names[k][0]]>0]
            #print(datatot[k])

        #type of merge, inner (common) or outer (all)
        dataJoin = datatot[0].merge(datatot[1], how=style, on=accP)
        #REPLACE ALL NaNs with 0's
        dataJoin = dataJoin.fillna(0)
        #filters out proteins with 0 of the PTM seen in either sample (for when looking at specific PTMs)
        dataJoin["sum"]=dataJoin[names[0][0]]+dataJoin[names[1][0]]
        dataJoin = dataJoin[dataJoin['sum']>0]
        #T-test for all proteins
        dataJoin['tStat']=((dataJoin[names[0][0]]-dataJoin[names[1][0]])/(np.sqrt((dataJoin[names[0][2]].pow(2)/float(len(file1)))+(dataJoin[names[1][2]].pow(2))/float(len(file2)))))
        #counting how many proteins are greater with each type
        for i in dataJoin.index:
            if(dataJoin[names[0][0]][i] > dataJoin[names[1][0]][i]):
                more1all+=1
            elif(dataJoin[names[0][0]][i] < dataJoin[names[1][0]][i]):
                more2all+=1
            elif(dataJoin[names[0][0]][i] == dataJoin[names[1][0]][i]):
                equalall+=1                                                                          
        
        #p-value and t-stat calcs
        alpha = 0.05
        t=scipy.stats.t.ppf(1 - alpha / 2, len(file1)+len(file2)-2)                                                                                    
        dataJoin['Keep?']=''
        #counting # of significant proteins
        for i in dataJoin.index:
            if(dataJoin['tStat'][i] >= t):
                more1+=1
                dataJoin.loc[i,"Keep?"] = True
            elif(dataJoin['tStat'][i] <= -t):
                more2+=1
                dataJoin.loc[i,"Keep?"] = True
            elif(dataJoin['tStat'][i] > -t and dataJoin['tStat'][i] < t):
                equal+=1
                dataJoin.loc[i,"Keep?"] = False
        #print(dataJoin)
        dataJoin['Abundance'] = dataJoin[names[0][0]] + dataJoin[names[1][0]]
        dataJoin2 = dataJoin
        dataJoin = dataJoin[dataJoin['Keep?']==True]
        dataJoin2 = dataJoin2[dataJoin2['Keep?']==False]
        dataJoin['rank'] = dataJoin['Abundance'].rank()*100/len(dataJoin['Abundance'])
                
        print('\033[1m' + tissue+" "+detail+":"+'\033[0m')
        print("More PTMs with"+names[0][1]+"separation: " + str(more1all) + " proteins ("+str(more1)+" are significant).")
        print("More PTMs with"+names[1][1]+"separation: " + str(more2all) + " proteins ("+str(more2)+" are significant).")
        print("Equal PTM detection: " + str(equalall) + " proteins ("+str(equal)+" are insignificant).")
        print("")
        plt.scatter(dataJoin2[names[0][0]], dataJoin2[names[1][0]], c='gray',alpha=0.2)
        plt.scatter(dataJoin[names[0][0]], dataJoin[names[1][0]],c=dataJoin['rank'], cmap='viridis')
        if PoI =="all":
            PoI = "PTM"
        plt.title(PoI + " count in " + tissue+" "+detail)
        plt.xlabel(names[0][0])
        plt.ylabel(names[1][0])
        cbar = plt.colorbar()
        cbar.set_label("Abundance Percentile")
        x = np.linspace(0, 1000, 1000)
        plt.plot(x, x + 0, linestyle='solid',color='r')
        if showGraphs:
            plt.show()
        dataJoin['rank'] = dataJoin['sum'].rank()
        print(dataJoin)
    dataJoin["Rank Percentile"] = dataJoin["rank"]
    dataJoin = dataJoin.filter(items=[names[0][0],names[1][0],"Rank Percentile"])
    dataJoin = dataJoin.melt(id_vars=["Rank Percentile"], var_name='Separation Type', value_name=names[0][3])
    if showGraphs:
        fig = px.histogram(dataJoin,x="Rank Percentile",y=names[0][3],color="Separation Type",histfunc="avg",barmode='group',nbins=10)
        fig.update_traces(xbins_start=0,xbins_end=100)
        fig.show()
    return None

if __name__ == "__main__":
    tissue = input("What tissue? ")
    style = input("Look at common proteins or all proteins? (enter c or a): ")
    comparison = input("Comparing Separations? or against control? (enter s or c): ")
    if(comparison == 'c'):
        sepType = input("Which type of separation? ('Pro' or 'Pep): ")
    toCheck = input("PTMs or Coverage or Pep Quant (enter p or c or q): ")
    if(toCheck=='p'):
        PoI = input("PTM of interest (type one or put 'all'): ")
    else:
        PoI = 'all'
    print("")
    if(tissue!="all"):
        words1 = [None]*3
        words2 = [None]*3
        showGraphs = True
        if(tissue=="Liver"):
            groups = int(input("3 or 4 groups? "))
            detail = "("+str(groups)+" pools)"
            if(comparison=='s'):
                words1[2]="/"+str(groups)
                words2[2]="/"+str(groups)
            else:
                words1[2]=""
                words2[2]="/"+str(groups)
        elif(tissue=="Serum"):
            groups = input("Wide or thin? (enter w or t): ")
            if(groups=='w'):
                groups = "Wide"
            elif(groups=='t'):
                groups = "Thin"
            detail = "("+groups+")"
            if(comparison=='s'):
                words1[2]="/"+str(groups)
                words2[2]=""
            else:
                words1[2]=""
                words2[2]=""
        else:
            words1[2]=""
            words2[2]=""
            detail = ""
        if(toCheck=='p'):
            words1[0] = 'PTMs'
            words2[0] = 'PTMs'
        elif(toCheck=='c' or toCheck=='q'):
            words1[0] = 'Coverage'
            words2[0] = 'Coverage'
        if(comparison=="c"):
            words1[1]="Control"
            words2[1]=sepType
        else:
            words1[1]="Pro"
            words2[1]="Pep"
        path1 = "Paper MS Data/"+words1[0]+"/"+tissue+"/"+words1[1]+words1[2]
        path2 = "Paper MS Data/"+words2[0]+"/"+tissue+"/"+words2[1]+words2[2]
        file1 = glob.glob(path1 + "/*.csv")
        file2 = glob.glob(path2 + "/*.csv")
        Scatter(file1,file2,comparison,toCheck,PoI,showGraphs,tissue,style,detail)
    elif(tissue=="all"):
        #To check something for all 5 tissue settings
        words1 = [None]*3
        words2 = [None]*3
        TissueTypes = ["Brain","Muscle"]
        LivGroups = [3,4]
        SerumGroups = ['Wide','Thin']
        showGraphs = False
        for i in LivGroups:
            tissue = "Liver"
            groups = i
            detail = "("+str(groups)+" pools)"
            if(comparison=='s'):
                words1[2]="/"+str(groups)
                words2[2]="/"+str(groups)
            else:
                words1[2]=""
                words2[2]="/"+str(groups)
            if(toCheck=='p'):
                words1[0] = 'PTMs'
                words2[0] = 'PTMs'
            elif(toCheck=='c' or toCheck=='q'):
                words1[0] = 'Coverage'
                words2[0] = 'Coverage'
            if(comparison=="c"):
                words1[1]="Control"
                words2[1]=sepType
            else:
                words1[1]="Pro"
                words2[1]="Pep"
            path1 = "Paper MS Data/"+words1[0]+"/"+tissue+"/"+words1[1]+words1[2]
            path2 = "Paper MS Data/"+words2[0]+"/"+tissue+"/"+words2[1]+words2[2]
            file1 = glob.glob(path1 + "/*.csv")
            file2 = glob.glob(path2 + "/*.csv")
            Scatter(file1,file2,comparison,toCheck,PoI,showGraphs,tissue,style,detail)
        for j in SerumGroups:
            tissue = "Serum"
            groups = j
            detail = "("+groups+")"
            if(comparison=='s'):
                words1[2]="/"+str(groups)
                words2[2]=""
            else:
                words1[2]=""
                words2[2]=""
            if(toCheck=='p'):
                words1[0] = 'PTMs'
                words2[0] = 'PTMs'
            elif(toCheck=='c' or toCheck=='q'):
                words1[0] = 'Coverage'
                words2[0] = 'Coverage'
            if(comparison=="c"):
                words1[1]="Control"
                words2[1]=sepType
            else:
                words1[1]="Pro"
                words2[1]="Pep"
            path1 = "Paper MS Data/"+words1[0]+"/"+tissue+"/"+words1[1]+words1[2]
            path2 = "Paper MS Data/"+words2[0]+"/"+tissue+"/"+words2[1]+words2[2]
            file1 = glob.glob(path1 + "/*.csv")
            file2 = glob.glob(path2 + "/*.csv")
            Scatter(file1,file2,comparison,toCheck,PoI,showGraphs,tissue,style,sig,detail)
        for i in TissueTypes:
            tissue = i
            detail = ""
            words1[2]=""
            words2[2]=""
            if(toCheck=='p'):
                words1[0] = 'PTMs'
                words2[0] = 'PTMs'
            elif(toCheck=='c' or toCheck=='q'):
                words1[0] = 'Coverage'
                words2[0] = 'Coverage'
            if(comparison=="c"):
                words1[1]="Control"
                words2[1]=sepType
            else:
                words1[1]="Pro"
                words2[1]="Pep"
            path1 = "Paper MS Data/"+words1[0]+"/"+tissue+"/"+words1[1]+words1[2]
            path2 = "Paper MS Data/"+words2[0]+"/"+tissue+"/"+words2[1]+words2[2]
            file1 = glob.glob(path1 + "/*.csv")
            file2 = glob.glob(path2 + "/*.csv")
            Scatter(file1,file2,comparison,toCheck,PoI,showGraphs,tissue,style,detail)