import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#**IMPORTANT**: SYS is the package for getting input from the termal
import sys
import csv
import glob
import math

#This grabs the filename from the csv
def get_file_name_from_path(file):
    import os
    return file.split(os.path.sep)[-1]

#this function does all of the work
def sort_vialsPro(groupNum,peaks_filename,vials_filename,signal_filenames,cutoffs,width_factor,ProteinTotal,mass_units,isTester,isColor):
    #Here we read in the peaks and fractions csv's. This code may vary for different HPLC system outputs
    #for the peaks file, it is essential that we have Time, Area, Height, Width, Start, and End
    #for the vials file, it is essential that we have vial #, start, end, and volume.
    peaks = pd.read_csv(peaks_filename,index_col=0,names=["Time","Type","Area","Height","Width","Start","End"])
    vials = pd.read_csv(vials_filename,index_col='vial',
                    names=['vial','time','label0','label1','idk1','start','end','volume','extra1','extra2','extra3','extra4','extra5','extra6','extra7','extra8','extra9','extra10','extra11'])
    vials = vials.drop(columns=['time','idk1','label0','extra1','extra2','extra3','extra4','extra5','extra6','extra7','extra8','extra9','extra10','extra11'])
    vialSepLineStyle = "dashed"
    if (not isTester):
        #create plot
        fig, ax = plt.subplots()
        font = {'weight' : 'bold','size'   : 20}
        matplotlib.rc('font', **font)
        #This grabs and plots each of the signal files if there are replicates
        for filename in signal_filenames:
            data = pd.read_csv(filename)
            x = data.iloc[:,0].to_numpy()
            y = data.iloc[:,1].to_numpy()
            slope = (y[-1]-y[0])/(x[-1]-x[0])
            y_new = y - x*slope
            title=get_file_name_from_path(filename)[:-4]
            ax.plot(x, y_new,label=title)
        #plot fine tuning
        ax.plot(label = "Concatenatino Scheme")
        peaksYGraph = peaks[peaks["Time"]>min(vials.start)]
        peaksYGraph = peaksYGraph[peaksYGraph["Time"]<max(vials.end)]
        plt.ylim(-0.2,max(peaks.Height)+1)
        plt.xlim=(peaksYGraph.iloc[0],peaksYGraph.iloc[-1])
        plt.xlabel('Time (min)',fontsize=16)
        plt.ylabel('Absorbance (mAU)',fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.title(str(groupNum)+" groups", fontsize=20)
        plt.legend()

    #create a Area Percent column in the vials dataframe
    vialAreas=[]
    for vialIndex in vials.index:
        AreaSum = 0.0
        for peakIndex in peaks.index:
            if (peaks.Time[peakIndex]>=(vials.start[vialIndex])) and (peaks.Time[peakIndex]<=(vials.end[vialIndex])):
                AreaSum += peaks.Area[peakIndex]
        vialAreas.append(AreaSum)
    vials['area'] = vialAreas
    vials['AreaPer'] = vials.area/(vials.area.sum())

    #create lists of to put all of the peaks and vials for each group
    peakGroups=[None]*groupNum
    vialGroups=[None]*groupNum
    isVials=[None]*groupNum
    colOptions = ['red','orange','yellow','green','blue','purple']
    shadings = ['\\','/', '-', '+', 'x', 'o', 'O', '.', '*']
    for num in range(groupNum):
        peakGroups[num] = []
        vialGroups[num] = []
        isVials[num] = [False]
    
    #create a 2D list of peak start and stop times for each peak, spread based on the width factor, classified by group
    for j in peaks.index:
        for m in range(len(cutoffs)):
            if (peaks.Area[j]>=cutoffs[m]):
                peakStart = 0.0
                peakEnd = 0.0
                peakStart = peaks.Time[j]-peaks.Width[j]*width_factor
                peakEnd = peaks.Time[j]+peaks.Width[j]*width_factor
                peakGroups[m].append([peakStart,peakEnd])
    if (isTester):
        for k in vials.index:
            match = False
            for num in range(groupNum):
                for l in range(len(peakGroups[num])):
                    if (vials.start[k]>=peakGroups[num][l][0] and vials.start[k]<=peakGroups[num][l][1]) or (vials.end[k]>=peakGroups[num][l][0] and vials.end[k]<=peakGroups[num][l][1]) or (vials.start[k]<=peakGroups[num][l][0] and vials.end[k]>=peakGroups[num][l][1]) or (vials.start[k]>=peakGroups[num][l][0] and vials.end[k]<=peakGroups[num][l][1]):
                        match = True
                        vialGroups[num].append(vials.label1[k])
                        break
                if match==True:
                    break
            if match==False:
                vialGroups[groupNum-1].append(vials.label1[k])
    if (not isTester):
        #For each vial, check to see if it overlaps with the highest abundance peaks, then if not, continue to check, in order of abudnance classification
        for k in vials.index:
            match = False
            for num in range(groupNum):
                for l in range(len(peakGroups[num])):
                    if (vials.start[k]>=peakGroups[num][l][0] and vials.start[k]<=peakGroups[num][l][1]) or (vials.end[k]>=peakGroups[num][l][0] and vials.end[k]<=peakGroups[num][l][1]) or (vials.start[k]<=peakGroups[num][l][0] and vials.end[k]>=peakGroups[num][l][1]) or (vials.start[k]>=peakGroups[num][l][0] and vials.end[k]<=peakGroups[num][l][1]):
                        match = True
                        vialGroups[num].append(vials.label1[k])
                        if isColor == "c":
                            plt.axvspan(vials.start[k],vials.end[k],color=colOptions[num],alpha=0.2)
                        else:
                            plt.axvspan(vials.start[k],vials.end[k],hatch=shadings[num],alpha=0)
                        if vials.start[k]<=peakGroups[num][l][0] and (vials.end[k]>=peakGroups[num][l][0] and vials.end[k]<=peakGroups[num][l][1]):
                            plt.vlines(vials.start[k],0,200,colors=colOptions[num], linestyles=vialSepLineStyle)
                            plt.vlines(vials.end[k],0,200,colors=colOptions[num], linestyles=vialSepLineStyle)
                        elif (vials.start[k]>=peakGroups[num][l][0] and vials.start[k]<=peakGroups[num][l][1]) and vials.end[k]>=peakGroups[num][l][1]:
                            plt.vlines(vials.start[k],0,200,colors=colOptions[num], linestyles=vialSepLineStyle)
                            plt.vlines(vials.end[k],0,200,colors=colOptions[num], linestyles=vialSepLineStyle)
                        elif (vials.start[k]<=peakGroups[num][l][0] and vials.end[k]>=peakGroups[num][l][1]):
                            plt.vlines(vials.start[k],0,200,colors=colOptions[num], linestyles=vialSepLineStyle)
                            plt.vlines(vials.end[k],0,200,colors=colOptions[num], linestyles=vialSepLineStyle)
                        elif (vials.start[k]>=peakGroups[num][l][0] and vials.end[k]<=peakGroups[num][l][1]):
                            plt.vlines(vials.start[k],0,200,colors=colOptions[num], linestyles=vialSepLineStyle)
                            plt.vlines(vials.end[k],0,200,colors=colOptions[num], linestyles=vialSepLineStyle)
                        break
                if match==True:
                    break
            if match==False:
                vialGroups[groupNum-1].append(vials.label1[k])

    #create lists of group names and vials in each group
    names = []
    groups = []
    for z in range(len(vialGroups)):
        names.append("Group "+str(z+1))
        groups.append([int(v[4:]) for v in vialGroups[z]]) 

    if(isTester):
        totals = [0]*groupNum
        for i in range(len(names)):
            if vials.loc[groups[i]].empty:
                total = 0
                totalVol=0
            else:
                totals[i] = vials.loc[groups[i]].AreaPer.sum()*ProteinTotal/100
        return totals      
    if (not isTester):
        dfs = {}
        results = []
        totals = []
        #This prints all of the info we want
        for i in range(len(names)):
            print(names[i]+": "+str(groups[i]))
            if vials.loc[groups[i]].empty:
                total = 0
                totalVol=0
            else:
                result = names[i] + ": " + str(groups[i])
                results.append(result)
                total = vials.loc[groups[i]].AreaPer.sum()*ProteinTotal
                totalVol = vials.loc[groups[i]].volume.sum()
            print(str("{:.2f}".format(total)), mass_units,"of",names[i],"protein in ",str("{:.2f}".format(totalVol/1000))," mL")
            result = str("{:.2f}".format(total))+ " " + mass_units+ " " +"of"+" " +names[i]+" protein in "+str("{:.2f}".format(totalVol/1000))+" mL"
            results.append(result)
            totals.append(total)
        #Gives us a plot showing the concatenation scheme.
        plt.show()
        return results

def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

def sort_vialsPep(protein_total,mass_units,frac_num,vials_filename,signal_filenames,baseline_filename,isColor):
    fig, ax = plt.subplots()
    for filename in signal_filenames:
        data = pd.read_csv(filename)
        x = data.iloc[:,0].to_numpy()
        y = data.iloc[:,1].to_numpy()

        data2 = pd.read_csv(baseline_filename,header = None)
        x2 = data2.iloc[:,0].to_numpy()
        y2 = data2.iloc[:,1].to_numpy()

        if len(y) > len(y2):
            y = y[0:len(y2)]
            x = x2
        if len(y2) > len(y):
            y2 = y2[0:len(y)]

        y_new=y-y2
    
        #title=get_file_name_from_path(signal_filenames)[:-11]
        title = "Peptide Separation"
        ax.plot(x, y_new,label=title) #linestyle = '--'
    #plt.vlines(10,0,100)
    #plt.plot(x1, y1, color = 'b',linewidth=1, label = signal_filename1;x2,y2,color = 'r',linewidth=1, label = signal_filename2)


    vials = pd.read_csv(vials_filename,index_col='vial',
                    names=['vial','time','label0','label1','idk1','start','end','volume','extra1','extra2','extra3','extra4','extra5','extra6','extra7','extra8','extra9','extra10','extra11'])
    vials = vials.drop(columns=['time','idk1','label0','extra1','extra2','extra3','extra4','extra5','extra6','extra7','extra8','extra9','extra10','extra11'])
    plt.xlabel('Time (min)', fontsize=16)
    plt.ylabel('Absorbance (mAU)',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.legend()
    colOptions = ['red','orange','yellow','green','blue','purple']
    shadings = ['\\','/', '-', '+', 'x', 'o', 'O', '.', '*']

    total_frac = len(vials)
    fracs_per = int(math.ceil(total_frac/frac_num))
    GroupAreas=[]
    GroupVols=[]
    GroupMasses=[]
    big_list=[]
    for i in range(frac_num):
        little_list = []
        for j in range(fracs_per):
            if (((i+1)+j*frac_num)>total_frac):
                break
            else:
                little_list.append(((i+1)+j*frac_num))
        big_list.append(little_list)
    for group in range(len(big_list)):
        AreaSum = 0.0
        VolSum=0.0
        for vial in range(len(big_list[group])):
            index=big_list[group][vial]
            plt.vlines(vials.start.iloc[index-1],ymax=150,ymin=0,colors='black', linestyles='dashed')
            if isColor=="c":
                plt.axvspan(vials.start.iloc[index-1],vials.end.iloc[index-1],color=colOptions[group],alpha=0.2)
            else:
                plt.axvspan(vials.start.iloc[index-1],vials.end.iloc[index-1],hatch=shadings[group],alpha=0)
            for k in range(find_nearest(x,vials.start.iloc[index-1]),find_nearest(x,vials.end.iloc[index-1])):
                x1v = k
                x2v = k+1
                AreaSum += (x[x2v]-x[x1v])*(y_new[x1v]+y_new[x2v])/2
            VolSum+=vials.volume.iloc[index-1]
        GroupAreas.append(AreaSum)
        GroupVols.append(VolSum)
    GroupMasses=(GroupAreas/sum(GroupAreas))*protein_total
    for i in range(len(GroupMasses)):
        print("Group "+str(i+1)+" vials: "+str(big_list[i]))
        print(str("{:.2f}".format(GroupMasses[i])), mass_units,"of peptide in ",str("{:.2f}".format(GroupVols[i]/1000))," mL")
    plt.show()

if __name__ == "__main__":
    #asks for all variables
    sepType = input("Protein or Peptide separation? (enter 'pro' or 'pep'): ")
    groupNum = int(input("How many groups are there? "))
    mass_total = float(input("How many µg of protein is there? "))
    mass_units = "µg"
    colors = input("Colors or shading? (c/s): ")

    path1 = "HPLC/Peaks"
    path2 = "HPLC/Fractions"
    path3 = "HPLC/Signal"
    path4 = "HPLC/Baseline"
    peaks_filename = glob.glob(path1 + "/*.csv")
    vials_filename = glob.glob(path2 + "/*.csv")
    signal_filenames = glob.glob(path3 + "/*.csv")
    baseline_filename = glob.glob(path4 + "/*.csv")

    if(sepType == "pro"):
        width_factor = float(input("Width factor: "))
        cutoffs = [0]*(groupNum-1)
        path1 = "HPLC/Pro/Peaks"
        path2 = "HPLC/Pro/Fractions"
        path3 = "HPLC/Pro/Signal"
        peaks_filename = glob.glob(path1 + "/*.csv")
        vials_filename = glob.glob(path2 + "/*.csv")
        signal_filenames = glob.glob(path3 + "/*.csv")
        #asks whether cutoffs will be input manually or automatically
        oops = False
        while oops == False:
            isAuto = str(input("Choose cutoffs automatically or manually? (type: 'a' or 'm') "))
            if (isAuto == 'a'):
                print("Oops, we don't know how to do that yet!")
                print("")
            if (isAuto == 'm'):
                oops = True
                isTester = False
                cutoffs=[None]*(groupNum-1)
                for i in range(len(cutoffs)):
                    cutoffs[i] = float(input("Cutoff "+str(i+1)+": "))
        sort_vialsPro(groupNum,peaks_filename[0],vials_filename[0],signal_filenames,cutoffs,width_factor,mass_total,mass_units,isTester,colors)
    if(sepType == "pep"):
        path2 = "HPLC/Pep/Fractions"
        path3 = "HPLC/Pep/Signal"
        path4 = "HPLC/Pep/Baseline"
        vials_filename = glob.glob(path2 + "/*.csv")
        signal_filenames = glob.glob(path3 + "/*.csv")
        baseline_filename = glob.glob(path4 + "/*.csv")
        sort_vialsPep(mass_total,mass_units,groupNum,vials_filename[0],signal_filenames,baseline_filename[0],colors)

