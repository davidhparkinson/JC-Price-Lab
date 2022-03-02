import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#**IMPORTANT**: SYS is the package for getting input from the termal
import sys
import csv

def get_file_name_from_path(file):
    import os
    return file.split(os.path.sep)[-1]

def tk_get_single_file(extension='*', prompt="Select file"):
   from tkinter import filedialog
   from tkinter import Tk
   import fnmatch
   import os
   if extension[0] != '*':
      extension = '*' + extension
   root = Tk()
   root.withdraw()
   if (extension == '*'):
      root.filename = filedialog.askopenfilenames(
         initialdir=os.getcwd(), title=prompt)
   else:
      root.filename = filedialog.askopenfilenames(
         initialdir=os.getcwd(), title=prompt,
         filetypes=((extension[2:] + " files", extension)))
   root.update()
   
   filename = root.filename
   return filename  # A string representing the file path

def plot_HPLC(signal_filenames):
    #print(len(signal_filenames))
    #fig, ax = plt.subplots()
    fig, ax = plt.subplots()
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}
    for filename in signal_filenames:
      data = pd.read_csv(filename)
      x = data.iloc[:,0].to_numpy()
      y = data.iloc[:,1].to_numpy()
      slope = (y[-1]-y[0])/(x[-1]-x[0])
      y_new = y - x*slope
      #x = x/3
      title=get_file_name_from_path(filename)[:-4]
      ax.plot(x, y_new,label=title) #linestyle = '--'
      #plt.vlines(10,0,100)
    #plt.plot(x1, y1, color = 'b',linewidth=1, label = signal_filename1;x2,y2,color = 'r',linewidth=1, label = signal_filename2)
    plt.xlabel('Time (min)', fontsize=16)
    plt.ylabel('Absorbance (mAU)',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.xlim([1,50])
    #plt.ylim([-0.7,20.5])
    plt.title('Combined Tissues')
    plt.legend(fontsize=16)
    plt.show()
    #plt.savefig(signal_filename1[:-4]+'.png')
    return None

plot_HPLC(tk_get_single_file(prompt="Select HPLC file"))
