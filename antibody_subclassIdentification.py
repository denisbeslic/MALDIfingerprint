# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 12:44:49 2020
@author: Denis Beslic
Processing of raw MS data
Determination of Antibody Subclass
"""


from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from BaselineRemoval import BaselineRemoval
import pandas as pd
import numpy as np
import glob, os
from scipy.signal import find_peaks
from scipy import stats
import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import mpld3
from mpld3 import plugins
import plotly.express as px
import plotly.graph_objects as go
import chart_studio.plotly as py
import chart_studio.tools as tls
from pathlib import Path

'''
preprocessing, peakdetection, peakfiltering
input: 
    - dataframe with mz and intensity
output: 
    - dataframe with mz, intensity and peak
parameters:
    - smoothwindow: windowsize of savitzky-golay filter
    - threshold1: Calculates the height for first threshold (up to 4500 m/z).
    Number is multiplied by global MAD value of the MS data. F.e. 6*mad is then
    used as a height threshold
    - threshold2: same as threshold 1, only for data over a mass of 4500

thresholds are seperated in low and high mass area (at 4500), so the user is 
able to choose a lower threshold for the high mass area, because the intensity
gets lower at a higher mass. Alternatively, the user can just choose to use the
same value for threshold1 and threshold2.
'''
def processing(msfile, smoothwindow = 41, threshold1 = 6, threshold2 = 4):
    
    msfile = msfile.rename(columns={0: "mz", 1: "intensity"})
    
    
    '''
    PREPROCESSING
    '''
    
    
    intensity = msfile["intensity"]
    # normalization with sqaure root
    intensity = np.sqrt(abs(intensity))
    # smoothing with savitsky-golay filter
    intensity = savgol_filter(intensity, smoothwindow, 3)
    # remove baseline with ZhangFit
    intensity = BaselineRemoval(intensity)
    intensity = intensity.ZhangFit()
    
    
    # remove old intensity column and add new preprocessed intensity column
    msfile = msfile.drop(columns = ["intensity"])
    msfile.insert(1, "intensity", intensity)
    
    
    '''
    PEAK DETECTION with scipy.signal_find.peaks
    '''
    
    # scipy.signal.find_peaks
    intensity = msfile["intensity"]
    mz = msfile["mz"]

    # pick two thresholds to detect peaks in high mass region

    mad = stats.median_absolute_deviation(msfile["intensity"])
    TH1= threshold1*mad
    TH2= threshold2*mad
    
    # border where the second threshold starts is defined at 4500 m/z
    detectionborder = 0
    for i in list(msfile["mz"]):
        if i > 4500:
            detectionborder = i
            break
    
    
    index1 = list(mz).index(detectionborder)
    peaks, _ = find_peaks(msfile["intensity"][:index1], distance = None, height = TH1)

    
    #second threshold
    peaks2, _ =  find_peaks(msfile["intensity"][index1:], distance = None, height= TH2)
    peaks2 = peaks2 + index1
    peaks = np.append(peaks, peaks2)
    
    # peaks gives only index at which a peak is detected
    # use peaks to create array of 0 (no peak) and 1 (peak)
    signals = np.zeros(len(msfile["intensity"]))
    for i in peaks:
        signals[i] = 1
        
    msfile.insert(2, "peak", signals)    

    

    '''
    PEAK FILTERING to pick only monoisotopic peaks
    '''
    
    # separate the peaklist to perform different postprocessing according to mass
    # reflectormode does not work above a certain mass
    
    peaksSmallMass = []
    peaksBigMass = []
    peakborder = 4500
    
    intensity = msfile["intensity"]
    mz = msfile["mz"]
    peak = msfile[msfile["peak"]==1]
    peak = peak["peak"]
    peak = peak.index.tolist()
    
    
    for j in peak:
        if mz[j] < peakborder:
            peaksSmallMass.append(j)
        else:
            peaksBigMass.append(j)
            
            
    '''
    Reflector Mode: Just pick smallest peak to get monoisotopic peak
    '''      
            
    # flip the list so you have sorted from highest to lowest peak
    revpeaklist = np.flip(peaksSmallMass)
    #FilteredPeaks will contain only index of  monoisotopic peaks
    FilteredPeakIndex = []
    # go through list and remove peaks which are distanced less than 1.5 dalton
    # next to each other. Smallest is the monoistopic peak which stays in the list
    for i in range(len(revpeaklist)-1):
        if mz[revpeaklist[i]] - mz[revpeaklist[i+1]] > 1.3:
            FilteredPeakIndex.append(revpeaklist[i])
            

    '''
    Linear Mode: Calculate Mean based on intensity for post-processing
    '''
            
            
    revBigpeaklist = np.flip(peaksBigMass)
    finallistHighMass = []
    if revBigpeaklist != []:
       finallistHighMass.append(revBigpeaklist[-1])
       
    # arithMean is a list to put all peaks together which are close to each other (<1.3)
    # to calculate the mean value 
    arithMean = []
    arithMeanSum = 0
    
    for k in range(len(revBigpeaklist)-1):
        # if next peak further away than 1.3 then just put peak in finallist
        if(abs(mz[revBigpeaklist[k]] - mz[revBigpeaklist[k+1]]) > 1.3):
            finallistHighMass.append(revBigpeaklist[k])
        # if next peak close then put peak in temporary list arithMean
        elif(abs(mz[revBigpeaklist[k]] - mz[revBigpeaklist[k+1]]) < 1.3):
            arithMean.append(revBigpeaklist[k])
        elif(abs(mz[revBigpeaklist[k]] - mz[revBigpeaklist[k+1]]) > 1.3) & (abs(mz[revBigpeaklist[k]] - mz[revBigpeaklist[k-1]]) < 1.3):
            for p in range(len(arithMean)):
                # calculate the mean based on the intesity of the different peaks in arithMean
                arithMeanSum = arithMeanSum + mz[arithMean[p]]*intensity[arithMean[p]]
            number = arithMeanSum/sum(intensity[arithMean])
            number = round(number,2)
            
            # add missing number which is not included in loop
            correctNumber=[]
            pop=0
            for cc in mz:
                pop = pop+1
                if cc==number:
                    correctNumber.append(pop-1)
                    break
                elif cc>number:
                    correctNumber.append(pop-1)
                    break


            finallistHighMass = finallistHighMass + correctNumber
            arithMeanSum=0
            arithMean=[]



    
    
    #finallist is index of all filtered peaks
    FilteredPeakIndex = finallistHighMass + FilteredPeakIndex
    FilteredPeakIndex.sort()
   
    print(str(len(FilteredPeakIndex)) + " peaks detected.")
    
    
    # create new list with only 0s
    newpeaks = [0]*len(mz)
    
    # fill this list with the correct peaks at the correct index
    for i in range(len(newpeaks)):
        if (FilteredPeakIndex.count(i) > 0):
            newpeaks[i] = 1
    
    # replace new peaklist with the unfiltered peaks
    msfile = msfile.drop(columns = ["peak"])
    msfile.insert(2, "peak", newpeaks)
    return msfile





'''
Determine Subclass. Comparing between file and library peaks. 
Generates summary table.
input: 
    - dataframe with mz, intensity, peak
    - path to directory with Masslists in .txt format f.e. IgG1.txt
      files need to have a header line with "m/z(mi)" column for masses
      columns need to be tab-seperated
      (Can also be used to create overlap of peaks with other MSdata)
output: 
    - table with overlap of peaks from your single msfile with masslists 
parameters:
    - tol: tolerance/ deviation which you accept for peak overlap
      tol = 0.3 means: peak1= 1.3, peak2= 1.5 are regarded as overlapping 
'''
def findSubclass(msfile, pathmasslist, tol = 0.3):
    # Pklist will contain mass of each masslist
    Pklist = []
    
    filenames = []
    os.chdir(pathmasslist)
    
    
    for file in glob.glob("*.txt"):
        
        # change header and sep if using other format
        listOfMasses2 = pd.read_table(file, sep="	", header = 0)
        start = 0
        end = file.find(".txt")
        file = file[start:end]
        filenames.append(file)
        Pklist.append(listOfMasses2)
    
    Masslists = pd.Series(Pklist, filenames)
    
    tablemz = msfile[msfile["peak"]==1]

    mz = tablemz["mz"] 
    
    # create table compare number of peaks between masslists and your MS file
    listeOverlap=[]
    for i in range(len(Masslists)):
        Masslists1 = np.array(Masslists[i]["m/z(mi)"])
        c = np.array(mz)[(np.abs(Masslists1[:,None] - np.array(mz)) < tol).any(0)]
        anzahl  = len(c)
        listeOverlap.append(anzahl)
        c = 0
        anzahl = 0
    overlapTable = pd.Series(listeOverlap, filenames)
    return overlapTable


'''
Determine Subclass for every peak. Comparing between file and library peaks. In contrast to
findSubclass, it does not produce a summary table, but gives back your data with
additional column with assigned subclass for each peak. 
input: 
    - dataframe with mz, intensity, peak
    - path to directory with Masslists in .txt format f.e. IgG1.txt
      files need to have a header line with "m/z(mi)" column for masses
      columns need to be tab-seperated
      Optionally, also include a "sequence" column for sequences of each peak
      (Can also be used to create overlap of peaks with other MSdata)
output: 
    - dataframe with mz, intensity, peak, subclass. Whereas subclass is a list
    of subclasses which share the same peak as your MS file.
parameters:
    - tol: tolerance / deviation which you accept for peak overlap
      tol = 0.3 means: peak1= 1.3, peak2= 1.5 are regarded as overlap 
'''

def PeakToClass(msfile, pathmasslist, outputname, tol = 0.3):
    # Pklist will contain mass of each masslist
    Pklist = []
    
    filenames = []
    os.chdir(pathmasslist)
    
    #depends on the dataformat: .txt files?
    for file in glob.glob("*.txt"):
        
        # change header and sep if using other format
        listOfMasses2 = pd.read_table(file, sep="	", header = 0)
        start = 0
        end = file.find(".txt")
        file = file[start:end]
        filenames.append(file)
        Pklist.append(listOfMasses2)
    
    Masslists = pd.Series(Pklist, filenames)
    

    tablemz = msfile[msfile["peak"]==1]
    mz = tablemz["mz"]

    
    # newsubclass will contain the name of the file which share the same peak
    newsubclass = [[]]*len(msfile["mz"])

 
        
        

    for i in range(len(Masslists)):
        # array format
        Masslists1 = np.array(Masslists[i]["m/z(mi)"])

        
        '''
         c gives mass peaks overlap from your input MS data 
         d vices mass peaks overlap from insilico library data
         f.e. if msfile = [1.1, 2.3, 2.9], libdata = [1.2, 2.1, 3.8]
         and tol = 0.3
         c gives out [1.1, 2.3]. d gives out [1.2, 2.1]
        '''
        c = np.array(mz)[(np.abs(Masslists1[:,None] - np.array(mz)) < tol).any(0)]

        
        # get index of your msfile and fill empty lists with names of subclass
        # where we got a peak
        for j in c:
            if(len(msfile.query("mz=="+str(j))) > 0):
                idxs = msfile.query("mz=="+str(j)).index[0]
                newsubclass[idxs] = newsubclass[idxs] + [filenames[i]]
                
    # append the subclass and sequence information to your file 
    msfile.insert(3, "subclass", newsubclass)
    return msfile
    


'''
Create Plot for MS data
input: 
    - dataframe with mz, intensity, peak
    - name of file for title and output
output: 
    - creates png file in result directory of MSfile  
'''
def plotMS(msfile, name):
    title = name.replace("_", " ")
    dots = msfile.query("peak == 1")

    plt.style.use(['science','grid'])
    dots = msfile.query("peak == 1")
    fig = plt.figure()
    
    fig.set_figheight(4)
    fig.set_figwidth(8)
    
    
    plt.plot(msfile["mz"], msfile["intensity"], zorder=1, linewidth=0.8, color="#314ca3")
    plt.xlabel('m/z')
    plt.ylabel('intensity')
    
    maxxtick = max(msfile["mz"])
    minxtick = min(msfile["mz"])
    x_ticks= np.arange(minxtick, maxxtick, 500)
    plt.xticks(x_ticks, rotation=90)
    
    plt.suptitle(title)
    plt.plot(dots["mz"], dots["intensity"], "x", markersize=4, zorder=2, color = "red", alpha=0.6)
    plt.savefig(name + '.png', dpi=600)

    
'''
Create interactive Plot for MS data
input: 
    - dataframe with mz, intensity, peak
    - name of file for title and output
output: 
    - creates html file in result directory of MSfile  
'''

def InteractivePlotMS(msfile, name):
    title = name.replace("_", " ")
    dots = msfile.query("peak == 1")
    overlap = msfile[msfile["subclass"].map(lambda d: len(d)) > 0]
    # how to get dataframe subset with peak ==1 & subclass nonempty
    
    fig = go.Figure()

    # create line for all datapoints
    fig.add_trace(go.Scatter(x=msfile["mz"], y=msfile["intensity"], name="",
                             hovertemplate='<i>mass</i>: %{x}'+'<br><i>intensity</i>: %{y}<br>' + '<extra></extra>'))
    
    # marks your detected peaks
    fig.add_trace(go.Scatter(x=dots["mz"], y=dots["intensity"],mode='markers', name='detected peak',
                             hovertemplate='<i>mass</i>: %{x}'+'<br><i>intensity</i>: %{y}<br>' + '<extra></extra>'))
    
    
    
    # marks overlap of your detected peaks with library
    fig.add_trace(go.Scatter(x=overlap["mz"], y=overlap["intensity"],mode='markers', 
                             text="<i>overlap with: </i> " + overlap["subclass"].astype(str) + "<br>"+
                             "<i>corresponding sequence: </i> " + overlap["sequence"].astype(str), name='database peak',
                             hovertemplate='<i>mass</i>: %{x}'+'<br><i>intensity</i>: %{y}<br>'+'%{text}'+'<extra></extra>',
                             textfont=dict(
                                 family="Computer Modern",
                                 size=18,
                                 color="LightSeaGreen"
                                 )
                             )
                  )    
    

    fig.update_layout(hovermode="closest")
    
    fig.update_layout(
    showlegend = True,
    title=name,
    xaxis_title="m/z",
    yaxis_title="intensity",
    font=dict(
        family="Computer Modern",
        size=18,
        color="black"
        ),
    xaxis = dict(
        tickmode = 'linear',
        tick0 = 0,
        dtick = 500
        ),
    xaxis_tickformat = '000',
    )


    #fig.savefig(name + '.png', dpi=600)
    fig.write_html(name +'.html', auto_open=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l","--library", help="source path to directory of theoretical spectra files for each subclass", type=str, required=True)
    parser.add_argument("-i","--input", help="source path to directory of experimental MS data", type=str, required=True)
    args = parser.parse_args()
    
    
    
    # list where the lists of the overlap of your MSdata and library will be stored
    tableSC = []

    # MSdirectory: specify directory with your MSdata you want to process and 
    # determine the subclass. Expected to be in .dat or .txt format with two columns
    # and correct headers m/z and intensity.
    #MSdir = "C:/Users/dbeslic/Desktop/Trypsin"
    MSdir = args.input

    # path to your library folder which f.e. contains specific Mass Lists of each subclass
    # these will be compared to your MSdata to find overlapping peaks
    # Expected to be in .txt format with headers m/z, intensity, sequence
    #libraryDir = "C:/Users/dbeslic/Desktop/Massenlisten_insilico"
    libraryDir = args.library
    
    os.chdir(MSdir)
    
    # create results directory
    Path(MSdir+"/results").mkdir(parents=True, exist_ok=True)
    
    # loop to read all files in a specified directory
    for file in glob.glob("*.dat") + glob.glob("*.txt"):
        print("Processing: ",file)
        start = 0
        end = file.find(".dat")
        if (end == -1): end = file.find(".txt")
        mstable = pd.read_table(file, sep=" ", header = None)    
    
        # remove fileending ".dat" or ".txt" in final table
        file = file[start:end]
    
        # preprocess and peakpicking.
        mstable_processed = processing(mstable)
    
        # findSubclass gives a list with overlap of peaks between your MS file and masslists
        subclassList = findSubclass(mstable_processed, libraryDir)
        subclassList.name = file
        # add that list to tableSC to have a table over your data if you have more than one MSfile
        tableSC.append(subclassList)
        os.chdir(MSdir) 
    
    
        msfile_processed = PeakToClass(mstable_processed, libraryDir, file)
        os.chdir(MSdir+"/results") 
        msfile_processed.to_csv(file + "_results.csv", index = False, header=True, sep = ";")
        # findSubclass changes directory for Comparison with Masslists
        

        #InteractivePlotMS(mstable_processed, file)
        os.chdir(MSdir) 
    # dataframe enables better representation of table than list 
    subclassTable = pd.DataFrame(tableSC)
    subclassTable.to_csv(MSdir + "/results/result_overview.csv", index = False, header=True, sep = ";")
    
if __name__ == "__main__":
    main()

