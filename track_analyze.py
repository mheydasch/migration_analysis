#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:05:31 2019

@author: max
"""
import sys
import os
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np
import pandas as pd
#%%
sys.path.append(os.path.realpath(__file__))
import track_import 

def parseArguments():
  # Define the parser and read arguments
  parser = argparse.ArgumentParser(description='analyze migration data)
  parser.add_argument('-d', '--dir', type=str, help='the experiment directory, e.g. /Volumes/imaging.data/Max/REF52/beta_pix/pix_4/', required=True)
  args = parser.parse_args()
  return(args)

#%% calculating distance traveled
def calc_dist():
    l=[]
    #looping through each unique cell
    for n in data.cells:
        #taking only the data of that cell
        single_track=data.tracks[data.tracks['unique_id']==n]
        #sorting the data by time
        single_track=single_track.sort_values(by=['Metadata_Timepoint'])
        #calculating x and y step sizes over the whole timecourse
        xd=single_track['Location_Center_X'][0:].reset_index()-single_track['Location_Center_X'][1:].reset_index()
        yd=single_track['Location_Center_Y'][0:].reset_index()-single_track['Location_Center_Y'][1:].reset_index()
        #calculating 2D stepsize over the whole time course
        stepsize=np.sqrt(xd['Location_Center_X']**2+yd['Location_Center_Y']**2)
        #last number wil always be nan, so drop it
        stepsize=stepsize.dropna()
        #summing up stepsize to get cumulative distance
        cumulative_distance=sum(stepsize)
        #calculating the distance between start and end point
        net_distance=np.sqrt((single_track['Location_Center_X'].iloc[-1]-
                                single_track['Location_Center_X'].iloc[0])**2)
        #adding the calculated values to a dictionary
        temp={'unique_id':n, 'cumulative_distance':cumulative_distance, 'net_distance':net_distance}
        #appending the dictionary to a list to keep it over the loops
        l.append(temp)
    #making data frame of the dictionaries    
    distances=pd.DataFrame(l)
    #calculating persistence
    distances['persistence']=distances['net_distance']/distances['cumulative_distance']
    #adding classifier
    distances.loc[distances.unique_id.str.contains('C'), 'Classifier']='Ctrl'
    distances.loc[distances.unique_id.str.contains('B'), 'Classifier']='pix'
    distances=distances.dropna()
    return distances



#%%making migration lineplots

if __name__ == '__main__':
    args=parseArguments()
    path=args.dir
    print(args)
    #initiating data object with given path
    data=track_import.Experiment_data(path)
    #performing loading of data, interpolation and normalization
    data.normalize_tracks()
    #calculate distances travelled for each cell
    distances=calc_dist()
    data.tracks.loc[data.tracks.Metadata_Well.str.contains('C'), 'Classifier']='Ctrl'
    data.tracks.loc[data.tracks.Metadata_Well.str.contains('B'), 'Classifier']='pix'
    #figure dimensions
    dims=(10,8)
    fig, ax=plt.subplots(figsize=dims)
    #creating figure
    migrationplot=sns.lineplot(data=data.tracks, x='Location_Center_X_Zeroed', y='Location_Center_Y_Zeroed', 
                 hue='Classifier', units="unique_id", estimator=None)    
    persistanceplot=sns.violinplot(x='Classifier', y='persistence', data=distances)
    distanceplot=sns.violinplot(x='Classifier', y='cumulative_distance', data=distances)
    #saving figures at same location as given path
    migrationplot.savefig('{}migrationplot.png'.format(path))
    persistanceplot.savefig('{}persistanceplot.png'.format(path))
    distanceplot.savefig('{}distanceplot.png'.format(path))