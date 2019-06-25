#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:34:49 2019

@author: max
"""

import os
import re
import functools
import pandas as pd
import numpy as np

path='/Volumes/imaging.data/Max/REF52/beta_pix/pix_2/'

class Experiment_data:
    def __init__(self, path):
        self.exp_path=path
    def find_csv(self):
        '''
        returns a list with all csv files in the KD folder
        '''
        #pattern: must end with '.csv'
        csv_find=re.compile('tracks\.csv$')
        #finds the directory with that specific name
        find_dir='trackXY'
        self.i_dirs=[]
        #Knockdown pattern, must match inserted variable self.__KD followed by '/'
        #and one or more digits
        for root, dirs, files in os.walk(self.exp_path):
            #looks in each folder from given path
            #if folder is matching the KD pattern, the find_dir pattern and
            #if there are csv files
            if find_dir in root and len([x for x in files if re.search(csv_find, x)])!=0:
                #finds the csv files in the folder and adds them to a list
                csv_files=[x for x in files if re.search(csv_find, x)]
                for csv in csv_files:
                    #each csv file is added to the path and appended to a list
                    self.i_dirs.append(os.path.join(root, csv))            
        return self.i_dirs
   
#distances.loc[distances.unique_id.str.contains('C'), 'Classifier']='Ctrl'

    def load_tracks(self):
        track_list=[]
        for i in self.i_dirs:
            temp=pd.read_csv(i, header=0)
            track_list.append(temp)
        self.tracks = pd.concat(track_list, axis=0, sort=True)
        self.tracks=self.tracks.reset_index()
        for line, i in enumerate(self.tracks['Location_Center_X']):
            self.tracks.loc[line, 'unique_id']=self.tracks.loc[line, 'Metadata_Well']+\
            '_'+str(self.tracks.loc[line, 'Metadata_Site_Num'])+'_'+str(self.tracks.loc[line, 'track_id'])
        return self.tracks
    
    def interpolate_tracks(self):
        #interpolating x
        for line ,i in enumerate(self.tracks['Location_Center_X']):
            #checking if line is nan
            if np.isnan(i):
                #checking if next value is not nan, if the next line is not above df boundaries and if next line is still
                #part of the same cell
                if np.isnan(self.tracks.loc[line+1, 'Location_Center_X'])==False and len(self.tracks)>line+1 \
                    and self.tracks.loc[line+1, 'unique_id']==self.tracks.loc[line, 'unique_id']:
                    upper=self.tracks.loc[line+1, 'Location_Center_X']
                if np.isnan(self.tracks.loc[line-1, 'Location_Center_X'])==False and line-1>=0 \
                and self.tracks.loc[line-1, 'unique_id']==self.tracks.loc[line, 'unique_id']:
                    lower=self.tracks.loc[line-1, 'Location_Center_X']
                self.tracks.loc[line, 'Location_Center_X']=(upper+lower)/2
        #interpolating y
        for line ,i in enumerate(self.tracks['Location_Center_Y']):
            #checking if line is nan
            if np.isnan(i):
                #checking if next value is not nan, if the next line is not above df boundaries and if next line is still
                #part of the same cell
                if np.isnan(self.tracks.loc[line+1, 'Location_Center_Y'])==False and len(self.tracks)>line+1 \
                    and self.tracks.loc[line+1, 'unique_id']==self.tracks.loc[line, 'unique_id']:
                    upper=self.tracks.loc[line+1, 'Location_Center_Y']
                if np.isnan(self.tracks.loc[line-1, 'Location_Center_Y'])==False and line-1>=0 \
                and self.tracks.loc[line-1, 'unique_id']==self.tracks.loc[line, 'unique_id']:
                    lower=self.tracks.loc[line-1, 'Location_Center_Y']
                self.tracks.loc[line, 'Location_Center_Y']=(upper+lower)/2

    def normalize_tracks(self):
        self.find_csv()
        self.load_tracks()
        self.interpolate_tracks()
        #getting a list of the unique ids
        self.cells=self.tracks['unique_id'].unique()
        #initiating the new column for normalized values
        self.tracks['Location_Center_Y_Zeroed']=0
        self.tracks['Location_Center_X_Zeroed']=0
        #looping through the unique ids
        for name in self.cells:
            single_track=self.tracks[self.tracks['unique_id']==name]
            single_track=single_track.sort_values(by=['Metadata_Timepoint'])
            first_index=single_track.index.values[0]
            first_y=single_track.loc[first_index, 'Location_Center_Y']
            first_x=single_track.loc[first_index, 'Location_Center_X']
            single_track['Location_Center_Y_Zeroed']=single_track['Location_Center_Y']-first_y
            single_track['Location_Center_X_Zeroed']=single_track['Location_Center_X']-first_x
            self.tracks['Location_Center_Y_Zeroed'][self.tracks['unique_id']==name]=single_track['Location_Center_Y_Zeroed']
            self.tracks['Location_Center_X_Zeroed'][self.tracks['unique_id']==name]=single_track['Location_Center_X_Zeroed']

        