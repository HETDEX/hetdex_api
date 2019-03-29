# -*- coding: utf-8 -*-
"""

Initiates the Detections class. 

An API may or may not be developed in the future.
We recommend just using pytables since this is only one file.

Created on 2019/01/28

@author: Erin Mentuch Cooper
"""

import tables as tb

import config


class Detections:
    def __init__(self, survey):
        '''
        Initialize the detection catalog class for a given data release

        Input
        -----
        survey : string
            Data release you would like to load, i.e., 'DR1' or 'Parallel'.
            This is case insensitive.
        '''
        survey_options = {'dr1': config.detecth5,
                          'parallel': 'PATHNAME'}
        if survey.lower() not in survey_options:
            print('survey not in survey options')
            print(survey_options)
            return None
        self.filename = survey_options[survey]
        self.hdfile = tb.open_file(self.filename, mode='r')
        colnames = self.hdfile.root.Info.Detections.colnames
        for name in colnames:
            setattr(self, name,
                    getattr(self.hdfile.root.Info.Detections.cols, name)[:])
