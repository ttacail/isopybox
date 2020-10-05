#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Reading Excel input file for ISOPYBOX.py

This file has 4 sheets (CONSTS/INITIAL/FLUXES/COEFS)

"""

# Module importation
from __future__ import print_function
import numpy as np
import pandas as pd
import os


class INPUT_file():
    
    # Initialization 
    def __init__(self, pnamefile):
        occ1 = str(pnamefile).find('INPUT.xlsx') # allow to accept both xls and xlsx files
        if occ1 != -1:
            self.namefile = str(pnamefile)
            if os.name == 'posix':
                occ2 = max([x for x, v in enumerate(self.namefile) if v == '/'])+1  # for Unix
            if os.name == 'nt':
                occ2 = 0#max([x for x, v in enumerate(self.namefile) if v == '\\'])+1  # for Windows, mod AH: fixé à 0 car incapable de découper correctement, implique d'avoir le fichier INPUT dans le même dossier que le script Isopybox
            self.prefix = self.namefile[occ2:occ1]
        else:
            print('\nWrong file or check file name\n')
    
    # Methods
    def read_CONSTS(self):
        consts_f = pd.read_excel(self.namefile, 'CONSTS')
        self.element = consts_f['CONSTS'][0]
        self.numerator = consts_f['CONSTS'][1]
        self.denominator = consts_f['CONSTS'][2]
        self.ratio_standard = consts_f['CONSTS'][3]
        self.time_max = consts_f['CONSTS'][4]
        self.nb_steps = consts_f['CONSTS'][5]
        self.time = np.linspace(0, self.time_max, self.nb_steps)
        print('\nConstants')
        print('--------------------------------------')
        print('Element: ', self.element)
        print('Numerator isotope: ', self.numerator)
        print('Denominator isotope: ', self.denominator)
        print('Standard ratio: ', self.ratio_standard)
        print('Maximum time: ', self.time_max)
        print('Number of time steps: ', self.nb_steps)
        print('--------------------------------------')
    
    def read_INITIAL(self):
        initial_f = pd.read_excel(self.namefile, 'INITIAL')
        self.boxes_id = np.asarray(initial_f['BOXES_ID'], np.dtype(str))
        self.boxes_nb = len(self.boxes_id)
        self.boxes_size = np.asarray(initial_f['SIZE_INIT'], np.dtype(float))
        self.delta = np.asarray(initial_f['DELTA_INIT'], np.dtype(float))
#        print('\nInitial')
#        print('--------------------------------------')
#        print(initial_f)
#        print('Number of boxes: ', self.boxes_nb)
#        print('Boxes: ', self.boxes_id)
#        print('Masses: ', self.boxes_size)
#        print('Deltas: ', self.delta)
#        print('--------------------------------------')
        
    def read_FLUXES(self):
        fluxes_f = pd.read_excel(self.namefile, 'FLUXES')
        temp = np.asarray(fluxes_f)
        self.fluxes = np.zeros((self.boxes_nb, self.boxes_nb))
        for ii in range(0, self.boxes_nb):
            for jj in range(0, self.boxes_nb):
                self.fluxes[ii][jj] = temp[ii][jj+2] # peut etre +1 si .xlsx avec mac... #mod AH: passé de +1 à +2 pour lecture correct de l'input file
#        print('\nFluxes')
#        print('--------------------------------------')
#        print('Fluxes: ', self.fluxes)
#        print('--------------------------------------')
    
    def read_COEFFS(self):
        coeffs_f = pd.read_excel(self.namefile, 'COEFFS')
        temp = np.asarray(coeffs_f)
        self.coeffs = np.zeros((self.boxes_nb, self.boxes_nb))
        for ii in range(0, self.boxes_nb):
            for jj in range(0, self.boxes_nb):
                self.coeffs[ii][jj] = temp[ii][jj+2] #mod AH: passé de +1 à +2 pour lecture correct de l'input file
#        print('\nCoeffs')
#        print('--------------------------------------')
#        print('Coeffs: ', self.coeffs)
#        print('--------------------------------------')