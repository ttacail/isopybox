#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

ISOPYBOX.py

"""

# Module importation
from __future__ import print_function
import numpy as np
import pandas as pd
import INPUT_reader
import os, sys, getopt
import scipy.integrate as sc_int
import matplotlib.pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import time
from shutil import copyfile
from openpyxl import load_workbook
import pydot

class IsotopicBoxModel():
    """
    This is the main engine that is used to create custom isotopic models
    """
    def __init__(self, input_file, outdir=None):
        print('\n\n     Welcome to ISOPYBOX\n')
        self.input_file = input_file
        if outdir==None:
            cwd = os.getcwd()
            run_dir = 'run_' + time.strftime("%Y-%m-%d_%Hh%Mm%Ss")
            if os.name == 'posix':
                os.makedirs(cwd + '/' + run_dir)
                self.outdir = cwd + '/' + run_dir + '/'
            if os.name == 'nt':
                os.makedirs(cwd + '\\' + run_dir)
                self.outdir = cwd + '\\' + run_dir + '/'
        else:
            cwd = outdir
            run_dir = 'run_' + time.strftime("%Y-%m-%d_%Hh%Mm%Ss")
            if os.name == 'posix':
                os.makedirs(cwd + '/' + run_dir)
                self.outdir = cwd + '/' + run_dir + '/'
            if os.name == 'nt':
                os.makedirs(cwd + '\\' + run_dir)
                self.outdir = cwd + '\\' + run_dir + '\\'

    def initial_state(self):
        print('\n**************************************')
        print('Reading input file')
        print('**************************************')
        self.INPUT = INPUT_reader.INPUT_file(self.input_file)
        self.INPUT.read_CONSTS()
        self.INPUT.read_INITIAL()
        self.INPUT.read_FLUXES()
        self.INPUT.read_COEFFS()
        
        # Minimal residence time and length of the time array
        self.tau = np.zeros(np.shape(self.INPUT.boxes_id))
        for bb in range(0, self.INPUT.boxes_nb):
            self.tau[bb] = self.INPUT.boxes_size[bb]/(np.sum(self.INPUT.fluxes[bb]))
        print('\n Residence time')
        print('--------------------------------------')
        for bb in range(0, self.INPUT.boxes_nb):
            print('%s: %.2e' % (self.INPUT.boxes_id[bb], self.tau[bb]))
    
        self.Boxes_size = np.zeros((len(self.INPUT.time),len(self.INPUT.boxes_size)))
        for tt in range(0, len(self.INPUT.time)):
            for ii in range(0, self.INPUT.boxes_nb):
                outflux = 0
                influx = 0
                for jj in range(0, self.INPUT.boxes_nb):
                    outflux += self.INPUT.fluxes[ii][jj]
                    influx += self.INPUT.fluxes[jj][ii]
                self.Boxes_size[tt][ii] = (influx - outflux)*tt + self.INPUT.boxes_size[ii]

    def compute_evolution(self, Delta, func=None):
        print('\n**************************************')
        print('Computing evolution')
        print('**************************************')
        if func is None:
            func = self.evol_ratio
        Ratio = [(delta / 1e3 + 1e0) * self.INPUT.ratio_standard for delta in Delta]
        Ratio = sc_int.odeint(func, Ratio, self.INPUT.time)   # contains time dimension
        Delta = ((Ratio / self.INPUT.ratio_standard) - 1.0) * 1000  # contains time dimension
        return Delta

    def evol_ratio(self, ratio, t):
        """ The evolution function that is used for isotopic ratio evolution"""
        rationew = np.zeros(ratio.size)
        
        # Element mass evolution       
        bsizenew = np.zeros(self.INPUT.boxes_size.size)
        for ii in range(bsizenew.size):
            outflux = 0
            influx = 0
            for jj in range(bsizenew.size):
                outflux += self.INPUT.fluxes[ii][jj]
                influx += self.INPUT.fluxes[jj][ii]
            bsizenew[ii] = (influx - outflux)*t + self.INPUT.boxes_size[ii]
            
        # Ratio evolution
        for ii in range(ratio.size):                                                                            # PREVIOUS LOCAL CALCULATION OF BOX SIZES EVOLUTION IS MISSING
            outflux = 0
            influx = 0
            outflux_bsize = 0
            influx_bsize = 0
            for jj in range(ratio.size):
                outflux = outflux + \
                          self.INPUT.fluxes[ii][jj] / bsizenew[ii] * \
                          self.INPUT.coeffs[ii][jj] * ratio[ii] - \
                          self.INPUT.fluxes[ii][jj] / bsizenew[ii] * ratio[ii]
                influx = influx + \
                          self.INPUT.fluxes[jj][ii] / bsizenew[ii] * \
                          self.INPUT.coeffs[jj][ii] * ratio[jj] - \
                          self.INPUT.fluxes[jj][ii] / bsizenew[ii] * ratio[ii]
                outflux_bsize += self.INPUT.fluxes[ii][jj]
                influx_bsize += self.INPUT.fluxes[jj][ii]
            rationew[ii] = influx - outflux
        return rationew
    
    def final_state(self, Delta_final, outdir=None, wr_evo=None):
        """ """
        print('\n**************************************')
        print('Writing evolution and final state')
        print('**************************************')
        if outdir is None:
            outdir = self.outdir
            
        # -------Writing evolution for delta -------
        if wr_evo != None:
            f = open(self.outdir + self.INPUT.prefix + 'Delta_t.txt', 'w')
            line1 = '# Time'
            for bb in range(0, len(self.INPUT.boxes_id)):
                line1 += '\t%s' % str(self.INPUT.boxes_id[bb])
            line1 += '\n'
            f.write(line1)
            for ii in range(0, len(self.INPUT.time)):
                line = '%s' % str(self.INPUT.time[ii])
                for bb in range(0, len(self.INPUT.boxes_id)):
                    line += '\t%s' % str(self.Delta[ii][bb])
                line += '\n'
                f.write(line)
            f.close()
        # -----------------------------------------
        
        # -------Writing evolution for size -------
        if wr_evo != None:
            f = open(self.outdir + self.INPUT.prefix + 'Size_t.txt', 'w')
            line1 = '# Time'
            for bb in range(0, len(self.INPUT.boxes_id)):
                line1 += '\t%s' % str(self.INPUT.boxes_id[bb])
            line1 += '\n'
            f.write(line1)
            for ii in range(0, len(self.INPUT.time)):
                line = '%s' % str(self.INPUT.time[ii])
                for bb in range(0, len(self.INPUT.boxes_id)):
                    line += '\t%s' % str(self.Boxes_size[ii][bb])
                line += '\n'
                f.write(line)
            f.close()
        # -----------------------------------------
        
        # -------Plotting evolution (delta) --------
        fig = plt.figure(1, (13,5))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.set_xlabel('Time')
        ax1.set_ylabel('$\delta^{%s}$%s' % (str(self.INPUT.numerator),str(self.INPUT.element)) )
        ax2.set_xlabel('Time')
        ax2.set_ylabel('$\delta^{%s}$%s' % (str(self.INPUT.numerator),str(self.INPUT.element)) )
        
        values = range(len(self.INPUT.boxes_id))
        color_map = plt.get_cmap('rainbow') #mod AH: passé de 'spectral' à 'rainbow' 
        cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=color_map)
        
        for bb in range(0, len(self.INPUT.boxes_id)):
            colorVal = scalarMap.to_rgba(values[bb])   
            ax1.plot(self.INPUT.time, self.Delta[:, bb], '-', color=colorVal, lw='1.5', 
                    label='%s' % str(self.INPUT.boxes_id[bb]))
            ax2.semilogx(self.INPUT.time, self.Delta[:, bb], '-', color=colorVal, lw='1.5', 
                    label='%s' % str(self.INPUT.boxes_id[bb]))
        ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        ax1.grid()
        ax2.grid()
        fig.subplots_adjust(right=0.8)
        plt.savefig(self.outdir + self.INPUT.prefix + 'Delta_t.pdf', dpi=200)
        plt.close()
        # ------------------------------
        
        # -------Plotting evolution (size)--------
        fig = plt.figure(2, (13,5))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Size of X')
        ax2.set_xlabel('Time')
        ax2.set_ylabel('Size of X')
        
       
        values = range(len(self.INPUT.boxes_id))
        color_map = plt.get_cmap('rainbow') 
        cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=color_map)
        
        for bb in range(0, len(self.INPUT.boxes_id)):
            colorVal = scalarMap.to_rgba(values[bb])   
            ax1.plot(self.INPUT.time, self.Boxes_size[:, bb], '-', color=colorVal, lw='1.5', 
                    label='%s' % str(self.INPUT.boxes_id[bb]))
            ax2.loglog(self.INPUT.time, self.Boxes_size[:, bb], '-', color=colorVal, lw='1.5', 
                    label='%s' % str(self.INPUT.boxes_id[bb]))
        
        ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#        ax1.set_ylim(-2, 0.2)
#        ax2.set_ylim(-2, 0.2)
        ax1.grid()
        ax2.grid()
        fig.subplots_adjust(right=0.8)
        plt.savefig(self.outdir + '/' + self.INPUT.prefix + 'Size_t.pdf', dpi=200)
        plt.close()
        # ------------------------------
        
        # -------Writing final state in _OUTPUT.xlsx-------   
        self.output_excel = self.outdir + self.INPUT.prefix + '_OUTPUT.xlsx'
        copyfile(self.INPUT.namefile, self.output_excel)
        
        book = load_workbook(self.output_excel)  
        writer = pd.ExcelWriter(self.output_excel, engine='openpyxl')
        writer.book = book
        writer.sheets = dict((ws.title, ws) for ws in book.worksheets)        
        
        df = pd.DataFrame(index=range(0,len(self.INPUT.boxes_id)), 
                          columns=['BOXES_ID', 'SIZE_FINAL', 'DELTA_FINAL'])
        df.BOXES_ID = self.INPUT.boxes_id  
        df.SIZE_FINAL = self.INPUT.boxes_size # to change if time dependent
        df.DELTA_FINAL = self.Delta[-1, :]
        df.to_excel(writer,'FINAL')
        writer.save()
        # ---------------------------------

       
    def plot_state(self, boxes, deltas, outdir=None):
        """ Make a graph of a given state """
        
        # Colors
        values = range(len(self.INPUT.boxes_id))
        color_map = plt.get_cmap('rainbow') 
        cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=color_map)
        colorList = []
        for c in range(0, len(self.INPUT.boxes_id)):
            colorList.append(scalarMap.to_rgba(values[c]))   

        # Graph type
        gr = pydot.Dot(graph_type='digraph')
        
        i_box = 0
        for box in boxes:
            gr.add_node(pydot.Node(box))
            i_box += 1
        
        # Creation of a dictionary
        self.dico_fluxes_f = {}
        for bb in range(0, len(self.INPUT.boxes_id)):
            self.dico_fluxes_f[self.INPUT.boxes_id[bb]] = {}
            for cc in range(0, len(self.INPUT.boxes_id)):
                self.dico_fluxes_f[self.INPUT.boxes_id[bb]][self.INPUT.boxes_id[cc]] = self.INPUT.fluxes[bb][cc]
      
        for i, origin_state in enumerate(self.INPUT.boxes_id):
            for j, destination_state in enumerate(self.INPUT.boxes_id):
                flux = self.INPUT.fluxes[i][j]  
                if flux > 0:
                    gr.add_edge(pydot.Edge(origin_state, destination_state, label="{:.02f}".format(flux)))

        if outdir is None:
            outdir = self.outdir
        gr.write_png(outdir + self.INPUT.prefix + 'final_state.png')
        gr.write_pdf(outdir + self.INPUT.prefix + 'final_state.pdf')
       
    def run(self, plot_graph):
        """ Run the model """
        self.initial_state()
        self.Delta = self.compute_evolution(self.INPUT.delta)
        self.final_state(self.Delta[-1, :], wr_evo=1)
        if plot_graph is not None:
            self.plot_state(self.INPUT.boxes_id, self.Delta[-1, :])
        print('\n**************************************')
        print('Done')
        print('**************************************')
    
#------------------------------------ Main ---------------------------------------  
        
def main(argv):
    out_dir = None
    plot_graph = None
    try:
        opts, arg = getopt.getopt(argv,"hi:o:g:",["iExcelfilename","ooutput_dir","gplot_graph"])
    except getopt.GetoptError:
        print('ISOPYBOX.py -i <Excelfilename> -o <output_dir> -g <plot_graph>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('Usage: python ISOPYBOX.py -i <INPUT_Excelfilename> -o <output_dir> -g <plot_graph>')
            print('ISOPYBOX needs both INPUT_Excelfilename and output_dir arguments to work, no default value for that.')
            print('INPUT_Excelfilename is the path of the name_INPUT.xlsx file')
            print('output_dir is the path of the output directory in which ISOPYBOX will generate the outputs')
            print('-g is optional: -g 1 for the box graph')
            sys.exit()
        elif opt in ("-i", "--iExcelfilename"):
            INPUT_XL_file = str(arg)
            print('ISOPYBOX will use the file ',INPUT_XL_file,' for input')
        elif opt in ("-o", "--ooutput_dir"):
            out_dir = str(arg)
        elif opt in ("-g", "--gplot_graph"):
            plot_graph = str(arg)

    if (os.path.isfile(INPUT_XL_file)): 
        model = IsotopicBoxModel(INPUT_XL_file, outdir=out_dir)
        model.run(plot_graph)
    else:
        print('Excel file ',INPUT_XL_file,'does not exist')
        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])            
