import re
import numpy as np
from math import factorial
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt


# For specs files
params_to_parse = {
    'acquisition_date': '# Acquisition Date:             ',
    'analyzer_slit':    '# Analyzer Slit:                ',
    'curves_scan':      '# Curves/Scan:                  ',
    'values_curve':     '# Values/Curve:                 ',
    'dwell_time':       '# Dwell Time:                   ',
    'Eph':              '# Excitation Energy:            ',
    'pass_energy':      '# Pass Energy:                  ',
    'bias_voltage':     '# Bias Voltage:                 ',
    'detector_voltage': '# Detector Voltage:             ',
    'eff_workfunction': '# Eff. Workfunction:            ',
    'number_of_scans':  '# Number of Scans: ',
}


class specs_file(object):

    def __init__(self, path: str, Eph:float=0):
        '''
        Input:
        ---
        path: str path to a specs datafile in *.txt format
        Eph: photon energy of the measurement in eV
        '''
        ## convert the whole file to string and split by regions
        with open(path, 'r') as f:
            self.body = re.split(r'# Region:                       ', f.read())[1:]
        
        self.regions = {}
        for spectrum in self.body:           
            region_id = int(re.search(
                r'# Spectrum ID:\s*?(.*?)\n',
                spectrum
            ).group(1))
           
            self.regions[region_id] = {} # create a separate dict for a region
            self.regions[region_id]['region_name'] = spectrum.split('\n')[0]
            
            # get spectrum with background
            self.regions[region_id]['averaged_spectrum'] = pd.read_csv( 
                StringIO(spectrum.split('counts/s\n#\n')[1].replace('#', '').rstrip()),
                sep='\s+',
                names=['energy', 'raw_counts'],
                on_bad_lines='skip'
            )

            self.regions[region_id]['averaged_spectrum'] = self.regions[region_id]['averaged_spectrum'][
                pd.to_numeric(self.regions[region_id]['averaged_spectrum']['raw_counts'], errors='coerce').notnull()].astype('float')
            
            # add the parameters from params_to_parse
            for key, value in params_to_parse.items():
                self.regions[region_id][key] = re.search(
                    r'{}(.*?)\n'.format(value),
                    spectrum
                ).group(1)
                
                # convert digits to float
                try:
                    self.regions[region_id][key] = float(self.regions[region_id][key])
                except ValueError:
                    pass
    
    
    def plot_averaged(self, region_id):
        '''
        plot raw counts for chosen region_id
        '''
        fig, ax = plt.subplots(figsize=(7, 3))
        self.regions[region_id]['averaged_spectrum'].plot(
            x='energy',
            y='raw_counts',
            label='counts',
            ax=ax)
        plt.legend()
        plt.show()
        
        return ax


class scienta_file(object):
    def __init__(self, path: str, Eph: float = 0):
        '''
        ---
        Input:
        ---
        path: str path to the scienta datafile in *.txt
        '''
        ## convert the whole file to string and split by regions
        with open(path, 'r') as f:
            self.body = re.split(r'\[Region \d\]', f.read())
        
        self.numOfRegions = int(re.search(r'Number of Regions=....',    # define number of regions
                                          self.body[0])[0].split('=')[1])
        self.regions = {}
        for region_id in range(1, self.numOfRegions+1):
            self.regions[region_id] = {}
            self.regions[region_id]['region_name'] = re.search(r'Region Name=(.*?)\n',self.body[region_id]).group(1)
            self.regions[region_id]['energy_scale'] = re.search(r'Energy Scale=(.*?)\n',self.body[region_id]).group(1)

            self.regions[region_id]['averaged_spectrum'] = pd.DataFrame()
            
            self.regions[region_id]['averaged_spectrum']['energy'] = [
                float(j) for j in re.search(r'Dimension 1 scale=(.*?)\n',self.body[region_id]).group(1).split(' ')
            ]
            
            if self.regions[region_id]['energy_scale'] == 'Kinetic':
                if Eph==0:
                    print('[ERR] Energy scale is kinetic, please add Eph=...')
                    print(f' Hint: {self.regions[region_id]["region_name"]}')
                    raise Exception(f'Photon energy is required')
                else:
                    self.regions[region_id]['averaged_spectrum']['energy'] = [Eph - i for i in self.regions[region_id]['averaged_spectrum']['energy']]
            
            ## if we have 3 dimensional dataset
            rawSpectra = re.split(r'\n\n\[Data \d:\d*?\]\n ', self.body[region_id])[1:]

            ## check if dimension size less than 3
            if not len(rawSpectra):
                rawSpectra = re.split(r'\n\n\[Data \d\]\n ', self.body[region_id])[1:]
            
            ## convert raw str spectra to float 2D matrix
            self.regions[region_id]['2D_spectra'] = []
            for dimension in rawSpectra:
                ## convert spectra to float and stuck them in one 2D matrix             
                lines = [
                    list(x) for x in zip(*[
                        [
                            float(j) for j in i.split('  ')
                        ] for i in dimension.strip().splitlines()])
                ][1:]
                
                ## append each spectra to matrix
                [self.regions[region_id]['2D_spectra'].append(line) for line in lines]
                
                self.regions[region_id]['averaged_spectrum']['raw_counts'] = np.sum(
                    self.regions[region_id]['2D_spectra'], axis=0
                    ) / len(self.regions[region_id]['2D_spectra'])
            
            del rawSpectra # clear from memory
        
        print(
            f'Found regions: {self.numOfRegions}\n',
            ', '.join([self.regions[region_id]['region_name'] for region_id in range(1, self.numOfRegions+1)]),
            'Energy scale:',
            ', '.join([self.regions[region_id]['energy_scale'] for region_id in range(1, self.numOfRegions+1)])
        )
        
    
    def reduce_dimension(self, region_id, exclude = []):
        """ Averaged spectrum is the raw data ready to be fitted
        it calls reduceDimension because we convert 2D matrix of sweeps to 1D array
        and can be used to update this array excluding some bad sweeps
        -----------
        Parameters:
        -----------
        region: int
        the scanned region to which function will be applied
        exclude: 1D array
        some spectra could be bad and thus excluded for better result.
        --------
        Returns:
        --------
        self.averagedSpectrum: np.array
        averaged spectrum
        """
        arr2D = self.regions[region_id]['2D_spectra'] # get 2D matrix with sweeps
        for i in exclude: # exclude bad spectra
            del arr2D[i]
        ## save spectrum to main class parameters
        self.regions[region_id]['averaged_spectrum']['raw_counts'] = np.sum(arr2D, axis=0) / len(arr2D)
        return self.regions[region_id]['averaged_spectrum']['raw_counts']
    
    
    def plot_heatmap(self, region_id):
        """region - region to plot
        plot heatmap using pixel
        """
        fig = px.imshow(
            self.regions[region_id]['2D_spectra'],
            x = self.regions[region_id]['averaged_spectrum']['energy'],
            width=700,
            height=300,
            title = 'Region: ' + self.regions[region_id]['region_name'],
            labels=dict(x='Binding or Kinetic energy',y='Sweep', color='Intensity'),
            aspect='auto',
        )
        fig.update_layout(
            margin={"t": 30, "b": 0, "r": 0, "l": 0, "pad": 0},
        )
        fig.show(renderer='jupyterlab')


    def plot_averaged(self, region_id):
        '''
        plot raw counts for chosen region_id
        '''
        fig, ax = plt.subplots(figsize=(7, 3))
        self.regions[region_id]['averaged_spectrum'].plot(
            x='energy',
            y='raw_counts',
            label='counts',
            ax=ax)
        plt.legend()
        plt.show()