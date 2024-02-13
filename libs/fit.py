# ToDo:
# 1. Variable names are not pytonic

import json
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve
from math import factorial


def smooth_Savitzky_Golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
    the values of the time history of the signal.
    window_size : int
    the length of the window. Must be an odd integer number.
    order : int
    the order of the polynomial used in the filtering.
    Must be less then `window_size` - 1.
    deriv: int
    the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
    the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
      Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    try:
        window_size = np.abs(int(window_size))
        order = np.abs(int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
    
    
def cut_range(file, region, cut_range=None):
    '''
    cut the x and y according to the interval given in cut_range
    return cut x and y
    '''
    x = np.array(file.regions[region]['averaged_spectrum']['energy'])
    y = np.array(file.regions[region]['averaged_spectrum']['raw_counts'])
    
    if cut_range:
        # find indexes of the nearest values
        indx0 = (np.abs(x - cut_range[0])).argmin()
        indx1 = (np.abs(x - cut_range[1])).argmin()

        # check if inverted and cut x and y
        if x[0] > x[-1]:
            x = x[indx1:indx0]
            y = y[indx1:indx0]
        else:
            x = x[indx0:indx1]
            y = y[indx0:indx1]
        
    return x, y
    
    
def get_shirley(x, y, tol=1e-5, maxit=10, window_size=5):
    """
    Calculate the best auto-Shirley background S for a dataset (x,y). Finds the biggest peak
    and then uses the minimum value either side of this peak as the terminal points of the
    Shirley background.
    The tolerance sets the convergence criterion, maxit sets the maximum number
    of iterations.

    ref.: https://github.com/kaneod/physics/blob/master/python/specs.py
    """

    # Make sure we've been passed arrays and not lists.
    x = np.array(x)
    y = smooth_Savitzky_Golay(np.array(y), window_size=window_size, order=3, deriv=0, rate=1)
    
    # Sanity check: Do we actually have data to process here?
    if not (x.any() and y.any()):
        print ("specs.shirley_calculate: One of the arrays x or y is empty. Returning zero background.")
        return np.zeros(x.shape)

    # Next ensure the energy values are *decreasing* in the array,
    # if not, reverse them.
    if x[0] < x[-1]:
        is_reversed = True
        x = x[::-1]
        y = y[::-1]
    else:
        is_reversed = False

    # Locate the biggest peak.
    maxidx = abs(y - np.amax(y)).argmin()

    # It's possible that maxidx will be 0 or -1. If that is the case,
    # we can't use this algorithm, we return a zero background.
    if maxidx == 0 or maxidx >= len(y) - 1:
        print('specs.shirley_calculate: Boundaries too high for algorithm: returning a zero background.')
        return np.zeros(x.shape)

    # Locate the minima either side of maxidx.
    lmidx = abs(y[0:maxidx] - np.amin(y[0:maxidx])).argmin()
    rmidx = abs(y[maxidx:] - np.amin(y[maxidx:])).argmin() + maxidx
    xl = x[lmidx]
    yl = y[lmidx]
    xr = x[rmidx]
    yr = y[rmidx]

    # Max integration index
    imax = rmidx - 1

    # Initial value of the background shape B. The total background S = yr + B,
    # and B is equal to (yl - yr) below lmidx and initially zero above.
    B = np.zeros(x.shape)
    B[:lmidx] = yl - yr
    Bnew = B.copy()

    it = 0
    while it < maxit:
        # Calculate new k = (yl - yr) / (int_(xl)^(xr) J(x') - yr - B(x') dx')
        ksum = 0.0
        for i in range(lmidx, imax):
            ksum += (x[i] - x[i + 1]) * 0.5 * (y[i] + y[i + 1]
                                               - 2 * yr - B[i] - B[i + 1])
        k = (yl - yr) / ksum
        # Calculate new B
        for i in range(lmidx, rmidx):
            ysum = 0.0
            for j in range(i, imax):
                ysum += (x[j] - x[j + 1]) * 0.5 * (y[j] +
                                                   y[j + 1] - 2 * yr - B[j] - B[j + 1])
            Bnew[i] = k * ysum
        # If Bnew is close to B, exit.
        if np.linalg.norm(Bnew - B) < tol:
            B = Bnew.copy()
            break
        else:
            B = Bnew.copy()
        it += 1

#         if it >= maxit:
#             print("specs.shirley_calculate: Max iterations exceeded before convergence.")
    if is_reversed:
        return (yr + B)[::-1]
    else:
        return yr + B
        
        
def get_tougaard(x: np.ndarray, y: np.ndarray, window_size = 30, tb=2866, tc=1643, tcd = 1, td=1, maxit=10):
    """
    Calculate the best auto-Tougaard background for a dataset (x,y). Maxit sets
    the maximum number of iterations, tb and tc are imperical coefficients.
    """
    x = np.array(x)
    y = smooth_Savitzky_Golay(np.array(y), window_size=window_size, order=3, deriv=0, rate=1)
    # Sanity check: Do we actually have data to process here?
    if not (any(x) and any(y)):
        print("One of the arrays x or y is empty. Returning zero background.")
        return np.zeros(x.shape)

	# KE in XPS or PE in XAS
    if x[0] < x[-1]:
        is_reversed = True
        x = x[::-1]
        y = y[::-1]
    else:
        is_reversed = False

    Btou = np.zeros(x.shape)

    it = 0
    while it < maxit:
        for i in range(len(y)-1, -1, -1):
            Bint = 0
            for j in range(len(y)-1, i-1, -1):
                Bint += (y[j] - y[len(y)-1]) * (x[0] - x[1]) * (x[i] - x[j]
                ) / ((tc + tcd * (x[i] - x[j])**2)**2 + td * (x[i] - x[j])**2)
            Btou[i] = Bint * tb

        Boffset = Btou[0] - (y[0] - y[len(y)-1])
        if abs(Boffset) < (0.000001 * Btou[0]) or maxit == 1:
            break
        else:
            tb = tb - (Boffset/Btou[0]) * tb * 0.5
        it += 1

    if is_reversed:
        return (y[len(y) - 1] + Btou)[::-1]
    else:
        return y[len(y) - 1] + Btou
    

def baseline_als(y, lam, p=0, niter=10):
    '''
    Baseline correction using asymmetrically reweighted penalized least squares smoothing
    based on doi.org/10.1039/C4AN01061B
    
    Params:
    p: assymetry, 0 by default. 0.001 ≤ p ≤ 0.1 is a good choice (for a signal with positive peaks)
    lam: smoothness. Should be approx 10^2 ≤ λ ≤ 10^9 
    
    Returns: baseline z 
    '''
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w) # Do not create a new matrix, just update diagonal values
        Z = W + D
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z


class numpy_array_encoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


class fit_results(object):
    def __init__(self, x_arr, y_arr, bkg, model, params, name, calibration=0):
        """
        """
        self.name = name               # ID of current spectrum
        self.bkg = bkg                 # save baseline
        self.x = x_arr + calibration   # energy calibration
        y_arr -= bkg                   # subtract background

        # fit spectrum
        output = model.fit(
            y_arr,
            params,
            x=self.x)
        
        # save report to check fit results if needed
        self.report = output.fit_report(
            sort_pars=True,
            show_correl=False).split("[[Variables]]")[-1]
        
        self.summary = output.summary()                 # save dict with the whole meta data
        self.dely = output.eval_uncertainty(sigma=3)    # evaluate uncertainties
        self.raw_spectrum = output.data                 # spectrum before fitting
        self.fitted_spectrum = output.best_fit          # spectrum after fitting
        self.residual = output.residual                 # residuals
        self.peaks = output.eval_components()           # save peaks separately

        # want to print fitting results as a table. get them first
        # best_values = output.best_values
        best_values = {i[0]: i[1] for i in self.summary['params']}

        # add missing values to output
        for key in self.peaks.keys():
            best_values[key + 'area'] = sum(self.peaks[key])

            # sometimes fitting error can not be calculated
            try:
                best_values[key + 'area_err'] = sum(output.dely_comps[key])
            except AttributeError:
                best_values[key + 'area_err'] = 0
                 # TODO find universal way to evaluate fitting error
                 
            # add position err
            best_values[key + 'position_err'] = output.params[key + 'center'].stderr
        
        # get peaks and their parameters
        unique_peaks = list(set([i.split('_')[0] for i in best_values.keys()]))

        # parameters to display
        unique_parameters = ['center', 'position_err', 'area', 'area_err', 'height','fwhm']
        
        ## fill the table with final results
        ## TODO it can be done for efficiently
        # create dataframe results table
        results_to_display = pd.DataFrame(
            columns=unique_parameters, 
            index=unique_peaks)
        
        # convert best velues to dataframe too
        best_values = pd.DataFrame(best_values, index=['value'])

        # fill in results_to_display from best_values
        for peak in unique_peaks:
            for parameter in unique_parameters:
                results_to_display[parameter][peak] = best_values[peak + '_' + parameter].loc['value']
        
        # save table to attribute
        self.displayed_results = results_to_display.sort_values(by='center')

        # display table, round numbers
        display(self.displayed_results.astype(float).round(2))
    
            
    def export_json(self, name=''):
        """
        Parameters:
            path:str
            path to the folder where the json file will be stored
            name:str
            name of the file if needs to be specified
        """
        data = self.__dict__
        data['displayed_results'] = data['displayed_results'].to_dict()

        with open('fits/' + name + '.json', 'w') as f:    
            json.dump(data, f, cls=numpy_array_encoder)
    
    
    def plot_results(self, with_background=True):
        """
        
        """
        if with_background:
            background = self.bkg
        else:
            background = 0.75 * min(self.bkg)
        
        fig = make_subplots(rows=2, cols=1, row_heights=[0.2, 0.8], subplot_titles=("Residuals", "Fitted spectrum"), vertical_spacing = 0.1)
        
        ## residuals baseline
        fig.add_trace(
            go.Scatter(
                x=self.x,
                y=np.zeros(len(self.x)),
                mode='lines',
                showlegend=False,
                line=dict(color='black')
            ),
            row=1,
            col=1
        )
        ## residuals
        fig.add_trace(
            go.Scatter(
                x=self.x,
                y=self.residual,
                mode='markers',
                name='Residuals',
                showlegend=False,
                line=dict(color='red', width=0.5)),
            row=1,
            col=1
        )
        ## raw spectrum
        fig.add_trace(
            go.Scatter(
                x=self.x,
                y=self.raw_spectrum + self.bkg,
                mode='markers',
                name='Raw spectrum',
                line=dict(color='red', width=3)
            ),
            row=2,
            col=1
        )
        ## fitted spectrum
        fig.add_trace(
            go.Scatter(
                x=self.x,
                y=self.fitted_spectrum + self.bkg,
                mode='lines',
                name='Fitted spectrum',
                line=dict(color='black', width=2)
            ),
            row=2,
            col=1
        )
        ## plot peaks
        for key in self.peaks.keys():
            fig.add_trace(
                go.Scatter(
                    x=self.x,
                    y=self.peaks[key] + background,
                    mode='lines',
                    name=key,
                    line=dict(width=2)
                ),
                row=2,
                col=1
            )
        ## plot background
        fig.add_trace(
            go.Scatter(
                x=self.x,
                y=self.bkg,
                mode='lines',
                name='Baseline',
                line=dict(color='black', width=3, dash='dash')
            ),
            row=2,
            col=1
        )

        fig.update_layout(margin=dict(l=5, r=5, t=20, b=5), height=400, width=800)
            
        fig.update_xaxes(range=[min(self.x), max(self.x)])
        fig.show()
        return fig