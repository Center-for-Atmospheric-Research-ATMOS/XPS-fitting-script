# 1. Introduction

The project is built on [lmfit](https://lmfit.github.io/lmfit-py/index.html) package and demonstrates how it can be used to fit XPS spectra.

# 2. Structure of the project

    |- data (input examples)
    |    |- excel (custom excel files)
    |    |- scienta (txt-files from SCIENTA)
    |    |- specs (xy-files from SPECS Prodigy)
    |- fits (fitted spectra in json format)
    |- libs (main code library)
    |    |- __init__.py (to recognize libs as Python library)
    |    |- datafile.py (classes to read data from SPECS and SCIENTA files)
    |    |- fit.py (fitting routines)

# 3. Installation

The current instruction is for Windows users. I assume the Python is already installed. 

It is recommended to create a virtual environment:

```python -m venv <environment_name>```

Then activate the environment:

```<environment_name>\Scripts\activate.bat```

And finally, install the required packages from **requirements.txt**:

```pip install -r requirements.txt```

After everything is installed, with the activated virtual environment type ```jupyter notebook``` to run the Jupyter.

# 4. Tutorials

Please check **[tutorial.ipynb](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-fitting-script/blob/main/tutorial.ipynb)** first to familiarise yourself with the code pipeline.

# 5. Documentation

The main code in **libs** consists of two modules:
1. **datafile** containing classes to parse files containing spectra
2. **fit** containing instruments for fitting

## datafile module

***

**Class name** ```specs_file```

This class represents a parsed xy file from SPECS Prodigy software, commonly used in XPS experiments.
It reads the file, organizes its content by region, and extracts key information like averaged spectra and associated parameters.
Additionally, it offers a method to plot the raw counts for a specific region.

**Attributes:**
* ```path```: (str) Path to the SPECS file.
* ```Eph```: (float, optional) Photon energy of the measurement in eV. Defaults to 0.
* ```body```: (list) Split the file content by region delimiters.
* ```regions```: (dict) Dictionary containing information for each region:
  * ```region_id```: (int) Unique identifier for the region.
  * ```region_name```: (str) Name of the region.
  * ```averaged_spectrum```: (pandas.DataFrame) DataFrame containing the averaged spectrum data points.
    * ```energy```: (Series) Binding energy scale.
    * ```raw_counts```: (Series) Raw counts associated with each energy.
  * ```parameters```: (dict) Dictionary containing extracted parameters from the region details.

**Methods**
*  ```__init__(self, path: str, Eph: float = 0)```:
    *  Initializes the ```specs_file``` object by reading the data file, splitting it by regions, and extracting information for each region.
*  ```plot_averaged(self, region_id)```:
    *  Creates a matplotlib plot for the raw averaged counts of a specified region.

***

**Class name** ```scienta_file```

This class represents a parsed Scienta txt file, commonly used in XPS experiments.
It reads the file, organizes its content by regions, extracts key information like averaged spectra and individual sweeps, and offers methods for plotting and data manipulation.

**Attributes:**
*  ```path```: (str) Path to the .txt format Scienta file containing spectra.
*  ```Eph```: (float, optional) Photon energy of the measurement in eV. Required for kinetic energy scale conversion.
*  ```body```: (list) Split the file content by region delimiters.
*  ```numOfRegions```: (int) Number of regions identified in the file.
*  ```regions```: (dict) Dictionary containing information for each region:
    *  ```region_id```: (int) Unique identifier for the region.
    *  ```region_name```: (str) Name of the region.
    *  ```energy_scale```: (str) Energy scale used in the region (Binding or Kinetic).
    *  ```averaged_spectrum```: (pandas.DataFrame) DataFrame containing the averaged spectrum data points.
        *  ```energy```: (Series) Binding energy scale.
        *  ```raw_counts```: (Series) Raw counts associated with each energy.
    *  ```2D_spectra```: (list) List of 2D arrays, each representing a sweep of data for the scanned region.

**Methods**
*  ```__init__(self, path: str, Eph: float = 0)```:
    *  Initializes the scienta_file object by reading the data file, splitting it by regions, extracting information for each region, and handling energy scale conversions.
*  ```reduce_dimension(self, region_id, exclude = [])```:
    *  Reduces the 2D spectra for a region to a 1D averaged spectrum, optionally excluding specific sweeps.
*  ```plot_heatmap(self, region_id)```:
    *  Creates a heatmap plot using Plotly to visualize the 2D spectra for a region.
*  ```plot_averaged(self, region_id)```:
    *  Creates a matplotlib plot for the raw counts of a specified region's averaged spectrum.
 
***

## fit module

***

**Class name** ```fit_results```

This class encapsulates the results of fitting a spectrum using a specified model. It stores key information like the fitted spectrum, residuals, peak parameters, and uncertainties, and offers methods for data export and visualization.

**Attributes:**
*  ```name```: (str) ID of the fitted spectrum.
*  ```bkg```: (array) Scattered background.
*  ```x```: (array) Binding energy scale.
*  ```y_arr```: (array) Raw counts in XPS spectrum (y-axis).
*  ```report```: (str) Text report of the fit results.
*  ```summary```: (dict) Dictionary containing detailed fit summary information.
*  ```dely```: (array) Uncertainties associated with the fit.
*  ```raw_spectrum```: (array) Raw counts of the spectrum before fitting.
*  ```fitted_spectrum```: (array) The fitted spectrum line.
*  ```residual```: (array) The difference between the raw and fitted spectra.
*  ```peaks```: (dict) Dictionary containing extracted peaks and their components.
*  ```displayed_results```: (pandas.DataFrame) DataFrame summarizing peak parameters for easy viewing.

**Methods**
*  ```__init__(self, x_arr, y_arr, bkg, model, params, name, calibration=0)```:
    *  Initializes the fit_results object by storing fitting data, performing the fit, and calculating/organizing results.
*  ```export_json(self, name='')```:
    *  Exports the fit results as a JSON file, optionally specifying a filename.
*  ```plot_results(self, with_background=True)```:
    *  Creates a visualization of the fit results, including raw data, fitted spectrum, residuals, and individual peaks (optional with background).

***

**Class name** ```numpy_array_encoder```

This class is a custom JSON encoder that enables the serialization of NumPy arrays into JSON format. It handles NumPy arrays by converting them to Python lists, which are compatible with JSON encoding.

**Methods**:
*  ```default(self, obj)```:
    *  Called by the JSON encoder to handle object serialization. Checks if the object is a NumPy array, if it's an array, converts it to a Python list. Otherwise, delegates to the default behavior of the base class json.JSONEncoder.

***

**Function name** ```smooth_Savitzky_Golay```

This function performs smoothing (and optionally differentiation) on a given data array using the Savitzky-Golay filter. This filtering technique effectively removes high-frequency noise while preserving the original shape and features of the signal compared to other methods like moving averages.

**Parameters**:
*  `y`: (array_like, shape (N,)) The data array to be smoothed or differentiated.
*  `window_size`: (int) The length of the smoothing window. Must be an odd positive integer.
*  ```order```: (int) The order of the polynomial used for fitting within the window. Must be less than ```window_size``` - 1.
*  ```deriv```: (int, optional) The order of the derivative to compute (default: 0 for smoothing only).

**Returns**:
*  ```ys```: (ndarray, shape (N)) The smoothed or differentiated data array.

**References**:
*  A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of Data by Simplified Least Squares Procedures. Analytical Chemistry, 1964, 36 (8), pp 1627-1639.
*  Numerical Recipes 3rd Edition: The Art of Scientific Computing, W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery, Cambridge University Press ISBN-13: 9780521880688

***

**Function name** ```cut_range```

This function narrows a binding energy range. It allows for selective analysis of desired range of the spectrum.

**Parameters**:
*  ```file```: (scienta_file or specs_file object) The loaded data file containing the spectrum data.
*  ```region```: (int) The identifier of the region within the file from which to extract data.
*  ```cut_range```: (list or None, optional) A two-element list defining the energy range (start and end) to be extracted. If None, the entire region is returned.

**Returns**:
*  ```x```: (array) The extracted binding energy values within the specified range.
*  ```y```: (array) The corresponding raw counts for the extracted energy range.

***

**Function name** ```get_shirley```

This function calculates the Shirley background for a given dataset represented by energy (x) and intensity (y) values. It finds the largest peak and uses the minimum values on either side as anchors to fit a background function. This background is then subtracted from the original data to reveal the true signal.

**Parameters**:
*  ```x```: (array) Array of energy values.
*  ```y```: (array) Array of corresponding intensity values.
*  ```tol```: (float, optional) Tolerance for convergence criterion (default: 1e-5).
*  ```maxit```: (int, optional) Maximum number of iterations (default: 10).
*  ```window_size```: (int, optional) Window size for smoothing with Savitzky-Golay filter (default: 5).

**Returns**:
*  ```S```: (array) The calculated Shirley background function values.

**Note**
*  The source code was taken from [here](https://github.com/kaneod/physics/blob/master/python/specs.py) with minor changes.

***

**Function name** ```get_tougaard```

This function calculates the Tougaard background for a given dataset represented by energy (x) and intensity (y) values. It iteratively fits a background function based on empirical coefficients and energy differences to model the continuum background in XPS or PE spectra.

**Parameters**:
*  ```x```: (np.ndarray) Array of energy values.
*  ```y```: (np.ndarray) Array of corresponding intensity values (counts).
*  ```window_size```: (int, optional) Window size for smoothing with Savitzky-Golay filter (default: 30).
*  ```tb```: (float, optional) Empirical coefficient for background intensity (default: 2866).
*  ```tc```: (float, optional) Empirical coefficient for energy dependence (default: 1643).
*  ```tcd```: (float, optional) Additional coefficient for energy dependence (default: 1).
*  ```td```: (float, optional) Additional coefficient for energy dependence (default: 1).
*  ```maxit```: (int, optional) Maximum number of iterations (default: 10).

**Returns**:
*  ```S```: (np.ndarray) The calculated Tougaard background function values.

***

**Function name** ```baseline_als```

This function performs background correction on a given data array (y) using the Asymmetrically Reweighted Penalized Least Squares (ARPLS) smoothing method. This technique effectively separates the true signal from the underlying baseline in various spectroscopy applications.

**Parameters**:
*  ```y```: (array_like) The data array containing the signal with baseline.
*  ```lam```: (float) Smoothness parameter controlling the balance between fitting and noise reduction. Typically, good values range from 10^2 to 10^9.
*  ```p```: (float, optional) Asymmetry parameter (default: 0). Set between 0.001 and 0.1 for signals with positive peaks, and between 0.9 and 0.999 for negative peaks.

**Returns**:
*  ```z```: (array_like) The background data array.

**Note**
*  The script is based on [this manuscript](doi.org/10.1039/C4AN01061B).

***
