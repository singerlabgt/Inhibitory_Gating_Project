# Code and figures of "Goal specific hippocampal inhibition gates learning". 
## System requirements
The code has been tested on the Windows 11 Operating system. No required non-standard hardware is needed.

## Installation guide
### 1. MATLAB Environment (for MATLAB-based figures)
- To run MATLAB-based figures, a MATLAB environment is needed. The scripts have been tested with **MATLAB R2023a**. 
    - **Required Add-On**: Statistics and Machine Learning Toolbox
    - **Estimated installation time**: The MATLAB installation process usually takes 5-30 minutes depending on internet speed and system resources.
### 2. Python Environment (for `PlotExtendedDataFigure03_C.ipynb` notebook)
- To run `PlotExtendedDataFigure03_C.ipynb`, a python-based jupyter environment is needed. The code was tested with a **Conda environment** with the following packages
    - python=3.11.5
    - jupyter=1.0.0
    - numpy=1.26.0
    - scipy=1.11.3
    - pandas=2.1.1
    - matplotlib=3.7.2
    - seaborn=0.12.2
    - statannot=0.2.3

    If you are using **Conda**, follow these steps:
    1. **Install Conda** (if not already installed). You can find the installation [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
    2. **Create the Conda environment** by running the following command in your terminal:
        ```
        conda create -n PlotExtendedDataFigure03_C python=3.11.5 jupyter numpy=1.26.0 scipy=1.11.3 pandas=2.1.1 matplotlib=3.7.2 seaborn=0.12.2 statannot=0.2.3 -c anaconda -c conda-forge -y
        ```
        **Estimated Installation Time:** This environment setup usually takes 5-10 minutes, depending on your internet speed and system resources.
    
    
### 3. R Environment (for Linear Mixed Model)
- To run linear-mixed model, an **R environment** is required. The scripts have been tested with **R version 4.2.2**. Based on the package requirements, an R version **>= 4.1.0** is recommended. The following packages are required (versions should not affect the result)
    - lme4 (v1.1.35.1)
    - lmerTest (v3.1.3)
    - emmeans (v1.8.9)

    To install the packages, run the following command in R:
    ```
    install.packages(c("lme4", "lmerTest", "emmeans"))
    ```
    **Estimated installation time:**: Installing each R package typically takes 1-2 minutes, though it may take longer based on internet speed and system resources.

## Demo/Instructions for use
### 1. Data Setup:  
- Download source data (doi: 10.6084/m9.figshare.28191596) are shared through [figshare](https://figshare.com/s/d61d3088abbea65dd6ad). 
- Unzip `Demo_Data` compressed file and put it at the same directory as `Demo_Code`. The code scripts will use the data from `Demo_Data`. The spreadsheet `VR_NoveltySpreadsheet.xlsx` contains metadata for each animal and session.
### 2. MATLAB-based figures:
- Navigate to `Demo_Code/{FIGURE_INDEX}`  
- Change `maindir` to manuscript folder 
- Run the script file. The output figures will be saved in `Demo_Figures` with the figure index in the name. 
  - **Estimated run time:**: For each figure the expected run time should be around 1-5 minutes. 
### 3. Jupyter notebook (for `PlotExtendedDataFigure03_C.ipynb`)
-  **Activate the environment:**
    ```
    conda activate PlotExtendedDataFigure03_C
    ```
- **Launch Jupyter notebook:**
    ```
    jupyter notebook
    ```
    Then, open and run all cells in `PlotExtendedDataFigure03_C.ipynb` to generate the figures. Output figures will be saved in the `Demo_Figures` folder with the figure index in the name.
### 4. Linear Mixed Model (LMM) in R
- Set R working directory to a figure folder under `Demo_Code/LMM_R`
- Run the R script. The resulting stat test files will be saved in the same folder.

## Intermediate files generation
The folder `Demo_Code/Helpers` contains scripts to generate key intermediate files in `Demo_Data`.
- `script_DecodeRippleContent.m`: related with `RipData_250ms.mat`. 
    - This script uses a single time window per sharp-wave ripple event to
    decode the most likely information the SWR has about the distance to the
    closest reward zone. Outputs a giant 'RipData' output structure with all
    ripple events and the bins with the highest spatial probability.
    This script requires access to the Singer Lab ProcessedData folder on
    server for the extracted 'ripples' and 'bestRippleChan' matfiles. Also
    requires 'placecodingoutput' matfile from script_PlaceCodingProperties.m

- `script_MultipleLinearRegression.m`: related with `allsess_raw_vs_residuals` data files. 
    - This script takes binned firing rate maps (binned based on distance or
    time) and performs multiple linear regression using binned speed and lick
    rate behavioral maps as two predictors. The main purpose of the script is
    to regress out potential effects of position/time-related changes in
    behavior on the observed firing rates, and to create a residual firing
    rate map that cannot be attributed to simple changes in behavior. The
    outputs are 1) trial-by-trial residual firing rates of the same size as
    original binned rate map structure, and 2) averaged residual firing rates
    across trials per unit, using at least 5 correct trials only. 
    
- `script_PlaceCodingProperties.m`: related with `placecodingout_09212024.mat`. 
    - This script takes the raw binned firing rate maps and shuffled rate maps to calculate spatial information and identify "spatially modulated" units. To determine spatially modulated cells, used 95th percentile or above the spatial information distribution from shuffled data. Takes a long time to run all 1000 shuffled data. 


## License
This code is covered under the MIT License. The dataset is covered under the CC BY 4.0 License. 

## GitHub link to Demo_Code and Demo_Figures
https://github.com/singerlabgt/Inhibitory_Gating_Project.git
