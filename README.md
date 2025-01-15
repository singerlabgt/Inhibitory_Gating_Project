# Code and figures of "Goal specific hippocampal inhibition gates learning". 
## System requirements
The code has been tested on the Windows 11 Operating system. No required non-standard hardware is needed.

## Installation guide
- To run MATLAB-based figures, a MATLAB environment is needed. The scripts have been tested with MATLAB R2023a. 
    - Required Add-On: Statistics and Machine Learning Toolbox
- To run `PlotExtendedDataFigure03_C.ipynb`, a python-based jupyter environment is needed. The code was tested with the conda environment with the following packages
    - python=3.11.5
    - jupyter=1.0.0
    - numpy=1.26.0
    - scipy=1.11.3
    - pandas=2.1.1
    - matplotlib=3.7.2
    - seaborn=0.12.2
    - statannot=0.2.3

    If you are using conda, the conda installation guide can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
    To create the environment, open the terminal and enter 
    ```
    conda create -n PlotExtendedDataFigure03_C python=3.11.5 jupyter numpy=1.26.0 scipy=1.11.3 pandas=2.1.1 matplotlib=3.7.2 seaborn=0.12.2 statannot=0.2.3 -c anaconda -c conda-forge -y
    ```
- To run linear-mixed model, an R environment is required. The scripts have been tested with R version 4.2.2. Based on the package requirements, an R version >= 4.1.0 is recommended. The following packages are required (versions should not affect the result)
    - lme4=1.1.35.1
    - lmerTest=3.1.3
    - emmeans=1.8.9

    To install the packages, use the command 
    ```
    install.packages("{PACKAGE_NAME}", version='{PACKAGE_VERSION}')
    ```

## Demo/Instructions for use
- The source data (doi: 10.6084/m9.figshare.28191596) are shared through [figshare](https://figshare.com/s/d61d3088abbea65dd6ad). Download and unzip `Demo_Data` file and put it at the same directory as Demo_Code. The code scripts will use the data from `Demo_Data`. The spreadsheet `VR_NoveltySpreadsheet.xlsx` contains metadata for each animal and session.
- To generate the MATLAB-based figures, go to `Demo_Code/{FIGURE_INDEX}`, change `maindir` to manuscript folder and run the script file. The output figures will be saved in `Demo_Figures` with the figure index in the name. For each figure the expected run time should be around 1-5 minutes. 
- To run `PlotExtendedDataFigure03_C.ipynb`, activate the conda environment with the command 
    ```
    conda activate PlotExtendedDataFigure03_C
    ```
    Jupyter notebook can be opened with the command
    ```
    jupyter notebook
    ```
    Then run all the cells to generate the figures. The output figures will be saved in `Demo_Figures` with the figure index in the name. 

- To run the linear mixed model stats, set R working directory to a figure folder under `Demo_Code/LMM_R`, run the R script. The resulting stat test files will be saved in the same folder.

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
This code is covered under the Apache 2.0 License. The dataset is covered under the CC BY 4.0 License. 

## GitHub link to Demo_Code and Demo_Figures
https://github.com/singerlabgt/Inhibitory_Gating_Project.git
