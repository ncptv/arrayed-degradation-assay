
# Inceptive Arrayed RNA Degradation Assay Analysis Code Documentation

The purpose of this code is to analyse the data for RNA stability assay.

# Requirements

First clone the code repository from git into your machine.

[Install Python](https://docs.python.org/3/using/windows.html), if it isn’t installed there yet.

This script requires several standard computational libraries like `numpy`, `pandas`, `sklearn` and `scipy`, as well as `plotly` and `PIL` for visualizations. You will find all of them in `requirements.txt` file.

In order to install them enter into a directory with your code and run:

```sh
pip install -r requirements.txt
```

# Usage

Main script responsible for running analysis is located in `arrayed_degradation/analyze_experiment.py`.

In order to list all arguments along with description just run:

```sh
python arrayed_degradation/analyze_experiment.py --help
```

## Arguments

Here is a list of most important arguments to run the stability script:

Obligatory arguments:

- `--data_dir_path` path to the directory with input data in a proper format; could be either absolute path or relative with respect to the directory from which you run the script;
- `--min_peak` left bound of the target peak in nucleotides;
- `--max_peak` right bound of the target peak in nucleotides; target peak will be searched within [`min_peak`, `max_peak`] interval;

Optional arguments:

- `--disable_control_peak` whether to use or not the control peak; if you don't have control peaks in your electropherogram please use this argument!
- `--control_peak_min` left bound of the control peak (for example p4p6) in nucleotides; default is 235 nucleotides;
- `--control_peak_max` right bound of the control peak (for example p4p6) in nucleotides; default is 310 nucleotides; control peak will be searched within [`control_min_peak`, `control_max_peak`] interval;
- `--remove_background` whether to remove background from target peak. By default script does not remove background. If this argument is set it will draw a line at the base of detected peak and only the area above that line will be considered for half life calculations. It is not advised to set it because of unstable behaviour and less repeatable results in some cases.
- `--width_param` parameter that controls the width of the target peak. It ranges from 0 to 1. Default is 0.85. The higher the value the wider bounds of the target peak. Do not exeed 0.95 as you will capture everything as peak. 0.75 is safe but quite conservative. 0.85 captures mostly what you expect but can grab some noise.

## Examples

### Example 1 - simple run without control peak and target peak between 1000 and 1500 nts

```sh
python arrayed_degradation/analyze_experiment.py --data_dir_path data --min_peak 1000 --max_peak 1500 --disable_control_peak
```

With the command above you will run the analysis; script will search for the input data within `data` folder (that is a relative path to the directory with the data; but you can also provide absolute path as well); target peak will be searched within **1000 - 1500** nucleotides window on electropherograms. Script **will not use control peak** for normalization.

### Example 2 - run with (1) control peak, (2) target peak background removal and (3) narrower peaks

```sh
python arrayed_degradation/analyze_experiment.py --data_dir_path data --min_peak 1000 --max_peak 1500 --remove_background --width_param 0.7
```

Difference here is that now script will search for control peak and normalize electropherograms traces with respect to them. Peak bounds are not provided so script will use default bounds (235, 310). Additionally, script will remove backgroud from the target peak by drawing a line a base of the peak (keep in mind that it might not work well for certain cases like migrating peaks). Finally, script will try to make peaks bounds more narrow by decreasing `width_param` from default 0.85 into 0.7.

# Input data format

All input data should be placed in the same directory and named following the convention described below (names of CSV files must not be changed anyhow!):

```sh
data_dir
    plate_map.csv
    epg_0.csv
    epg_1.csv
    ...
```

There are two obligatory input files in order to run the script:

- `plate_map.csv` - it should be a CSV file containg well positions of samples being analyzed. It should be universal across replicates, meaning that each of your replicate should be a copy of each other in terms of positioning on a plate. In other words - replicates of the same RNA - timepoint pair must be on a physically different plates at the same position, so that one plate map can describe all of them.

    Here is a simplified example of what a proper plate_map file should look like:

    ```
    001_TP0,002_TP0,003_TP0,001_TP3,002_TP3,003_TP3
    004_TP0,005_TP0,006_TP0,004_TP0,005_TP3,006_TP3
    ```

    In this example in well `A3` we have a sample with label `003_TP0` which represent RNA having ID `003` and run at timepoint `0`h. At well `A6` there is a sample labels ad `003_TP3` which is the same RNA but run at timepoint `3`h. The sample labels should have this format `{rna_id}_TP{timepoint}`, where `rna_id` can be only digits (no letters, no special characters!), and `timepoint` can be intiger or float number (examples: `0`, `1`, `3.5`, `3.333`) that determines timepoints used **in hours**. Do not add any additional prefixes or suffixes as the script might not recognize RNA and timepoint properly.

    Missing wells within plate map are acceptable, so if you don't want to include specifc RNA - timepoint pair in the analysis just remove it from the platemap and leave cell empty. This will affect all replicates though! Do not remove any timepoint=`0`h.

    Do not add any rows or column names in `plate_map.csv`; those will be added internaly by the script during data processing.

    Script will inform you how many unique RNA ids, replicates and timepoints have been detected. If that number is different than expected it will raise a warning and you should check you data integrity. Expected number of samples (and as a consequence electropherograms) is `number of RNA` x `number of replicates` x `number of timepoints`. When can it happen? For instance when you drop some timepoints for a given RNA or when some of the sample labels were not recognized.

- `epg_0.csv`, `epg_1.csv`, ... - these should be a set of Fragment Analyzer electropherogram files (result of capilary electrophoresis); each of them representing one replicate of the same plate map.

    Files needs to be named with this convention `epg_{i}.csv` where i are following numbers indicating replicate number. So for 4 replicates you should have `epg_0.csv`, `epg_1.csv`, `epg_2.csv`, `epg_3.csv`

    Here is a simplified example of what a proper plate_map file should look like:

    ```csv
    Size (nt),"A1: SampA1","A2: SampA2","A3: SampA3"
    1.16,6.93,3.43,8.81
    1.89,9.54,4.29,12.33
    2.62,12.63,6.80,16.83
    ...
    ```

    First column must always be `Size (nt)` indicating nucleotide length obtained from fitting ladder (which is done within Fragment Analyzer software). Other columns names must contain well identifier like `A3` or `F12`.

# Output files

After the analysis is run you will be provided with two result files that will be placed in the same directory as your input files.

- `result.csv` - CSV containing a table with RNA ids, half life, decay rate, standard deviations and other metrics.
- `results.pdf` - PDF containing plots with analysis details like sample traces, decay curve, summary plots etc..
