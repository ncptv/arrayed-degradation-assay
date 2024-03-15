
# Inceptive Arrayed RNA Degradation Assay Analysis Code Documentation

The purpose of this code is to analyse the data for RNA stability assay.

# Requirements

First clone the code repository from github into your machine via git. If you don't have git you can download it from [here](https://github.com/git-guides/install-git). If you don't want to install it, just download the .zip file containing code for this repo (you can do it by clicking on blue `<>Code` button on main page of this repo, and then `Download ZIP`).

[Install Python 3.10](https://docs.python.org/3.10/using/windows.html), if it isn’t installed there yet.

This script requires several standard computational libraries like `numpy`, `pandas`, `sklearn` and `scipy`, as well as `plotly` and `PIL` for visualizations. You will find all of them in `requirements.txt` file.

In order to install them enter into a directory with your code and install it via pip. Here are commands to do that:

```sh
cd arrayed-degradation-assay
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
- `--time_unit` time unit used in `plate_map.csv` timepoint labels. `m` for minutes, `h` for hours, `d` for days. Defaults to minutes.
- `--remove_background` whether to remove background from target peak. By default script does not remove background. If this argument is set it will draw a line at the base of detected peak and only the area above that line will be considered for half life calculations. It is not advised to set it because of unstable behaviour and less repeatable results in some cases.
- `--width_param` parameter that controls the width of the target peak. It ranges from 0 to 1. Default is 0.85. The higher the value the wider bounds of the target peak. Do not exeed 0.95 as you will capture everything as peak. 0.75 is safe but quite conservative. 0.85 captures mostly what you expect but can grab some noise.

## Examples

### Example 1 - simple run without control peak and target peak between 1500 and 2000 nts

```sh
python arrayed_degradation/analyze_experiment.py --data_dir_path data --min_peak 1500 --max_peak 2000 --disable_control_peak
```

With the command above you will run the analysis; script will search for the input data within `data` folder (that is a relative path to the directory with the data; but you can also provide absolute path as well); target peak will be searched within **1500 - 2000** nucleotides window on electropherograms. Script **will not use control peak** for normalization.

### Example 2 - run with (1) control peak, (2) target peak background removal and (3) narrower peaks

```sh
python arrayed_degradation/analyze_experiment.py --data_dir_path data --min_peak 1500 --max_peak 2000 --remove_background --width_param 0.7
```

Difference here is that now script will search for control peak and normalize electropherograms traces with respect to them. Peak bounds are not provided so script will use default bounds (235, 310) for the control sequence. Additionally, script will remove backgroud from the target peak by drawing a line at base of the peak (keep in mind that it might not work well for certain cases like migrating peaks, so use this argument with caution and double check if areas are correctly determined in results.pdf). Finally, script will try to make peaks bounds more narrow because of decreasing `width_param` from default 0.85 into 0.7.

# Input data format

## Directory structure

Input data should be split across plates that were inserted into Fragment Analyzer. Each plate should have a distinct directory in which plate map and electropherogram files must be placed.

Here is an example of a proper directory structure:

```sh
my_experiment_dir
├── my_plate_1
│   ├── epg.csv
│   └── plate_map.csv
├── my_plate_2
│   ├── epg.csv
│   └── plate_map.csv
└── my_plate_3
    ├── epg.csv
    └── plate_map.csv
```

In this example we have a separate directory for our experiment called `my_experiment_dir` (you can of course name it differently). Within that directory we have 3 directories for 3 of our plates (`my_plate_1`, `my_plate_2`, etc..; again names can be different) that were used in the experiment. Within each plate's directory we **must have** 2 files named **exactly** `plate_map.csv`, `epg.csv`.

Please do not create any additional directories within experiment folder as the script will break. Additional files (like the ones with results) are okay.

## Input files

There are two obligatory input files in order to read data for a given plate:

- `plate_map.csv` - it should be a CSV file containg well positions of samples being analyzed.

    Here is a simplified example of what a proper plate_map file should look like:

    ```
    001_TP0,002_TP0,003_TP0,001_TP180,002_TP180,003_TP180
    004_TP0,005_TP0,006_TP0,004_TP180,005_TP180,006_TP180
    ```

    In this example in well `A3` we have a sample with label `003_TP0` which represents RNA having ID `003` and run at timepoint `0`h. At well `A6` there is a sample labeled as `003_TP180` which is the same RNA but run at timepoint `180`min. The sample labels should have this format **`{rna_id}_TP{timepoint}`**, where `rna_id` is your custom design identifier, and `timepoint` can be integer or float number that determines timepoints used (by default  **in minutes**; could be changed to hours or days with `--time_unit` argument). Do not add any additional suffixes as the script might not recognize RNA and timepoint properly. For each RNA and each technical replicate you must have sample with timepoint 0. Replicates of the same RNA - timepoint pair within one plate or across plates are acceptable. There are no strict layout requirements or assumptions for plating your samples, meaning that you can order them as you wish. Plate maps files do not need to be the same for all of your plates. However, it is advised to think of physically distinct plates as technical replicates (so exact copies of each other) and design your experiment to follow this pattern.

    **Missing wells** within plate map are acceptable, so if you don't want to include specifc RNA - timepoint pair on a particular plate in the analysis just remove it from the platemap and leave cell empty. Keep in mind though that you shoud not drop timepoint 0 as it might lead to improper results (we use timepoint 0 as a reference/baseline point for remaining timepoints). If you really want to remove timepoint 0 please remove as well all the others timepoint that are paired with it.

    **Do not** add any other labels within platemap besides actual samples. If your platemap contains ladder please remove it and leave the cell empty.
    
    **Do not** add any rows or column names in `plate_map.csv`; those will be added internally by the script during data processing. Adding that will cause incorrect well assignment by frameshifting and nonsense results as a consequence.

    Script will inform you in the logs how many unique RNA ids, replicates and timepoints have been detected. Please double check if the numbers match your expectations.

- `epg.csv` - this should be a Fragment Analyzer electropherogram file (result of capillary electrophoresis) with intensity traces and the corresponding nucleotide sizing.

    Here is a simplified toy example of what a proper plate_map file should look like:

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

# Toy example

In the repo you can find a working toy example to run the code. You can also crosscheck whether your experiment that you want to run has consistent format.
Data for running as well as expected results are within `data` directory. Electropherograms contain control peaks for normalization.

You can run it by executing:

```sh
python arrayed_degradation/analyze_experiment.py --min_peak 600 --max_peak 2500 --data_dir_path data 
```
