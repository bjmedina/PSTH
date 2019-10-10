# Plots for _Neuropixels_ Data
---

These files assume the following hierarchy: 

```misc
~/
  data/ (code and data go here) [in code, put path in DIRECTORY]
    plots/ [in code, put path in SUMMARY_PLOTS_DIRECTORY]
      cells/ (PSTH plots for every cell in every probe) [put path in CELL_PLOTS_DIRECTORY]
        probeA/ [ you just have to make these directories, with those exact names ]
        probeB/
        ...
        probeF/
      percentile/ (Plots for percentiles of every mouse)
      probes/ (PSTH plots averaged across all cells and trials) [put path in PROBE_PLOTS_DIRECTORY]
      variations/ (directories of class names) [in code, put path in VAR_DIREC]
    trial_data/ (data calculated for every trial)
```


## Libraries Needed
---
``` library >= version ``` || ``` pip install library ```

```rpy2 >= 3.0.4``` || ```pip install rpy2``` 

```numpy >= 1.14.2``` || ```pip install numpy```

```pickle >= 4.0``` || ```pip install pickle```

```h5py >= 2.8.0``` || ```pip install h5py```

```matplotlib >= 3.0.2``` || ```pip install matplotlib```


## Steps for Running
---

1. Make sure you have all the directories and folders set up. Follow the ```[...put path in (VARIABLE_NAME_HERE)]``` comments to help you figure out where to put the paths in the code (Code could be written to automatically create the directories...).

2. Run ```nwb_plots.py```. If everything goes well, the output before the program ends should be ```run again```.

3. Run ```nwb_plots.py``` again. Repeat 2 and 3 for all mice (you can change the mouse you're working on by looking for the variable ```MOUSE_ID``` in the code). You should now see plots in ```SUMMARY_PLOTS_DIRECTORY``` and ```CELL_PLOTS_DIRECTORY```. These two steps should take the longest.

4. Run ```nwb_plots_percentile.py``` for each mouse. These plots will be saved in ```percentile```.

5. Run ```nwb_dropout.py``` for each mouse. These plots will be saved in ```probes```.

__NOTE__: ```VAR_DIREC``` is no longer needed. Neither is ```nwb_trials.py```
