## README.md

This is the peak finder for NV centers at the Harvard Yao group. 

Here's the general troubleshooting flowchart for the algorithm:
![Untitled](https://github.com/AnmayG/peak_finder/assets/84810366/541b21ce-59de-41ec-af73-d486acfdcdb7)

In order to quickly find the right settings for the peak finder, I suggest using the ```global_peaks_gui.mlapp``` file to quickly and easily find settings according to the above flowchart.

In general, the peak finder is split across 5 files:
- ```find_peaks_at_point.m``` contains the logic for the first pass of peak finding
- ```partner_peak.m``` contains the logic for peak pairing based on the split, including phantom peak finding
- ```peaks_to_seed.m``` transforms the given peak pairings and points into physically interpretable variables such as D_111 and E_111
- ```find_baseline.m``` finds the baseline of any given pixel by binning all of the data and taking the bottom n%.
- ```peak_find_function.m``` wraps all of these into one function for easy use.

Testing wise, I suggest using the ```global_peaks_gui.mlapp``` to debug anything and everything related to the peak finder itself, 
but if you'd prefer to run it on the command line the `test.m` file and the `full_peak_find.m` file are both decent debugging tools. 

The `test.m` file returns the signal from a set of random points in a grid, with the size of each cell in the grid decided by the `chunk_size` parameter. 
The `full_peak_find.m` file conducts peak finding on the entire image and returns a set of plots depicting the entire shape of the data.

Have fun!
