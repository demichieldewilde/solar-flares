**Integrated full disc Sun-as-a-star flare observations synthesised
from a small field of view**

This code is used for the data analysis of the paper "Integrated full disc Sun-as-a-star flare observations synthesised
from a small field of view" By Michiel De Wilde, Alexander G.M. Pietrow , Malcolm Keith Druett , and Adur Pastor Yabar


-----------------------------------------------------------------------------------------------

Only the observational data is not included for memory and copyright reasons. 
Some of these data is free available for example for the flares 
    -2011-08-06
    -2014-09-06
    -2015-06-24

For those flare the full analysis thus can be repeated yourself. 

Schematically what we have done is the following:

    1. Gauge the quiet sun to a NESSI profile for each line in each flare. These are all notebooks, that can be found under a date
    2. Create sun-as-a-star datasets from this gauge (most .npy, .npz files in the same folder) for further investigation
    3. Fitted a voight profile to the difference in Field Of View (FOV) en quiet sun, This difference are the contrastprofiles. We also calculated the residues. (in folder voight_fitting). 
    4. Plotting the full analysis of the flare (in folder full_analysis)


When you want to take a look I advise the 2014-09-06 and the 2015-06-24, these are the best structured, 
and therefore the most easy to understand. 
