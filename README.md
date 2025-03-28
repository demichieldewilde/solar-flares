**Integrated full disc Sun-as-a-star flare observations synthesised
from a small field of view**

This code is used for the data analysis of the paper "Integrated full disc Sun-as-a-star flare observations synthesised
from a small field of view" By M. De Wilde, A. G.M. Pietrow , M. k. Druett , A. P. Yabar, 
J. Koza, O. Andriienko, A. R. Brunvoll, J. de la Cruz Rodríguez, J. T. Faber, R. Joshi, D. Kuridze,
D. Nóbrega-Siverio, L. H. M. Rouppe van der Voort, J. Rybák, E. Scullion, A. M. Silva,
Z. Vashalomidze, A. Vicente Arévalo, R. Yadav, T. V. Zaqarashvili, J. Zbinden, and E. S. Øyre

Videos of the observations dicussed in the paper can be found under data/animations. 

-----------------------------------------------------------------------------------------------

The raw observational data is not included for memory and copyright reasons. 
Some of these data is free available for example for the flares 

    -2011-08-06
    -2014-09-06
    -2015-06-24

via (https://dubshen.astro.su.se/sst_archive/observations/281)

For those flare the full analysis thus can be repeated yourself. 

Schematically what we have done is the following:

    1. Gauge the quiet sun to a NESSI profile for each line in each flare. These are all notebooks, that can be found under a date
    2. Create sun-as-a-star datasets from this gauge (most .npy, .npz files in the same folder) for further investigation
    3. Fitted a voight profile to the difference in Field Of View (FOV) en quiet sun, This difference are the contrastprofiles. We also calculated the residues. (in folder voight_fitting). 
    4. Plotting the full analysis of the flare (in folder full_analysis)


When you want to take a look I advise the 2014-09-06 and the 2015-06-24, these are the best structured, 
and therefore the most easy to understand. 


