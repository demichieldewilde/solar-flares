**Integrated full disc Sun-as-a-star flare observations synthesised
from a small field of view**

This repo accompanies the paper "Integrated full disc Sun-as-a-star flare observations synthesised
from a small field of view" By M. De Wilde, A. G.M. Pietrow , M. k. Druett , A. P. Yabar, 
J. Koza, O. Andriienko, A. R. Brunvoll, J. de la Cruz Rodríguez, J. T. Faber, R. Joshi, D. Kuridze,
D. Nóbrega-Siverio, L. H. M. Rouppe van der Voort, J. Rybák, E. Scullion, A. M. Silva,
Z. Vashalomidze, A. Vicente Arévalo, R. Yadav, T. V. Zaqarashvili, J. Zbinden, and E. S. Øyre

The code for all data analysis can be found in the data folder which is structured to have an 
individual folder extracting data in a uniform format for each event, and then a folder per 
stage of the data analysis (voight fitting, analysis of each event in folder full analysis, 
and scale law). 

The Figures used in the paper are assembled in Figures.

The videos of the observations dicussed in the paper can be found in animations. 

-----------------------------------------------------------------------------------------------

The raw observational data is not included for memory and copyright reasons. 
Some of these data is freely available for example for the flares 

    -2011-08-06
    -2014-09-06
    -2015-06-24

via (https://dubshen.astro.su.se/sst_archive/observations/281)

For those flare the full analysis can be repeated yourself. 

Schematically what we have done is the following:

    1. Gauge the quiet sun to a NESSI profile for each line in each flare. These are all notebooks, that can be found under a date
    2. Create sun-as-a-star datasets from this gauge (most .npy, .npz files in the same folder) for further investigation
    3. Fitted a voight profile to the difference in Field Of View (FOV) en quiet sun, This difference are the contrastprofiles. We also calculated the residues. (in folder voight_fitting). 
    4. Plotting the full analysis of the flare (in folder full_analysis)


When you want to take a look I advise the 2014-09-06 and the 2015-06-24, these are the best structured, 
and therefore the most easy to understand. 


