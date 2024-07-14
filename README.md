# Bat - A tool for post-processing experimental impedance tube data
### Renan Liupekevicius Carnielli [r.liupekevicius.carnielli@tue.nl](mailto::r.liupekevicius.carnielli@tue.nl)
Script adapted from L. Schijff (https://research.tue.nl/en/studentTheses/a-numerical-experimental-analysis-of-metafoams-and-their-design), A. Karunarathne, and N.A. van de Straat.

## What is Bat?
Bat is a tool to compute the following acoustic indicators for a sample 
- Transmission Loss [dB]
- Trasmission coefficient [-];
- Reflection coefficient [-];
- Absorption coefficient [-].

The Two-Load Method is implemented following [1] and [2].

## How to use
 - 1. Add sample data to the folder \Bat\data (5-mic measures/load). Example: create folder \Bat\data\sample_example
 - 2. Run `startup.m` to start the analysis;
 - 3. Indicate in `two_load_exp_postprocessing.m` the sample path & set parameters highlighted by the sandwitch of '%%%%'
 - 4. Run `two_load_exp_postprocessing.m` & access the plots in the corresponding sample folder ( e.g. \Bat\data\sample_example)
## Note 
 - 1. Do not rename the folders within \Bat. You are allowed to create folders within \Bat\data.
 
## References 

[1] Bolton, J. Stuart, Taewook Yoo, and Oliviero Olivieri. "Measurement of normal incidence transmission loss and other acoustical properties of materials placed in a standing wave tube." * *Brüel & Kjær Technical Review* * 1 (2007): 1-44.

[2] Pispola, Giulio, Kirill V. Horoshenkov, and Francesco Asdrubali. "Transmission loss measurement of consolidated granular media (L)." * *The Journal of the Acoustical Society of America* * 117.5 (2005): 2716-2719.
