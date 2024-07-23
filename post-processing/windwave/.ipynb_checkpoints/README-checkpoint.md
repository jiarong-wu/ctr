This is the code for wind wave post-processing.

Install by Install by running

pip install -e . 


Post-processing procedures:
1. dump files (pretty frequently) if storage can afford
2. run ./field_outputs TIME NSLICE to output slices and span-wise averaged results
3. Read in slices and convert to compressed netcdf files
4. (Delete the dump files and slice files?)
5. Analyse!

Dependence:
tqdm