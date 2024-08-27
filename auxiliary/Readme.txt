This code is used to generate Tables 1, 4-5 and Figure 3 in the 
supplement.

To generate the tables, first place the data file `quartly.mat` into 
this folder. Then, modify the R script `distribution_matrix_svd.R` to 
specify the choice of surrogates used to calculate the conditional 
distribution matrix P(S1 | S2, A, U), following the instructions in 
lines 5-14 of that R script. Finally, run the command `R --no-save --args < distribution_matrix_svd.R`

To generate the figure, run `R --no-save --args < bridge_existence_figure.R`
