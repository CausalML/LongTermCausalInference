This code is used to generate Table 1 in the main text, 
as well as Tables 2 and 5-9 in the Supplement.

To generate Table 1, first place the data file `quartly.mat`
(available upon request) into this folder and run `R --no-save --args < pipeline_data.R` 
to create the semi-synthetic data. Then, run `./run_all.sh` 
to execute each method. Finally, run `R --no-save --args < test_final.R` 
to obtain the table.

To generate Table 2 in the Supplementary Material, 
follow the same steps as above, but in the last step, 
run `R --no-save --args < test_raw.R`.

To generate Tables 5-8 in the Supplementary Material, 
modify the code in `pipeline_data.R` according to the 
instructions in the comments (lines 5-12) in that R 
script to generate the semi-synthetic data, then
follow the same steps as for Table 1.

To generate Table 9, follow the same procedure as for 
Table 1, but change the code in the R scripts `main.R`, 
`main_ridge_n1.R` and `main_ridge_n3.R` according to 
the comments in those scripts. See lines 91-94 in `main.R`, 
lines 102-105 in `main_ridge_n1.R`, and lines 96-99 in 
`main_ridge_n3.R` for those comments.
