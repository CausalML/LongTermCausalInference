This is code is used to generate Table 1 in the main text 
as well as Tables 2 and 5-9 in the supplement.

To generate Table 1, first put the data file "quartly.mat" 
into this folder and run "R --no-save --args < pipeline_data.R"
to create the processed data, then run "./run_all.sh" to 
execute each method, finally run "R --no-save --args < test_final.R"
to obtain all data for the table.

To generate Table 2 in the Supplementary Material, use the same 
steps for Table 2, except that run "R --no-save --args < test_raw.R" 
in the last step.

To generate Tables 5-8 in the Supplementary Material, change the 
code "pipeline_data.R" according to the instructions that are 
appeared as the comment in the source file "pipeline_data.R" to 
regenerate the processed data, then follow the same steps as in 
generating Table 1 to generate Tables 5-8.

To generate Table 9, you need to follow the same procedure as 
generating Table 1, except that you need to change the code in 
R scripts "main.R", "main_ridge_n1.R", "main_ridge_n3.R" according 
to the comments in those scripts.
