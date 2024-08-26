This is code is used to generate Tables 1, 4-5 in the 
supplement.

To generate each table, first put the data file "quartly.mat" 
into this folder and change the code to specify the choice of 
Surrogates used to calculate the conditional distribution matrix
P(S1 | S2, A, U). Then run "R --no-save --args < get_rank.R"