This is code is used to generate Table 2 in the main text 
as well as Table 10 in the supplement.

To generate data in each cell, first edit the R script "data.R"
to specify the dimension of X, S1/S3, U as well as the degree
of nonlinearity q, then edit the python script "surrogate.py"
and "ridge.R" to specify the dimension of X, S1/S3 as well as 
the number of neurons in the first hidden layer of the neural 
network. After that, first run "./run_all.sh", then run 
"R --no-save --args < evaluate.R"
to generate the result to fill in the table.
