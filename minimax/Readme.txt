This code is used to generate Table 2 in the main text 
as well as Table 10 in the supplement.

To generate the result for each cell in the table, first edit 
the R script `data.R` to specify the dimension of X, S1/S3, U, 
as well as the degree of nonlinearity q. Next edit the python 
script `surrogate.py' and the R script `ridge.R` to specify the 
dimensions of X, S1/S3, and the number of neurons in the first 
hidden layer of the neural network. After making these edits, 
run `R --no-save --args < data.R` first, then run `./run_all.sh`, 
finally execute `R --no-save --args < evaluate.R' to generate the 
results for that cell.
