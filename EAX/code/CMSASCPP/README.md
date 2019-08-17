##Construct, Merge, Solve & Adapt (CMSA) for the SCPP

###User input arguments
* `.txt` SCPP instance file must be the first argument (`PAn_xx.txt` or `PRn_xx.txt`)
* `-e`: time limit for CMSA. Default = 3600s.
* `-t`: tau, minimum scoring distance value. Default: 70mm.
* `-W`: strip width/bin capacity. Default:2500mm.
* `-p`: number of new solutions generated in each iteration. Default = 5.
* `-a`: maximum age. Default = 5.
* `-s`: seed value. Default = 1.
* `-xc`: time limit for Exact Cover. Default = 600s.

In all tests, seed value -s is set equal to instance number.

###Results
####File Name: `CAnWpa.txt` or `CRnWpa.txt`
* `CA`: Artificial instance type
* `CR`: Real instance type
* `n`: number of items / 100
* `W`: bin capacity/1000 (rounded down to nearest integer)
* `p`: number of new solutions generated in each iteration.
* `a`: maximum age.

####File Format
Line 1: File name
Line 2: `n`, instance type (1 = artificial, 2 = real), `-W`, `-p`, `-a`, `-e`, `-xc`
Line 3: Instance number, theoretical minimum, number of bins in best solution found, solution quality, fitness value, number of iterations of CMSA
Line 4 onwards: results for each of the 50 instances.


