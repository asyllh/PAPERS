##Construct, Merge, Solve & Adapt (CMSA) for the SCPP

###User input arguments
* `.txt` SCPP instance file must be the first argument (`PAn_xx.txt` or `PRn_xx.txt`)
* `-e`: time limit for CMSA. Default = 3600s.
* `-t`: tau, minimum scoring distance value. Default = 70mm.
* `-W`: strip width/bin capacity. Default = 2500mm.
* `-p`: number of new solutions generated in each iteration. Default = 5.
* `-a`: maximum age. Default = 5.
* `-s`: seed value. Default = 1.
* `-xc`: time limit for Exact Cover. Default = 600s.

In all tests, seed value -s is set equal to instance number.

###Results
####File Name: `CMSAnWpa.txt` or `CMSRnWpa.txt`
* `CMSA`: Artificial instance type
* `CMSR`: Real instance type
* `n`: number of items / 100
* `W`: bin capacity/1000 (rounded down to nearest integer)
* `p`: number of new solutions generated in each iteration.
* `a`: maximum age.


####File Format
Each line in the output file contains the following results for each instance (in order):
* Instance number
* Theoretical minimum (t)
* Number of bins in the best solution found (#S)
* Solution quality - number of bins in best solution/theoretical minimum (q)
* Fitness value of the best solution found
* Number of CMSA iterations performed in time limit

####File Name: `CTAnWpa_inst.txt` or `CTRnWpa_inst.txt`
A second file is also created for each instance that records when a better solution is found.
* Line 1:
    * Instance number
    * Number of items
    * Instance type (1 = Artificial, 2 = Real)
    * Bin size
    * Number of new solutions generated in each iteration
    * Maximum age
    * Time limit for CMSA
    * Time limit for Exact Cover
* Line 2 contains column headings
* Lines 3 onwards contains the iteration number, the number of bins in the best-so-far solution, the fitness of the best-so-far solution, and the time at which the solution was found

####File Name: `XCAnWpa_inst.txt` or `XCRnWpa_inst.txt`
A third file is also created for each instance that records results after every run of the exact cover procedure.
* Line 1:
    * Instance number
    * Number of items
    * Instance type (1 = Artificial, 2 = Real)
    * Bin size
    * Number of new solutions generated in each iteration
    * Maximum age
    * Time limit for CMSA
    * Time limit for Exact Cover
* Line 2 contains column headings
* Lines 3 onwards contains the following in order:
    * Iteration number (of CMSA)
    * Number of bins in the set B
    * Number of bins in the best solution found by the exact cover procedure
    * Fitness value of the best solution found by the exact cover procedure
    * Number of solutions found by the exact cover procedure
    * Time taken to complete recursive search procedure (up to time limit) within the exact cover procedure.

####Spreadsheets
* CMSAResults.xlsx is made up of 12 sheets, each containing results from 50 instances for each of the different instance classes.



