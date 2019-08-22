##Evolutionary Algorithm (EA) for the SCPP

###User input arguments
* `.txt` SCPP instance file must be the first argument (`PAn_xx.txt` or `PRn_xx.txt`)
* `-e`: time limit for EA. Default = 600s.
* `-t`: tau, minimum scoring distance value. Default = 70mm.
* `-W`: strip width/bin capacity. Default = 2500mm.
* `-p`: number of solutions in population. Default = 25.
* `-x`: recombination operator. 1 = GGA, 2 = AGX, 3 = AGX'. Default = 1.
* `-s`: seed value. Default = 1.

In all tests, seed value -s is set equal to instance number.

###Results
####File Name: `EAnWx.txt` or `ERnWx.txt`
* `EA`: Artificial instance type
* `ER`: Real instance type
* `n`: number of items / 100
* `W`: bin capacity/1000 (rounded down to nearest integer)
* `x`: recombination operator. 1 = GGA, 2 = AGX, 3 = AGX'.

####File Format
Each line in the output file contains the following results for each instance (in order):
* Instance number
* Theoretical minimum (t)
* Number of bins in the best solution found (#S)
* Solution quality - number of bins in best solution/theoretical minimum (q)
* Fitness value of the best solution found
* Number of EA iterations performed in time limit
* Proportion of feasible solutions found from all instances of the SubSCP
* Time in which the best solution was found

####File Name: `ETAnWx_inst.txt` or `ETRnWx_inst.txt`
Another file is also created for each instance that records when a better solution is found during the EA.
* Line 1:
    * Instance number
    * Number of items
    * Instance type ( 1 = Artificial, 2 = Real)
    * Bin size
    * Recombination operator (1 = GGA, 2 = AGX, 3 = AGX')
    * Theoretical minimum
* Line 2 contains column headings
* Lines 3 onwards contains the number of bins in the best-so-far solution, the fitness of the best-so-far solution, and the time at which the solution was found

####Spreadsheets
* EAResults1 contains 6 spreadsheets each containing 6 sheets with results from EA experiments run for 10 minutes on 1000 instances.
* EAResults2 contains 6 spreadsheets each containing 6 sheets with results from EA experiments run for an extended period of 60 minutes on 50 instances.


