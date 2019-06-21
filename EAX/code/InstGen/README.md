##Problem Instance Generator for the SCPP

###User input arguments
* `-i`: instance number
* `-t`: tau, minimum scoring distance value. Default: 70mm.
* `-n`: number of items. Default: 100.
* `-a`: minimum score width. Default: 1mm.
* `-b`: maximum score width. Default: 70mm.
* `-m`: minimum item width. Default: 150mm.
* `-M`: maximum item width. Default: 1000mm.
* `-c`: instance type, 1 = artificial, 2 = real. Default = 1.
* `-s`: seed value. Default = 1.

Example .bat file provided, with seed set equal to instance number.

###Output File
####Filename: PAn_i.txt or PRn_i.txt
* PA: Artificial instance type
* PR: Real instance type
* n: number of items / 100
* i: instance number

####File Format
* Lines 9 to 13 are only produced for real instance types.
1. Instance number
2. Number of items n
3. Instance type: 1 = artificial, 2 = real
4. AllScores - score widths in increasing order, size = 2n
5. Partners - score width indices and corresponding partner index, size = 2n
6. Item widths - score width indices and corresponding item width, size = 2n
7. Item number - score width indices and corresponding item number, size = 2n
8. Total width of all items
9. Number of item types
10. Number of items of each type
11. Value of smaller score width for each type
12. Value of larger score width for each type
13. Item width for each type
14. Type number - score width indices and corresponding type number


