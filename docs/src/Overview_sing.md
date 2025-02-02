@doc raw"""
proc printMat(matrix M)

**USAGE**   :  printMat(M); M matrix@*
**ASSUME**  :  M is a matrix.
**THEORY**  :  This is the print function used by Singular to print a matrix.
**KEYWORDS**: matrix
```singular
ring R=0,(x),lp;
matrix M[2][3]=1,243,3,4,522222,6;
printMat(M);
````
"""
