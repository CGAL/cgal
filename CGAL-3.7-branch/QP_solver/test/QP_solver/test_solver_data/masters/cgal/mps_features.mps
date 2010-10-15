* Number-type: rational
* Description: 
* min x^2 + y^2 + 1/2
*    1/3   <= x + y  = 1     (lower bound from RANGES) 
*    1/4   <= x     <= 1     (lower bound from RANGES)
*    0     <= y     <= 3/4   (upper bound from RANGES)
*
*    x, y  >= 0
NAME 	file with:  complex name,  constant term,  ranges  
ROWS
 N  obj
 E  xplusy
 L  xonly
 G  yonly
COLUMNS
    x	xplusy	1/1	xonly 1/1	
    y	xplusy 	1/1	yonly 1/1
RHS
    rhs	obj    -1/2
    rhs xplusy	1/1
    rhs xonly   1/1
    rhs yonly   0/1
RANGES
    range xplusy -2/3
    range xonly	  3/4
    range yonly   3/4  	
DMATRIX
    x  x	1/1
    y  y 	1/1    	
ENDATA
