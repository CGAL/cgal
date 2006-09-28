* Number-type: integer
* Description: Unbounded problem with asymmetric D matrix
* Generated-by: 
NAME QPE_solver_example_bug
ROWS
 N obj
 G c1
 G c2
COLUMNS
  x1  obj  -64
  x2  obj  33
  x3  obj  -8
  x1  c1  -4
  x2  c1  2
  x3  c1  0
  x1  c2  0
  x2  c2  1
  x3  c2  1
RHS
  rhs  c1  -8
  rhs  c2  7
BOUNDS
  LO  BND  x1  1
  LO  BND  x2  2
  LO  BND  x3  3
DMATRIX
  x1  x1  64
  x1  x2  -16
  x2  x1  -16
  x2  x2  4
ENDATA
