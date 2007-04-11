* Number-type: integer
* Description: Shifted instance of original file
* Generated-by: master_mps_to_derivatives-create_shifted_instance
NAME QPE_solver_example_bug
ROWS
  N obj
  G c0
  G c1
COLUMNS
  x0  obj  -128
  x0  c0  -4
  x1  obj  49
  x1  c0  2
  x1  c1  1
  x2  obj  -8
  x2  c1  1
RHS
  rhs c0  -8
  rhs c1  12
BOUNDS
  LO  BND  x0  2
  LO  BND  x1  4
  LO  BND  x2  6
QMATRIX
  x0  x0  128
  x1  x0  -32
  x0  x1  -32
  x1  x1  8
ENDATA
