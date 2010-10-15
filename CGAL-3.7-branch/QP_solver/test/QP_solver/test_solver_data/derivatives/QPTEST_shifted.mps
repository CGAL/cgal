* Number-type: floating-point
* Description: Shifted instance of original file
* Generated-by: master_mps_to_derivatives-create_shifted_instance
NAME QP example
ROWS
  N obj
  G c0
  L c1
COLUMNS
  x0  obj  -10.5
  x0  c0  2
  x0  c1  -1
  x1  obj  -24
  x1  c0  1
  x1  c1  2
RHS
  rhs c0  6
  rhs c1  9
BOUNDS
  LO  BND  x0  1
  UP  BND  x0  21
  LO  BND  x1  2
QMATRIX
  x0  x0  8
  x1  x0  2
  x0  x1  2
  x1  x1  10
ENDATA
