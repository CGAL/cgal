* Number-type: floating-point
* Description: Freed instance of original file
* Generated-by: master_mps_to_derivatives-create_free_instance
NAME QP example
ROWS
  N obj
  G c0
  L c1
  G c2
  L c3
  G c4
COLUMNS
  x0  obj  1.5
  x0  c0  2
  x0  c1  -1
  x0  c2  1
  x0  c3  1
  x1  obj  -2
  x1  c0  1
  x1  c1  2
  x1  c4  1
RHS
  rhs c0  2
  rhs c1  6
  rhs c3  20
BOUNDS
  MI  BND  x0
  MI  BND  x1
QMATRIX
  x0  x0  8
  x1  x0  2
  x0  x1  2
  x1  x1  10
ENDATA
