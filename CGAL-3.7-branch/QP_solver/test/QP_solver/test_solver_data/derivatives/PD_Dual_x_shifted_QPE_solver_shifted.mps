* Number-type: integer
* Description: Shifted instance of original file
* Generated-by: master_mps_to_derivatives-create_shifted_instance
NAME PD_Dual_x_shifted_QPE_solver
ROWS
  N obj
  L c0
  L c1
  L c2
  L c3
  L c4
  L c5
COLUMNS
  x0  obj  4
  x0  c0  -1
  x0  c1  -1
  x0  c2  4
  x1  obj  4
  x1  c0  -4
  x1  c1  1
  x1  c2  1
  x2  obj  -4
  x2  c3  3
  x2  c4  -3
  x3  obj  -4
  x3  c3  -2
  x3  c4  -1
  x3  c5  1
RHS
  rhs c0  -17
  rhs c1  3
  rhs c2  23
  rhs c3  8
  rhs c4  -32
  rhs c5  11
BOUNDS
  LO  BND  x0  1
  LO  BND  x1  2
  LO  BND  x2  3
  LO  BND  x3  4
QMATRIX
  x0  x0  2
  x1  x1  2
  x2  x0  -2
  x0  x2  -2
  x2  x2  2
  x3  x1  -2
  x1  x3  -2
  x3  x3  2
ENDATA
