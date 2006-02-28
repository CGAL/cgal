* Number-type: integer
* Description: Freed instance of original file
* Generated-by: master_mps_to_derivatives-create_free_instance
NAME PD_Dual_Eq_QPE_solver
ROWS
  N obj
  E c0
  L c1
  L c2
  L c3
  L c4
  L c5
  G c6
  G c7
  G c8
  G c9
COLUMNS
  x0  obj  0
  x0  c0  -1
  x0  c1  -1
  x0  c2  4
  x0  c3  0
  x0  c4  0
  x0  c5  0
  x0  c6  1
  x0  c7  0
  x0  c8  0
  x0  c9  0
  x1  obj  0
  x1  c0  -4
  x1  c1  1
  x1  c2  1
  x1  c3  0
  x1  c4  0
  x1  c5  0
  x1  c6  0
  x1  c7  1
  x1  c8  0
  x1  c9  0
  x2  obj  0
  x2  c0  0
  x2  c1  0
  x2  c2  0
  x2  c3  3
  x2  c4  -3
  x2  c5  0
  x2  c6  0
  x2  c7  0
  x2  c8  1
  x2  c9  0
  x3  obj  0
  x3  c0  0
  x3  c1  0
  x3  c2  0
  x3  c3  -2
  x3  c4  -1
  x3  c5  1
  x3  c6  0
  x3  c7  0
  x3  c8  0
  x3  c9  1
RHS
  rhs c0  -10
  rhs c1  0
  rhs c2  25
  rhs c3  13
  rhs c4  -25
  rhs c5  7
  rhs c6  0
  rhs c7  0
  rhs c8  0
  rhs c9  0
BOUNDS
  MI  BND  x0
  MI  BND  x1
  MI  BND  x2
  MI  BND  x3
QMATRIX
  x0  x0  2
  x0  x2  -2
  x1  x1  2
  x1  x3  -2
  x2  x0  -2
  x2  x2  2
  x3  x1  -2
  x3  x3  2
ENDATA
