* Number-type: integer
* Description: Freed instance of original file
* Generated-by: master_mps_to_derivatives-create_free_instance
NAME QP_Bounded_QPE_solver
ROWS
  N obj
  G c0
  G c1
  G c2
COLUMNS
  x0  obj  0
  x0  c0  1
  x0  c1  1
  x0  c2  0
  x1  obj  -1
  x1  c0  -1
  x1  c1  0
  x1  c2  1
RHS
  rhs c0  2
  rhs c1  0
  rhs c2  0
BOUNDS
  MI  BND  x0
  MI  BND  x1
QMATRIX
  x0  x0  2
ENDATA
