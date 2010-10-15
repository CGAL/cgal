* Number-type: integer
* Description: Shifted instance of original file
* Generated-by: master_mps_to_derivatives-create_shifted_instance
NAME QP_Bounded_QPE_solver
ROWS
  N obj
  G c0
COLUMNS
  x0  obj  -2
  x0  c0  1
  x1  obj  -1
  x1  c0  -1
RHS
  rhs c0  1
BOUNDS
  LO  BND  x0  1
  LO  BND  x1  2
QMATRIX
  x0  x0  2
ENDATA
