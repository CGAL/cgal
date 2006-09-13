* Number-type: integer
* Description: Shifted instance of original file
* Generated-by: master_mps_to_derivatives-create_shifted_instance
NAME LP_KM3_Dual_QPE_solver
ROWS
  N obj
  G c0
  G c1
  G c2
COLUMNS
  x0  obj  1
  x0  c0  1
  x1  obj  100
  x1  c0  20
  x1  c1  1
  x2  obj  10000
  x2  c0  200
  x2  c1  20
  x2  c2  1
RHS
  rhs c0  741
  rhs c1  72
  rhs c2  4
BOUNDS
  LO  BND  x0  1
  LO  BND  x1  2
  LO  BND  x2  3
ENDATA
