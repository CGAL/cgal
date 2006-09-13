* Number-type: integer
* Description: Freed instance of original file
* Generated-by: master_mps_to_derivatives-create_free_instance
NAME LP_KM3_Dual_QPE_solver
ROWS
  N obj
  G c0
  G c1
  G c2
  G c3
  G c4
  G c5
COLUMNS
  x0  obj  1
  x0  c0  1
  x0  c3  1
  x1  obj  100
  x1  c0  20
  x1  c1  1
  x1  c4  1
  x2  obj  10000
  x2  c0  200
  x2  c1  20
  x2  c2  1
  x2  c5  1
RHS
  rhs c0  100
  rhs c1  10
  rhs c2  1
BOUNDS
  MI  BND  x0
  MI  BND  x1
  MI  BND  x2
ENDATA
