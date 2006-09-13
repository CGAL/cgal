* Number-type: integer
* Description: Freed instance of original file
* Generated-by: master_mps_to_derivatives-create_free_instance
NAME LP_KM3_QPE_solver
ROWS
  N obj
  L c0
  L c1
  L c2
  G c3
  G c4
  G c5
COLUMNS
  x0  obj  -100
  x0  c0  1
  x0  c1  20
  x0  c2  200
  x0  c3  1
  x1  obj  -10
  x1  c1  1
  x1  c2  20
  x1  c4  1
  x2  obj  -1
  x2  c2  1
  x2  c5  1
RHS
  rhs c0  1
  rhs c1  100
  rhs c2  10000
BOUNDS
  MI  BND  x0
  MI  BND  x1
  MI  BND  x2
ENDATA
