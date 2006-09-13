* Number-type: integer
* Description: Freed instance of original file
* Generated-by: master_mps_to_derivatives-create_free_instance
NAME LP_KM4_QPE_solver
ROWS
  N obj
  L c0
  L c1
  L c2
  L c3
  G c4
  G c5
  G c6
  G c7
COLUMNS
  x0  obj  -1000
  x0  c0  1
  x0  c1  20
  x0  c2  200
  x0  c3  2000
  x0  c4  1
  x1  obj  -100
  x1  c1  1
  x1  c2  20
  x1  c3  200
  x1  c5  1
  x2  obj  -10
  x2  c2  1
  x2  c3  20
  x2  c6  1
  x3  obj  -1
  x3  c3  1
  x3  c7  1
RHS
  rhs c0  1
  rhs c1  100
  rhs c2  10000
  rhs c3  1000000
BOUNDS
  MI  BND  x0
  MI  BND  x1
  MI  BND  x2
  MI  BND  x3
ENDATA
