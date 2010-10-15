* Number-type: integer
* Description: Freed instance of original file
* Generated-by: master_mps_to_derivatives-create_free_instance
NAME QP_Unbounded_QPE_solver
ROWS
  N obj
  L c0
  G c1
  G c2
COLUMNS
  x0  c0  1
  x0  c1  1
  x1  obj  -1
  x1  c0  -1
  x1  c2  1
RHS
BOUNDS
  MI  BND  x0
  MI  BND  x1
QMATRIX
  x0  x0  2
ENDATA
