* Number-type: integer
* Description: Freed instance of original file
* Generated-by: master_mps_to_derivatives-create_free_instance
NAME QP_Degenerate_rank_zero_QPE_solver
ROWS
  N obj
  E c0
  E c1
  E c2
  G c3
  G c4
  G c5
COLUMNS
  x0  obj  -1
  x0  c3  1
  x1  obj  -3
  x1  c4  1
  x2  obj  4
  x2  c5  1
RHS
BOUNDS
  MI  BND  x0
  MI  BND  x1
  MI  BND  x2
QMATRIX
  x0  x0  2
  x1  x1  2
  x2  x2  2
ENDATA
