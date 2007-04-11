* Number-type: integer
* Description: Shifted instance of original file
* Generated-by: master_mps_to_derivatives-create_shifted_instance
NAME PD_Primal_QPE_solver
ROWS
  N obj
  E c0
  E c1
COLUMNS
  x0  obj  616
  x0  c0  1
  x1  obj  1108
  x1  c0  1
  x2  obj  1540
  x2  c0  1
  x3  obj  -1712
  x3  c1  1
  x4  obj  -2476
  x4  c1  1
  x5  obj  -1996
  x5  c1  1
RHS
  rhs c0  7
  rhs c1  16
BOUNDS
  LO  BND  x0  1
  LO  BND  x1  2
  LO  BND  x2  3
  LO  BND  x3  4
  LO  BND  x4  5
  LO  BND  x5  6
QMATRIX
  x0  x0  16
  x1  x0  28
  x0  x1  28
  x1  x1  74
  x2  x0  40
  x0  x2  40
  x2  x1  70
  x1  x2  70
  x2  x2  100
  x3  x0  -44
  x0  x3  -44
  x3  x1  -92
  x1  x3  -92
  x3  x2  -110
  x2  x3  -110
  x3  x3  130
  x4  x0  -64
  x0  x4  -64
  x4  x1  -122
  x1  x4  -122
  x4  x2  -160
  x2  x4  -160
  x4  x3  182
  x3  x4  182
  x4  x4  260
  x5  x0  -52
  x0  x5  -52
  x5  x1  -86
  x1  x5  -86
  x5  x2  -130
  x2  x5  -130
  x5  x3  140
  x3  x5  140
  x5  x4  206
  x4  x5  206
  x5  x5  170
ENDATA
