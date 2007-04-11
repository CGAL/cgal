* Number-type: integer
* Description: Freed instance of original file
* Generated-by: master_mps_to_derivatives-create_free_instance
NAME PD_Primal_QPE_solver
ROWS
  N obj
  E c0
  E c1
  G c2
  G c3
  G c4
  G c5
  G c6
  G c7
COLUMNS
  x0  c0  1
  x0  c2  1
  x1  c0  1
  x1  c3  1
  x2  c0  1
  x2  c4  1
  x3  c1  1
  x3  c5  1
  x4  c1  1
  x4  c6  1
  x5  c1  1
  x5  c7  1
RHS
  rhs c0  1
  rhs c1  1
BOUNDS
  MI  BND  x0
  MI  BND  x1
  MI  BND  x2
  MI  BND  x3
  MI  BND  x4
  MI  BND  x5
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
