* Description: Case that was sent in by William Taysom. It used to exibit
* a segfault in expel_artificial_variables (special artificial is linked to a
* constraint that is not active).
NAME MY_MPS
ROWS
  N obj
  E c0
  G c1
  G c2
  G c3
COLUMNS
  x0  obj  1
  x0  c0  1
  x0  c1  1
  x0  c2  1
  x0  c3  1
  x1  c0  1
  x1  c1  1
  x1  c3  1
  x2  c0  1
  x2  c1  1
  x2  c2  1
RHS
  rhs c0  4
  rhs c1  4
  rhs c2  3
  rhs c3  3
BOUNDS
  UP  BND  x0  2
ENDATA
