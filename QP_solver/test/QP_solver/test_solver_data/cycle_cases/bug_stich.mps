* Description: Case that was sent in by Sebastian Stich. It used to exibit
* a segfault in expel_artificial_variables (special artificial is linked to a
* constraint that is not active).
* Derivatives: none
NAME MY_MPS
ROWS
  N obj
  G c0
  G c1
  E c2
  G c3
  G c4
  G c5
  G c6
  G c7
  G c8
  E c9
  G c10
  G c11
  G c12
  E c13
  G c14
  G c15
  E c16
  G c17
  G c18
  G c19
  E c20
COLUMNS
  x0  c0  1
  x0  c1  1
  x0  c2  1
  x0  c3  1
  x0  c4  1
  x1  c5  1
  x1  c6  1
  x1  c7  1
  x1  c8  1
  x1  c9  1
  x1  c10  1
  x1  c11  1
  x2  c0  1
  x2  c5  1
  x2  c12  1
  x2  c13  1
  x3  c1  1
  x3  c6  1
  x3  c12  1
  x3  c14  1
  x3  c15  1
  x3  c16  1
  x4  c7  1
  x4  c13  1
  x4  c17  1
  x4  c18  1
  x4  c19  1
  x5  c2  1
  x6  c3  1
  x6  c8  1
  x6  c14  1
  x6  c17  1
  x6  c20  1
  x7  c9  1
  x8  c4  1
  x8  c10  1
  x8  c15  1
  x8  c18  1
  x8  c20  1
  x9  c11  1
  x9  c16  1
  x9  c19  1
RHS
  rhs c0  1
  rhs c1  1
  rhs c2  1
  rhs c3  1
  rhs c4  1
  rhs c5  1
  rhs c6  1
  rhs c7  1
  rhs c8  1
  rhs c9  1
  rhs c10  1
  rhs c11  1
  rhs c12  1
  rhs c13  1
  rhs c14  1
  rhs c15  1
  rhs c16  1
  rhs c17  1
  rhs c18  1
  rhs c19  1
  rhs c20  1
BOUNDS
  UP  BND  x0  1
  UP  BND  x1  1
  UP  BND  x2  1
  UP  BND  x3  1
  UP  BND  x4  1
  UP  BND  x5  1
  UP  BND  x6  1
  UP  BND  x7  1
  UP  BND  x8  1
  UP  BND  x9  1
ENDATA
