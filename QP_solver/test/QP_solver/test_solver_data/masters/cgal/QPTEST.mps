* Description: from the benchmarks at http://www.doc.ic.ac.uk/~im/
NAME          QP example
ROWS
 N  obj
 G  r1
 L  r2
COLUMNS
    c1        r1                 2.0   r2                -1.0
    c1        obj                1.5
    c2        r1                 1.0   r2                 2.0
    c2        obj               -2.0
RHS
    rhs1      r1                 2.0   r2                 6.0
BOUNDS
 UP bnd1      c1                20.0
QUADOBJ
    c1        c1                 8.0
    c1        c2                 2.0
    c2        c2                10.0
ENDATA
