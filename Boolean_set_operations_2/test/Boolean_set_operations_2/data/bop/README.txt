                          input files for test_bop (test boolean operations on two polygons)
                          ------------------------------------------------------------------

strcutre of file:
----------------

# first polygon
T1  # polygon's type : 0 is Polygon_2, 1 is Polygon_with_holes_2
P1  # the polygon

# second polygon
T2  # polygon's type
P2  # the polygon

# P1 + P2 (union)
BOOL  # Do P1 and P2 overlap? (can be 0 - no, 1 - yes)
NUION_RES #if BOOL is 1,NUION_RES is a Polygon_with_holes_2
           #else (if BOOL is 0), NUION_RES is empty (no result)


# P1 * P2 (intersection)
N # number of polygon with holes
POLYGON_WITH_HOLES # (Polygon_with_holes_2)
    |
    |
    |
POLYGON_WITH_HOLES # (Polygon_with_holes_2)

# P1 - P2 (difference of P1 and P2)
DIFF_RES # the same as in intersection

# P2 - P1 (difference of P2 and P1)
DIFF_RES # the same as in intersection

# P1 ^ P2 (symmetric difference)
SYMM_DIFF_RES # the same as in intersection

# ~P1  (complement of p1)
COMP_RES  # the same as in intersection

# ~P2 (complement of p2)
COMP_RES # the same as in intersection

# (oriented_side of P1 and P2)

example:
--------

# Polygon_2
0
3
0/1 0/1
10/1 0/1
0/1 10/1

# Polygon_with_holes_2
1
0
1
4
1/1 1/1
2/1 1/1
2/1 2/1
1/1 2/1

# union result
1  # P1 and P2 intersect
# the result is Polygon_with_holes_2 (representing the entire plane)
0
0

# intersection result
1
3
0/1 0/1
10/1 0/1
0/1 10/1
1
4
1/1 1/1
2/1 1/1
2/1 2/1
1/1 2/1

# difference between P1 and P2 result
1
4
1/1 1/1
2/1 1/1
2/1 2/1
1/1 2/1
0

# difference between P2 and P1 result
1
0
1
3
0/1 0/1
10/1 0/1
0/1 10/1

# symmetric difference result
2
0
1
3
0/1 0/1
10/1 0/1
0/1 10/1
4
1/1 1/1
2/1 1/1
2/1 2/1
1/1 2/1
0

# complement of P1
1
0
1
3
0/1 0/1
10/1 0/1
0/1 10/1

# complement of P2
1
4
1/1 1/1
2/1 1/1
2/1 2/1
1/1 2/1
0

# oriented side of P1 and P2
1
