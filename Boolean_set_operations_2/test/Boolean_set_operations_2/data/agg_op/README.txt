                          input files for test_agg_op (test aggregate operations on many polygons)
                          ------------------------------------------------------------------------

structure of file:
-----------------

# polygons (range of Polygon_2)
N # number of polygons
POLYGON
   |
   |
POLYGON

#polygons with holes (range of Polygon_with_holes_2)
N # number of polygons with holes
POLYGON_WITH_HOLES
   |
   |
POLYGON_WITH_HOLES



# aggregate union on Polygon_2 range
N # number of polygon with holes
POLYGON_WITH_HOLES # (Polygon_with_holes_2)
    |
    |
    |
POLYGON_WITH_HOLES # (Polygon_with_holes_2)

# aggregate union on Polygon_with_holes_2 range
RESULT # the same as above

# aggregate union on Polygon_2 and  Polygon_with_holes_2 ranges
RESULT # the same as above

# aggregate intersection on Polygon_2 range
RESULT  # the same as above

# aggregate intersection on Polygon_with_holes_2 range
RESULT  # the same as above

# aggregate intersection on Polygon_2 and  Polygon_with_holes_2 ranges
RESULT  # the same as above

# aggregate symmetric difference on Polygon_with_holes_2 range
RESULT  # the same as above

# aggregate symmetric difference on Polygon_2 and  Polygon_with_holes_2 ranges
RESULT  # the same as above

# aggregate symmetric difference on Polygon_2 and  Polygon_with_holes_2 ranges
RESULT  # the same as in union


example:
--------

# range of Polygon_2
2 # two polygons
4
0/1 0/1
2/1 0/1
2/1 2/1
0/1 2/1

4
1/1 1/1
3/1 1/1
3/1 3/1
1/1 3/1


# range pf Polygon_with_holes_2
1 # one polygon with holes
0
1/1 1/1
2/1 1/1
2/1 2/1
1/1 2/1


# union result(1)
1
8
0/1 0/1
2/1 0/1
2/1 1/1
3/1 1/1
3/1 3/1
1/1 3/1
1/1 2/1
0/1 2/1
0

#union result(2)
1 # one polygon with holes
0
1/1 1/1
2/1 1/1
2/1 2/1
1/1 2/1

#union result(3)
1
0
0

#intersection(1)
1
4
1/1 1/1
2/1 1/1
2/1 2/1
1/1 2/1
0

#intersection(2)
1 # one polygon with holes
0
1/1 1/1
2/1 1/1
2/1 2/1
1/1 2/1

#intersection(3)
1
8
0/1 0/1
2/1 0/1
2/1 1/1
3/1 1/1
3/1 3/1
1/1 3/1
1/1 2/1
0/1 2/1
4
1/1 1/1
1/1 2/1
2/1 2/1
2/1 1/1

#symm diff result(1)
1
8
0/1 0/1
2/1 0/1
2/1 1/1
3/1 1/1
3/1 3/1
1/1 3/1
1/1 2/1
0/1 2/1
4
1/1 1/1
1/1 2/1
2/1 2/1
2/1 1/1


#symm diff result(2)
1 # one polygon with holes
0
1/1 1/1
2/1 1/1
2/1 2/1
1/1 2/1

#symm diff result(3)
1
0
8
0/1 0/1
0/1 2/1
1/1 2/1
1/1 3/1
3/1 3/1
3/1 1/1
2/1 1/1
2/1 0/1

