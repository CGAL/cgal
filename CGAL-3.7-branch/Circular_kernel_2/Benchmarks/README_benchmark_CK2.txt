The usage is:
./example ${alpha} ${beta}

where alpha:
1: means to bench the BBox(CK) with Vartraits
2: means to bench the Lazy(CK) with Vartraits
3: means to bench the CK with Vartraits
4: means to bench the Bbox(Lazy(CK)) with Vartraits
5: means to bench the BBox(CK) with Circulartraits 
6: means to bench the Lazy(CK) with Circulartraits
7: means to bench the CK(CK) with Circulartraits
8: means to bench the Bbox(Lazy(CK)) Circulartraits
(le 0 est interne)

beta:
0: Compute the arrangement of DXF/51.dxf with the kernel ${alpha}
1: Compute the arrangement of DXF/cad_l1.dxf with the kernel ${alpha}
2: Compute the arrangement of DXF/cad_l2.dxf
3: Compute the arrangement of DXF/che_mod1.dxf
4: Compute the arrangement of DXF/CIOnZDraw.dxf
5: Compute the arrangement of DXF/mask1.dxf
6: Compute the arrangement of DXF/elekonta.dxf
7: Compute the arrangement of DXF/netlist_signal_1.dxf
8: Compute the arrangement of DXF/painttrack.dxf
9: Compute the arrangement of ${Scenario1}
a: Compute the arrangement of ${Scenario2}
b: Compute the arrangement of ${Scenario3}
c: Compute the arrangement of ${Scenario4}
d: Compute the arrangement of ${Scenario5}

${Scenario1}:
- Circles with center in [0,10]x[0,10], 
- 0.5 of distance between each circle
- unitary radius

${Scenario2}:
- Circles with center in [0,0.2]x[0,0.2],
- 0.01 of distance between each circle
- unitary radius

${Scenario3}:
Only one circle

${Scenario4}
100 Random Circles with
- Center in [0,1]x[0,1]
- Radius in [0,1]

${Scenario5}
Lattice, like:
o o o o o o o o o o o o o o o o
 o o o o o o o o o o o o o o o 
o o o o o o o o o o o o o o o o
 o o o o o o o o o o o o o o o 
o o o o o o o o o o o o o o o o
 o o o o o o o o o o o o o o o 
with no intersection.

The output:

In std::cout :
The number of elements to compute the arrangement
The number of circles and polygons (wich the side may be circular arcs)
The time needed to compute it,
The number of Vertices, Edges and Faces of the arrangement

In std::cerr :
Only the time needed to compute it. (it is useful to benchmark a lot of cases and redirect it on a .txt)

ATTENTION:
1) dont use ./example a b
with 5 <= a <= 8  and 0 <= b <= 8, we cannot use the Circulartraits to handle the files
2) The files have to be put on a folder name DXF where the program is located

Compile with -DCGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES if you want to benchmark with
the same kernel, but with an additional map for supporting circles.


