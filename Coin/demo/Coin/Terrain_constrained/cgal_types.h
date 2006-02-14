#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Vector_3.h>
#include <CGAL/point_generators_3.h> 
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Bbox_3.h> 
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>


typedef CGAL::Simple_cartesian<double>  Rp;
typedef CGAL::Point_2<Rp>               Point_2;
typedef CGAL::Segment_2<Rp>             Segment_2;
typedef CGAL::Triangulation_euclidean_traits_xy_3<Rp>  Gt;
typedef CGAL::Constrained_Delaunay_triangulation_2<Gt> Delaunay;
typedef Delaunay::Finite_edges_iterator Finite_edges_iterator;
typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;  
typedef Delaunay::Finite_vertices_iterator Finite_vertices_iterator;
typedef Delaunay::Vertex_handle         Vertex_handle;
typedef Rp::Point_3   TPoint_3;
typedef CGAL::Bbox_3  Bbox;


#include <list>
#include <iostream>
#include <istream>
#include <fstream>

