#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Vector_3.h>
#include <CGAL/point_generators_3.h> 
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/Triangulation_2.h>

typedef CGAL::Simple_cartesian<double>  Rp;
typedef CGAL::Triangulation_euclidean_traits_xy_3<Rp>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
//typedef CGAL::Triangulation_2<Gt> Delaunay;

typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;  

typedef Rp::Point_3   TPoint_3;



#include <list>
#include <iostream>
#include <istream>
#include <fstream>

