#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Vector_3.h>
#include <CGAL/point_generators_3.h> 
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>


#include <list>
#include <iostream>
#include <istream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>      Kernel;

typedef CGAL::Plane_3<Kernel>               Plane_3;
typedef CGAL::Point_3<Kernel>               Point_3;
typedef CGAL::Vector_3<Kernel>              Vector_3;
typedef CGAL::Tetrahedron_3<Kernel>         Tetrahedron_3;
typedef CGAL::Aff_transformation_3<Kernel>  Aff;

typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<double> > my_K;
struct  K : public my_K {};
typedef K::Point_3		Point;
typedef K::Segment_3	Segment;

typedef CGAL::Delaunay_triangulation_3<K>   Triangulation;
typedef Triangulation::Cell_handle		      Cell_handle;
typedef Triangulation::Vertex_handle	      Vertex_handle;
typedef Triangulation::Locate_type		      Locate_type;
typedef Triangulation::Edge_iterator	      Edge_iterator;
typedef Triangulation::Vertex_iterator	    Vertex_iterator;
typedef Triangulation::Facet_iterator       Facet_iterator;
typedef Triangulation::Finite_vertices_iterator
                                            Finite_vertices_iterator;
typedef Triangulation::Finite_cells_iterator
                                            Finite_cells_iterator;
typedef Triangulation::Geom_traits          Traits;


