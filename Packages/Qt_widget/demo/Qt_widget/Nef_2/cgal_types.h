#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_2.h>



//Cartesian types
typedef CGAL::Cartesian<CGAL::Gmpq>   Rp;
typedef CGAL::Partition_traits_2<Rp>
                                  Traits;
typedef Rp::Point_2               Cartesian_point_2;
typedef Rp::Line_2                Cartesian_line_2;
typedef Traits::Polygon_2         Cartesian_polygon_2;
typedef Cartesian_polygon_2::Vertex_iterator
                                  Vertex_iterator;

typedef CGAL::Cartesian<double>   RP_double;
typedef CGAL::Partition_traits_2<RP_double>
                                  Traits_double;
typedef Traits_double::Polygon_2  Polygon_2_double;
typedef Polygon_2_double::Vertex_iterator
                                  Vertex_iterator_double;

//The Nef_Polyhedron types
typedef CGAL::Gmpz                RT;
typedef CGAL::Filtered_extended_homogeneous<RT>
                                  Extended_kernel;
typedef CGAL::Nef_polyhedron_2<Extended_kernel>
                                  Nef_polyhedron;
typedef Nef_polyhedron::Point     Point;
typedef Nef_polyhedron::Line      Line;

