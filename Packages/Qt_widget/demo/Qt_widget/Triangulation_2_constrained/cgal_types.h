#ifndef CGAL_CGAL_TYPES_H
#define CGAL_CGAL_TYPES_H

#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/point_generators_2.h>


typedef double Coord_type;
typedef CGAL::Simple_cartesian<Coord_type>      K1;
typedef CGAL::Filtered_kernel<K1>               K2;
struct Rep : public K2{} ;

typedef Rep::Point_2                            Point;
typedef Rep::Segment_2                          Segment;
typedef Rep::Line_2                             Line;
typedef Rep::Triangle_2                         Triangle;
typedef Rep::Circle_2                           Circle;

typedef CGAL::Triangulation_vertex_base_2<Rep>  Vb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vb>  Vb1;
typedef CGAL::Constrained_triangulation_face_base_2<Rep>
                                                Fb1;
typedef CGAL::Triangulation_data_structure_2<Vb1, Fb1>
                                                TDS;
typedef CGAL::Exact_predicates_tag              Itag;

typedef CGAL::Constrained_Delaunay_triangulation_2<Rep, TDS, Itag>
                                                CT;
typedef CGAL::Triangulation_hierarchy_2<CT>     CDT;
typedef CDT::Vertex_iterator                    Vertex_iterator;
typedef CDT::Constraint                         Constraint;
typedef CDT::Vertex_handle                      Vertex_handle;
typedef CGAL::Partition_traits_2<Rep>           Traits;
typedef Traits::Polygon_2                       Cgal_Polygon;

#endif

