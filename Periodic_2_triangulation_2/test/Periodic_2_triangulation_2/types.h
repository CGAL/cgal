#ifndef P2T2_UNIT_TEST_TYPES_H
#define P2T2_UNIT_TEST_TYPES_H

#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>

using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
//struct K : public Exact_predicates_inexact_constructions_kernel {};
//typedef Simple_cartesian<Gmpq> K;
//typedef Simple_cartesian<double> K;

typedef Periodic_2_triangulation_traits_2<K> Gt;

typedef Gt::Point_2                                 Point;
typedef Gt::Vector_2                                Vector;
typedef Gt::Segment_2                               Segment;
typedef Gt::Triangle_2                              Triangle;
typedef Gt::Iso_rectangle_2                         Iso_rectangle;

#ifdef CGAL_USE_DELAUNAY
typedef Periodic_2_Delaunay_triangulation_2<Gt>     Triangulation;
#else
typedef Periodic_2_triangulation_2<Gt>              Triangulation;
#endif // CGAL_USE_DELAUNAY

typedef Triangulation::Offset                       Offset;
typedef Triangulation::Vertex_circulator            Vertex_circulator;
typedef Triangulation::Vertex_handle                Vertex_handle;
typedef Triangulation::Face_handle                  Face_handle;
typedef Triangulation::Vertex_iterator              Vertex_iterator;
typedef Triangulation::Face_iterator                Face_iterator;
typedef Triangulation::Periodic_point_iterator      Periodic_point_iterator;
typedef Triangulation::Periodic_segment_iterator    Periodic_segment_iterator;
typedef Triangulation::Periodic_triangle_iterator   Periodic_triangle_iterator;

typedef Periodic_2_Delaunay_triangulation_2<Gt>     Delaunay_triangulation;

typedef CGAL::Creator_uniform_2<double,Point>             Creator;
typedef CGAL::Random_points_in_square_2<Point, Creator>   Random_points_in_square;
typedef CGAL::Random_points_on_circle_2<Point, Creator>   Random_points_on_circle;
#endif // P2T2_UNIT_TEST_TYPES_H
