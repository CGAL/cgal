// Copyright (c) 2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

#ifndef P2T2_UNIT_TEST_TYPES_H
#define P2T2_UNIT_TEST_TYPES_H

#define _USE_MATH_DEFINES
#include <math.h>

#include <CGAL/Periodic_2_triangulation_2.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>

#include <CGAL/Periodic_2_triangulation_hierarchy_2.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/use.h>

using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
//struct K : public Exact_predicates_inexact_constructions_kernel {};
//typedef Simple_cartesian<Gmpq> K;
//typedef Simple_cartesian<double> K;

typedef Periodic_2_triangulation_traits_2<K> P2TT;
typedef Periodic_2_Delaunay_triangulation_traits_2<K> P2DTT;

#ifdef CGAL_USE_DELAUNAY
typedef Periodic_2_Delaunay_triangulation_2<P2DTT>    Triangulation;
#else
typedef Periodic_2_triangulation_2<P2TT>              Triangulation;
#endif // CGAL_USE_DELAUNAY

typedef Triangulation::Point                          Point;
typedef Triangulation::Vector                         Vector;
typedef Triangulation::Segment                        Segment;
typedef Triangulation::Triangle                       Triangle;
typedef Triangulation::Iso_rectangle                  Iso_rectangle;

typedef Triangulation::Offset                         Offset;
typedef Triangulation::Vertex_circulator              Vertex_circulator;
typedef Triangulation::Vertex_handle                  Vertex_handle;
typedef Triangulation::Face_handle                    Face_handle;
typedef Triangulation::Vertex_iterator                Vertex_iterator;
typedef Triangulation::Face_iterator                  Face_iterator;
typedef Triangulation::Periodic_point_iterator        Periodic_point_iterator;
typedef Triangulation::Periodic_segment_iterator      Periodic_segment_iterator;
typedef Triangulation::Periodic_triangle_iterator     Periodic_triangle_iterator;

typedef Periodic_2_Delaunay_triangulation_2<P2DTT>    Delaunay_triangulation;

typedef Periodic_2_triangulation_vertex_base_2<P2DTT>             PTH_Vbb;
typedef Triangulation_hierarchy_vertex_base_2<PTH_Vbb>            PTH_Vb;
typedef Periodic_2_triangulation_face_base_2<P2DTT>               PTH_Fb;
typedef Triangulation_data_structure_2<PTH_Vb, PTH_Fb>            PTH_Tds;
typedef Periodic_2_Delaunay_triangulation_2<P2DTT, PTH_Tds>       PTH_Dt;
typedef Periodic_2_triangulation_hierarchy_2<PTH_Dt>              Delaunay_triangulation_hierarchy;

typedef Creator_uniform_2<double, Point>              Creator;
typedef Random_points_in_square_2<Point, Creator>     Random_points_in_square;
typedef Random_points_on_circle_2<Point, Creator>     Random_points_on_circle;

#endif // P2T2_UNIT_TEST_TYPES_H
