#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Timer.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>

#include <fstream>
#include <list>
#include <cassert>


typedef double coord_type;

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;
typedef CGAL::Periodic_3_triangulation_traits_3<EPIC> K;

typedef K::FT         FT;
typedef K::Point_3    Point_3;
typedef K::Vector_3   Vector_3;
typedef K::Segment_3  Segment_3;
typedef K::Ray_3      Ray_3;
typedef K::Line_3     Line;
typedef K::Triangle_3 Triangle_3;
typedef K::Iso_cuboid_3 Iso_cuboid_3;
typedef K::Tetrahedron_3 Tetrahedron_3;

typedef CGAL::Creator_uniform_3<double,Point_3> Creator;

typedef CGAL::Periodic_3_Delaunay_triangulation_3<K> P3DT3;

typedef P3DT3::Cell  Cell;
typedef P3DT3::Vertex Vertex;
typedef P3DT3::Edge Edge;
typedef P3DT3::Facet Facet;
typedef P3DT3::Cell_handle  Cell_handle;
typedef P3DT3::Vertex_handle Vertex_handle;

typedef P3DT3::Cell_circulator  Cell_circulator;

typedef P3DT3::Locate_type Locate_type;

typedef P3DT3::Cell_iterator  Cell_iterator;
typedef P3DT3::Vertex_iterator  Vertex_iterator;
typedef P3DT3::Edge_iterator  Edge_iterator;
typedef P3DT3::Facet_iterator Facet_iterator;

typedef P3DT3::Periodic_point_iterator Periodic_point_iterator;
typedef P3DT3::Periodic_triangle_iterator Periodic_triangle_iterator;

typedef P3DT3::Periodic_tetrahedron_iterator Periodic_tetrahedron_iterator;

typedef CGAL::Triangulation_2<EPIC> T2;
typedef EPIC::Point_2 Point_2;
typedef EPIC::Triangle_2 Triangle_2;

typedef CGAL::Timer Timer;

#endif
