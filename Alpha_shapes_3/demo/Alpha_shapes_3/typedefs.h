#ifndef TYPEDEFS_H
#define TYPEDEFS_H


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Timer.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>



typedef double coord_type;


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3    Point_3;
typedef K::Vector_3   Vector_3;
typedef K::Segment_3  Segment_3;
typedef K::Ray_3      Ray_3;
typedef K::Line_3     Line;
typedef K::Triangle_3 Triangle_3;
typedef K::Iso_cuboid_3 Iso_cuboid_3;

typedef CGAL::Alpha_shape_vertex_base_3<K> Vb;
typedef CGAL::Alpha_shape_cell_base_3<K>   Fb;

typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds> Triangulation_3;

typedef CGAL::Alpha_shape_3<Triangulation_3>  Alpha_shape_3;

typedef Alpha_shape_3::Cell  Cell;
typedef Alpha_shape_3::Vertex Vertex;
typedef Alpha_shape_3::Edge Edge;
typedef Alpha_shape_3::Facet Facet;
typedef Alpha_shape_3::Cell_handle  Cell_handle;
typedef Alpha_shape_3::Vertex_handle Vertex_handle;

typedef Alpha_shape_3::Cell_circulator  Cell_circulator;

typedef Alpha_shape_3::Locate_type Locate_type;

typedef Alpha_shape_3::Cell_iterator  Cell_iterator;
typedef Alpha_shape_3::Vertex_iterator  Vertex_iterator;
typedef Alpha_shape_3::Edge_iterator  Edge_iterator;


typedef Alpha_shape_3::Coord_type Coord_type;
typedef Alpha_shape_3::Alpha_iterator Alpha_iterator;

typedef CGAL::Timer Timer;

struct Scene {

  Alpha_shape_3 alpha_shape;

  std::list<Point_3> points;
};



#endif
