#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Timer.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>



typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;

typedef EPIC::Vector_3 Vector_3;
typedef EPIC::Plane_3  Plane_3;
typedef EPIC::Line_3   Line_3;
typedef EPIC::Iso_cuboid_3   Iso_cuboid_3;

typedef CGAL::Projection_traits_xy_3<EPIC>   K;

typedef CGAL::Triangulation_vertex_base_2<K>                     Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDTbase;
typedef CGAL::Constrained_triangulation_plus_2<CDTbase> Terrain;

typedef Terrain::Point                    Point_3;
typedef Terrain::Vertex_handle            Vertex_handle;
typedef Terrain::Finite_vertices_iterator Finite_vertices_iterator;
typedef Terrain::Finite_faces_iterator    Finite_faces_iterator;
typedef Terrain::Finite_edges_iterator    Finite_edges_iterator;
typedef Terrain::Vertex_circulator        Vertex_circulator;
typedef Terrain::Vertex_iterator          Vertex_iterator;
typedef Terrain::Face_handle              Face_handle;

typedef  CGAL::Timer Timer;
typedef  CGAL::Bbox_3 Bbox_3;

struct Scene {

  Terrain terrain;
  Bbox_3 bbox;

};

#endif  // TYPEDEFS_H
