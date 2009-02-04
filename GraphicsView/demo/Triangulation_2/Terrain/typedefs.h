#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
//#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Timer.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>


typedef CGAL::Simple_cartesian<double> Construction_kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel EPEC;
typedef CGAL::Cartesian_converter<EPEC, Construction_kernel> Converter;

typedef EPEC::Vector_3 Vector_3;
typedef EPEC::Plane_3  Plane_3;
typedef EPEC::Line_3   Line_3;
typedef EPEC::Iso_cuboid_3   Iso_cuboid_3;

//typedef CGAL::Projection_traits_xy_3<EPEC>   K;
typedef CGAL::Triangulation_2_projection_traits_3<EPEC>   Traits;

typedef CGAL::Triangulation_vertex_base_2<Traits>                Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Traits>      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
// typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Exact_intersections_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS,
                                                   Itag>         CDTbase;
typedef CGAL::Constrained_triangulation_plus_2<CDTbase> Terrain;
// typedef CDTbase Terrain;

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
  Scene()
    : terrain(Traits(EPEC::Vector_3(0.,0.,1.))) 
  {
  }

  void setNormal(Vector_3 v) {
    terrain = Terrain(Traits(v));
    refresh();
  }

  void refresh() {
    std::cerr << "Recompute the terrain...\n";
    terrain.clear();
    terrain.insert(points.begin(), points.end());
    for(unsigned int line_id = 0; line_id < polylines.size(); ++line_id) 
    {
      Vertex_handle vh = terrain.insert(polylines[line_id][0]);
      for(unsigned int i = 1; i < polylines[line_id].size(); ++i)
      {
        Vertex_handle wh = terrain.insert(polylines[line_id][i], vh->face());
        terrain.insert_constraint(vh, wh);
        vh = wh;
      }
    }
  }


  std::vector<Point_3> points;
  std::vector<std::vector<Point_3> > polylines;
  Terrain terrain;
  Bbox_3 bbox;
};

#endif  // TYPEDEFS_H
