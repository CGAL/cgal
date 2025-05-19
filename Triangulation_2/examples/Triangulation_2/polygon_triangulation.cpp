#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_triangulation_2.h>
#include "polygon_triangulation_graphics_scene_options.h"

#include <iostream>
#include <unordered_map>

#include <boost/property_map/property_map.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>            Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CDT::Face_handle                                          Face_handle;
typedef CDT::Point                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;

int main( )
{
  //construct two non-intersecting nested polygons
  Polygon_2 polygon1;
  polygon1.push_back(Point(0,0));
  polygon1.push_back(Point(2,0));
  polygon1.push_back(Point(2,2));
  polygon1.push_back(Point(1,1.75));
  polygon1.push_back(Point(0,2));
  Polygon_2 polygon2;
  polygon2.push_back(Point(0.5,0.5));
  polygon2.push_back(Point(1.5,0.5));
  polygon2.push_back(Point(1.5,1.5));
  polygon2.push_back(Point(0.5,1.5));

  //Insert the polygons into a constrained triangulation
  CDT cdt;
  cdt.insert_constraint(polygon1.vertices_begin(), polygon1.vertices_end(), true);
  cdt.insert_constraint(polygon2.vertices_begin(), polygon2.vertices_end(), true);
  cdt.insert_constraint(Point(0.25, 0.25), Point(0.25, 1.75));

  std::unordered_map<Face_handle, bool> in_domain_map;
  boost::associative_property_map< std::unordered_map<Face_handle,bool> >
    in_domain(in_domain_map);

  //Mark facets that are inside the domain bounded by the polygon
  CGAL::mark_domain_in_triangulation(cdt, in_domain);

  unsigned int count=0;
  for (Face_handle f : cdt.finite_face_handles())
  {
    if ( get(in_domain, f) ) ++count;
  }

  std::cout << "There are " << count << " faces in the domain." << std::endl;
  assert(count > 0);
  assert(count < cdt.number_of_faces());

  Polygon_triangulation_gs_options<CDT> gso(in_domain);
  CGAL::draw(cdt, gso);
  return 0;
}
