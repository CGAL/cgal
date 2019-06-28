#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>              Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>        CDT;

typedef CDT::Point                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>                  Mesh_2_criteria;

using namespace CGAL;
int main()
{
 // Generated points are in that vector
  std::vector<Point> points;

  //Construct two non-intersecting nested polygons
  ::Polygon_2 polygon1;
  polygon1.push_back(Point(0,0));
  polygon1.push_back(Point(2,0));
  polygon1.push_back(Point(2,2));
  polygon1.push_back(Point(0,2));
  ::Polygon_2 polygon2;
  polygon2.push_back(Point(4.0,-2.0));
  polygon2.push_back(Point(4.0,2.0));
  polygon2.push_back(Point(6.0,0.0));

  //Insert the polygons into a constrained triangulation
  CDT cdt;
  cdt.insert_constraint(polygon1.vertices_begin(), polygon1.vertices_end(), true);
  cdt.insert_constraint(polygon2.vertices_begin(), polygon2.vertices_end(), true);

  // Refine the triangulation (and mark the faces as inside/outside)
  CGAL::refine_Delaunay_mesh_2(cdt, Mesh_2_criteria(0.125, 0.5));

  // Create the generator, input is the Triangulation_2 cdt
  Random_points_in_triangle_mesh_2<Point, CDT> g(cdt);

  // Get 100 random points in cdt
  std::copy_n(g, 100, std::back_inserter(points));

  // Check that we have really created 100 points.
  assert(points.size() == 100);

  // print the first point that was generated
  std::cout << points[0] << std::endl;

  return 0;
}
