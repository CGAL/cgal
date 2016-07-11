#include <CGAL/point_generators_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <fstream>

using namespace CGAL;
void test_volume_mesh()
{
  typedef Simple_cartesian<double>                           R;
  typedef R::Point_3                                         Point;
  typedef R::FT                                              FT;
  typedef Surface_mesh<Point>                                Surface_mesh;

  std::vector<Point> points;
  Surface_mesh sm;
  std::ifstream in("../../../Polyhedron/demo/Polyhedron/data/star.off");
  in >> sm;
  CGAL_assertion(in && !sm.is_empty());


  Random_points_on_triangle_mesh_3<Point, Surface_mesh>
      g(sm);
  CGAL::cpp11::copy_n( g, 3000, std::back_inserter(points));
  for (std::size_t i = 0; i<points.size(); ++i)
    std::cerr<<points[i]<<std::endl;
}


typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>              Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>        CDT;
typedef CDT::Point                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;



void test_T2()
{
  std::vector<Point> points;
//construct two non-intersecting nested polygons
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

Random_points_on_triangle_mesh_2<Point, CDT>
    g(cdt);
CGAL::cpp11::copy_n( g, 300, std::back_inserter(points));
for (std::size_t i = 0; i<points.size(); ++i)
  std::cerr<<points[i]<<std::endl;


}
int
main( )
{
  test_volume_mesh();
  //test_T2();

  return 0;
}
