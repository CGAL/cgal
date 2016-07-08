#include <CGAL/internal/Generic_random_point_generator.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/boost/graph/property_maps.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <iostream>
#include <fstream>

using namespace CGAL;

int
main( )
{
  typedef Simple_cartesian<double>                           R;
  typedef R::Point_3                                         Point;
  typedef R::FT                                              FT;
  typedef Polyhedron_3<R>                                    Polyhedron;
  typedef Surface_mesh<Point>                                Surface_mesh;
  typedef boost::property_map<Polyhedron,
      vertex_point_t>::type                                  Vertex_point_pmap;
  typedef boost::property_map<Surface_mesh,
      vertex_point_t>::type                                  Vertex_point_pmap_sm;
  typedef Triangle_from_face_descriptor_map<
      Polyhedron,Vertex_point_pmap>                          Generator;
  typedef Triangle_from_face_descriptor_map<
      Surface_mesh,Vertex_point_pmap_sm>                     Generator_SM;
  typedef Random_points_in_triangle_3<Point>                 Creator;
  typedef boost::graph_traits<Polyhedron>::face_descriptor   face_iterator;
  typedef boost::graph_traits<Surface_mesh>::face_descriptor face_iterator_sm;



  std::vector<Point> points;
  Polyhedron poly;
  Surface_mesh sm;
  std::ifstream in("../../../Polyhedron/demo/Polyhedron/data/star.off");
  in >> sm;
  CGAL_assertion(in && !sm.is_empty());

  poly.make_tetrahedron(Point(0.0,0.2,0.0), Point(-0.2,0.0,0.0), Point(0.2,0.0,0.0), Point(0.0,0.0,0.2));

  Random_points_on_triangle_mesh_3<Point, Surface_mesh>
      g(sm);
  CGAL::cpp11::copy_n( g, 3000, std::back_inserter(points));
  for (std::size_t i = 0; i<points.size(); ++i)
    std::cerr<<points[i]<<std::endl;
  return 0;
}
