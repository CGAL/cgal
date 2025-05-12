#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_3/predicates.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef K::Point_3                                               Point_3;
typedef CGAL::Surface_mesh<Point_3>                              Mesh;
typedef Mesh::Property_map<Mesh::Vertex_index, Point_3>          PointMap;

int main(int argc, char* argv[])
{
  const std::string f1 = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");
  const std::string f2 = (argc>2) ? argv[2] : CGAL::data_file_path("meshes/sphere.off");
  Mesh sm1, sm2;
  if(!CGAL::IO::read_polygon_mesh(f1, sm1))
  {
    std::cerr<< "Cannot open " << f1 <<std::endl;
    return 1;
  }
  if(!CGAL::IO::read_polygon_mesh(f2, sm2))
  {
    std::cerr<< "Cannot open " << f2 <<std::endl;
    return 1;
  }

  //This will contain the extreme vertices
  std::vector<Mesh::Vertex_index> extreme_vertices1, extreme_vertices2;

  //call the function with the traits adapter for vertices
  CGAL::extreme_points_3(vertices(sm1), std::back_inserter(extreme_vertices1),
                         CGAL::make_extreme_points_traits_adapter(sm1.points()));
  CGAL::extreme_points_3(vertices(sm2), std::back_inserter(extreme_vertices2),
                         CGAL::make_extreme_points_traits_adapter(sm2.points()));

  //print the number of extreme vertices
  std::cout << "There are " << extreme_vertices1.size() << " and "
                            << extreme_vertices2.size()
                            << " extreme vertices in the input meshes.\n";

  bool res = CGAL::Convex_hull_3::do_intersect<K>(extreme_vertices1, extreme_vertices2,
                                               CGAL::Convex_hull_3::predicates_impl::Spherical_disjoint_traits_with_point_maps<K, PointMap>(sm1.points(), sm2.points()));

  std::cout << "do convex hulls intersect? " << res << "\n";

  return 0;
}

/*
typedef typename K::Point_3                                               Point_3;
  typedef CGAL::Surface_mesh<Point_3>                              Mesh;

  std::vector<typename K::Point_3> pts1, pts2;
  std::vector<boost::container::small_vector<std::size_t, 3>> trs1, trs2;
  if (!CGAL::IO::read_polygon_soup(f1, pts1, trs1))
  {
    std::cerr << "Cannot read " << f1 << "\n";
  }
  if (!CGAL::IO::read_polygon_soup(f2, pts2, trs2))
  {
    std::cerr << "Cannot read " << f2 << "\n";
  }

  CGAL::Real_timer t;
  t.start();
  Mesh hull1, hull2;
  CGAL::convex_hull_3(pts1.begin(), pts1.end(), hull1);
  CGAL::convex_hull_3(pts2.begin(), pts2.end(), hull2);
  t.stop();
  std::cout << "Computing convex hulls: " << t.time() << " sec" << std::endl;
  std::cout << "Convex hull size:" << vertices(hull1).size() << ", " << vertices(hull2).size() << "\n" << std::endl;
  t.reset();

  t.start();
  CGAL::Convex_hull_3::do_intersect<K>(hull1, hull2);
  t.stop();
  std::cout << "Do intersect: " << t.time() << " sec\n" << std::endl;
  t.reset();

  t.start();
  std::array<Point_3, 8> obb1, obb2;
  CGAL::oriented_bounding_box(hull1, obb1, CGAL::parameters::use_convex_hull(false));
  CGAL::oriented_bounding_box(hull2, obb2, CGAL::parameters::use_convex_hull(false));
  t.stop();
  std::cout << "Computing Obbs: " << t.time() << " sec\n" << std::endl;
  t.reset();

  t.start();
  CGAL::Convex_hull_3::do_intersect<K>(obb1, obb2);
  t.stop();
  std::cout << "Do intersect with Obbs: " << t.time() << " sec\n" << std::endl;
*/