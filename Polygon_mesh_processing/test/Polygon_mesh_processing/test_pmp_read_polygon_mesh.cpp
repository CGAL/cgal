#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

#include <CGAL/Simple_cartesian.h>

#include <fstream>
#include <vector>

typedef CGAL::Simple_cartesian<double>                       Kernel;

template<class Mesh>
struct Custom_VPM
{
  typedef Custom_VPM<Mesh>                                      Self;

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

  typedef vertex_descriptor                                     key_type;
  typedef Kernel::Point_3                                       value_type;
  typedef value_type&                                           reference;
  typedef boost::lvalue_property_map_tag                        category;

  Custom_VPM(std::map<key_type, value_type>& points) : points(points) { }

  friend void put(const Self& m, const key_type& k, const value_type& v) { m.points[k] = value_type(v.x(), v.y(), v.z()); }
  friend reference get(const Self& m, const key_type& k) { return m.points[k]; }

  std::map<key_type, value_type>& points;
};

template< class Mesh>
void test(const std::string& filename, const bool is_pm)
{
  std::cout << "Test " << filename << " with Mesh = " << typeid(Mesh).name() << " is PM? " << is_pm << std::endl;

  Mesh g;
  std::cout<<"0"<<std::endl;
  bool success = CGAL::IO::read_polygon_mesh(filename, g, CGAL::parameters::verbose(true));
  std::cout<<"1"<<std::endl;
  assert(is_pm == success); // if it's a pm, BGL reader should be enough

  g.clear();
  success = CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, g,
                                                                 CGAL::parameters::verbose(true)
                                                                                  .erase_all_duplicates(true));
  std::cout<<"2"<<std::endl;
  assert(success);
  assert(is_valid(g));
std::cout<<"3"<<std::endl;
  // Test VPM NP
  typename boost::property_map<Mesh,CGAL::vertex_point_t>::type vpm = get(CGAL::vertex_point, g);
  std::map<typename boost::graph_traits<Mesh>::vertex_descriptor, Kernel::Point_3> cpoints;
  Custom_VPM<Mesh> cvpm(cpoints);
std::cout<<"4"<<std::endl;
  Mesh g2;
  success = CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, g2,
                                                                 CGAL::parameters::vertex_point_map(cvpm)
                                                                                  .erase_all_duplicates(true));
  std::cout<<"5"<<std::endl;
  assert(success);
  assert(num_vertices(g) == num_vertices(g2));

  auto it = vertices(g).begin(), it2 = vertices(g2).begin();
  std::cout<<"6"<<std::endl;
  while(it != vertices(g).end() && it2 != vertices(g2).end())
  {
    assert(get(vpm, *it++) == get(cvpm, *it2++));
  }
  std::cout<<"7"<<std::endl;
}

int main()
{
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>                       Polyhedron;

  typedef Kernel::Point_3                                                                    Point;
  typedef CGAL::Surface_mesh<Point>                                                          SM;

  typedef CGAL::Linear_cell_complex_traits<3, Kernel>                                        LCC_traits;
  typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper<2, 3, LCC_traits>::type LCC;

  test<SM>(CGAL::data_file_path("meshes/pig.off"), true);
  test<Polyhedron>(CGAL::data_file_path("meshes/pig.off"), true);
  test<LCC>(CGAL::data_file_path("meshes/pig.off"), true);

  test<SM>("data_polygon_soup/nm_vertex_and_edge.off", false);
  test<Polyhedron>("data_polygon_soup/nm_vertex_and_edge.off", false);
  test<LCC>("data_polygon_soup/nm_vertex_and_edge.off", false);

  test<SM>("data_polygon_soup/incompatible_orientation.off", false);
  test<Polyhedron>("data_polygon_soup/incompatible_orientation.off", false);
  test<LCC>("data_polygon_soup/incompatible_orientation.off", false);

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
