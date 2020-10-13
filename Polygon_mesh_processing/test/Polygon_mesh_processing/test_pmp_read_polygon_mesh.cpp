#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

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
  bool success = CGAL::read_polygon_mesh(filename, g, CGAL::parameters::verbose(true));
  CGAL_assertion(is_pm == success); // if it's a pm, BGL reader should be enough

  clear(g);
  success = CGAL::Polygon_mesh_processing::read_polygon_mesh(filename, g,
                                                             CGAL::parameters::verbose(true)
                                                                              .erase_all_duplicates(true));
  CGAL_assertion(success);
  CGAL_assertion(is_valid(g));

  // Test VPM NP
  typename boost::property_map<Mesh,CGAL::vertex_point_t>::type vpm = get(CGAL::vertex_point, g);
  std::map<typename boost::graph_traits<Mesh>::vertex_descriptor, Kernel::Point_3> cpoints;
  Custom_VPM<Mesh> cvpm(cpoints);

  Mesh g2;
  success = CGAL::Polygon_mesh_processing::read_polygon_mesh(filename, g2,
                                                             CGAL::parameters::vertex_point_map(cvpm)
                                                                              .erase_all_duplicates(true));
  CGAL_assertion(success);
  CGAL_assertion(num_vertices(g) == num_vertices(g2));

  auto it = vertices(g).begin(), it2 = vertices(g2).begin();
  while(it != vertices(g).end())
  {
    CGAL_assertion(get(vpm, *it++) == get(cvpm, *it2++));
  }
}

int main()
{
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>                       Polyhedron;

  typedef Kernel::Point_3                                                                    Point;
  typedef CGAL::Surface_mesh<Point>                                                          SM;

  typedef CGAL::Linear_cell_complex_traits<3, Kernel>                                        LCC_traits;
  typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper<2, 3, LCC_traits>::type LCC;

  test<SM>("data/pig.off", true);
  test<Polyhedron>("data/pig.off", true);
  test<LCC>("data/pig.off", true);

  test<SM>("data_polygon_soup/nm_vertex_and_edge.off", false);
  test<Polyhedron>("data_polygon_soup/nm_vertex_and_edge.off", false);
  test<LCC>("data_polygon_soup/nm_vertex_and_edge.off", false);

  test<SM>("data_polygon_soup/incompatible_orientation.off", false);
  test<Polyhedron>("data_polygon_soup/incompatible_orientation.off", false);
  test<LCC>("data_polygon_soup/incompatible_orientation.off", false);

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
