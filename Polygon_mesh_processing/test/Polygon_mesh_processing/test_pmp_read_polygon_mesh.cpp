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
void do_test()
{
  std::string name = "data_polygon_soup/bad_cube.off";
  Mesh g;
  CGAL_assertion(!CGAL::read_polygon_mesh(name, g));
  CGAL_assertion(CGAL::Polygon_mesh_processing::read_polygon_mesh(name, g));
  CGAL_assertion(is_valid(g));
  typedef typename boost::property_map<Mesh,CGAL::vertex_point_t>::type VertexPointMap;

  VertexPointMap vpm = get(CGAL::vertex_point, g);

  std::map<typename boost::graph_traits<Mesh>::vertex_descriptor, Kernel::Point_3> cpoints;
  Custom_VPM<Mesh> cvpm(cpoints);
  //test pmap param
  Mesh g2;
  CGAL_assertion(CGAL::Polygon_mesh_processing::read_polygon_mesh(name, g2, CGAL::parameters::vertex_point_map(cvpm)));
  CGAL_assertion(num_vertices(g2)==12);

  auto it = vertices(g2).begin(),
      it2 = vertices(g).begin();

  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));
  CGAL_assertion(get(cvpm, *(it++)) == get(vpm, *(it2++)));

}

int main()
{
  typedef Kernel::Point_3                                      Point;

  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

  typedef CGAL::Surface_mesh<Point>                            SM;

  typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
  typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
      <2, 3, MyTraits>::type LCC;


  do_test<Polyhedron>();
  do_test<SM>();
  do_test<LCC>();
  return 0;
}
