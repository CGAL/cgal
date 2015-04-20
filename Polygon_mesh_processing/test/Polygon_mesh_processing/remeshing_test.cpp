

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <boost/foreach.hpp>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef boost::graph_traits<Mesh>::edge_descriptor   edge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

template<typename VPmap>
void test_edge_lengths(const Mesh& pmesh,
                       const double& high,
                       const VPmap& vpmap)
{
  double sqhigh = high * high;
  BOOST_FOREACH(edge_descriptor e, edges(pmesh))
  {
    vertex_descriptor v1 = target(halfedge(e, pmesh), pmesh);
    vertex_descriptor v2 = source(halfedge(e, pmesh), pmesh);
    double sql = CGAL::squared_distance(vpmap[v1], vpmap[v2]);
    if (sqhigh < sql)
      std::cout << "sqhigh = " << sqhigh << "\t sql = " << sql << std::endl;
    CGAL_assertion(sqhigh >= sql);
  }
}

int main()
{
  std::ifstream input("data/U.off");
  Mesh m;

  if (!input || !(input >> m)){
    std::cerr << "Error: can not read file.\n";
    return 1;
  }

  double target_edge_length = 0.01;
  double low = 4. / 5. * target_edge_length;
  double high = 4. / 3. * target_edge_length;

  CGAL::Polygon_mesh_processing::incremental_triangle_based_remeshing(m,
    faces(m),
    target_edge_length,
    CGAL::Polygon_mesh_processing::parameters::number_of_iterations(5));

  boost::property_map<Mesh, boost::vertex_point_t>::const_type vpmap
    = boost::get(CGAL::vertex_point, m);

//  test_edge_lengths(m, high, vpmap);

  std::ofstream out("U_remeshed.off");
  out << m;
  out.close();

  return 0;
}