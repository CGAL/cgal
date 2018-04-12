#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/stitch_holes.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;


void test_connected_components(const char* fname)
{
  std::ifstream input(fname);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  typedef typename boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;


  std::vector<std::set<halfedge_descriptor> > connected_components;

  CGAL::Polygon_mesh_processing::extract_connected_components(mesh,
                             std::back_inserter(connected_components));

  std::cout << "# cc = " << connected_components.size();

  for(auto set : connected_components)
  {
    std::cout << "of size= " << set.size() << std::endl;
  }
}


void test_merge_points(const char* fname)
{
  std::ifstream input(fname);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  typedef typename boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
  std::vector<vertex_descriptor> verts(vertices(mesh).begin(), vertices(mesh).end());

  vertex_descriptor v_rm = verts[1];
  vertex_descriptor v_keep = verts[2];

  CGAL::Polygon_mesh_processing::merge_identical_points(mesh, v_keep, v_rm);

  std::ofstream out("/tmp/result.off");
  out << mesh;
  out.close();
}


int main()
{


  test_connected_components("data/small_ex.off");
  test_merge_points("data/merge_points.off");



  return 0;
}
