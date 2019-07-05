#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>
#include <CGAL/Shortest_noncontractible_cycle.h>
#include <CGAL/squared_distance_3.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;

struct Weight_functor {
  using Weight_t = double;
  Weight_functor(const Mesh& mesh) : m_mesh(mesh) {}
  double operator()(Mesh::Halfedge_index he) const {
    Point A = m_mesh.point(m_mesh.vertex(m_mesh.edge(he), 0));
    Point B = m_mesh.point(m_mesh.vertex(m_mesh.edge(he), 1));
    return CGAL::sqrt(CGAL::squared_distance(A, B));
  }
private:
  Mesh m_mesh;
};

using SNC = CGAL::Surface_mesh_topology::Shortest_noncontractible_cycle<Mesh, Weight_functor>;

int main(int argc, char* argv[]) {
  std::cout << "Program edgewidth_surface_mesh started.\n";
  Mesh sm;
  std::ifstream in ((argc > 1) ? argv[1] : "data/3torus-smooth.off");
  in >> sm;
  std::cout << "File loaded. Running the main program...\n";
  Weight_functor wf(sm);
  SNC snc(sm, wf);
  SNC::Path cycle;
  SNC::Distance_type x;
  std::cout << "Finding edge-width of the mesh...\n";
  snc.edge_width(cycle, &x);
  if (cycle.size() == 0) {
    std::cout << "  Cannot find edge-width. Stop.\n";
    return 0;
  }
  
  std::cout << "  Number of edges in cycle: " << cycle.size() << std::endl;
  std::cout << "  Cycle length: " << x << std::endl;
  std::cout << "  Root: " << sm.point(sm.vertex(sm.edge(cycle[0]), 0)) << std::endl;

  return EXIT_SUCCESS;
}