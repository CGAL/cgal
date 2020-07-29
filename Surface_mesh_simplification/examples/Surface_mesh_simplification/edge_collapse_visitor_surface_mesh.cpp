#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                                  Kernel;
typedef Kernel::Point_3                                                 Point_3;

typedef CGAL::Surface_mesh<Point_3>                                     Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor          halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor            vertex_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;

typedef SMS::Edge_profile<Surface_mesh>                                 Profile;

// The following is a Visitor that keeps track of the simplification process.
// In this example the progress is printed real-time and a few statistics are
// recorded (and printed in the end).
//
struct Stats
{
  std::size_t collected = 0;
  std::size_t processed = 0;
  std::size_t collapsed = 0;
  std::size_t non_collapsable = 0;
  std::size_t cost_uncomputable = 0;
  std::size_t placement_uncomputable = 0;
};

struct My_visitor : SMS::Edge_collapse_visitor_base<Surface_mesh>
{
  My_visitor(Stats* s) : stats(s) {}

  // Called during the collecting phase for each edge collected.
  void OnCollected(const Profile&, const boost::optional<double>&)
  {
    ++(stats->collected);
    std::cerr << "\rEdges collected: " << stats->collected << std::flush;
  }

  // Called during the processing phase for each edge selected.
  // If cost is absent the edge won't be collapsed.
  void OnSelected(const Profile&,
                  boost::optional<double> cost,
                  std::size_t initial,
                  std::size_t current)
  {
    ++(stats->processed);
    if(!cost)
      ++(stats->cost_uncomputable);

    if(current == initial)
      std::cerr << "\n" << std::flush;
    std::cerr << "\r" << current << std::flush;
  }

  // Called during the processing phase for each edge being collapsed.
  // If placement is absent the edge is left uncollapsed.
  void OnCollapsing(const Profile&,
                    boost::optional<Point> placement)
  {
    if(!placement)
      ++(stats->placement_uncomputable);
  }

  // Called for each edge which failed the so called link-condition,
  // that is, which cannot be collapsed because doing so would
  // turn the surface mesh into a non-manifold.
  void OnNonCollapsable(const Profile&)
  {
    ++(stats->non_collapsable);
  }

  // Called after each edge has been collapsed
  void OnCollapsed(const Profile&, vertex_descriptor)
  {
    ++(stats->collapsed);
  }

  Stats* stats;
};

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const char* filename = (argc > 1) ? argv[1] : "data/cube.off";
  std::ifstream is(filename);
  if(!is || !(is >> surface_mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // In this example, the simplification stops when the number of undirected edges
  // drops below xx% of the initial count
  const double ratio = (argc > 2) ? std::stod(argv[2]) : 0.1;
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);

  Stats stats;
  My_visitor vis(&stats);

  // The index maps are not explicitelty passed as in the previous
  // example because the surface mesh items have a proper id() field.
  // On the other hand, we pass here explicit cost and placement
  // function which differ from the default policies, ommited in
  // the previous example.
  int r = SMS::edge_collapse(surface_mesh, stop, CGAL::parameters::visitor(vis));

  std::cout << "\nEdges collected: "  << stats.collected
            << "\nEdges proccessed: " << stats.processed
            << "\nEdges collapsed: "  << stats.collapsed
            << std::endl
            << "\nEdges not collapsed due to topological constraints: "  << stats.non_collapsable
            << "\nEdge not collapsed due to cost computation constraints: "  << stats.cost_uncomputable
            << "\nEdge not collapsed due to placement computation constraints: " << stats.placement_uncomputable
            << std::endl;

  std::cout << "\nFinished!\n" << r << " edges removed.\n"
            << surface_mesh.number_of_edges() << " final edges.\n";

  std::ofstream os(argc > 3 ? argv[3] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
