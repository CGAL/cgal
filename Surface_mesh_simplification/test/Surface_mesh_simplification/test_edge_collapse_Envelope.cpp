#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Timer.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>


#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Polyhedral_envelope_filter.h>

//bbox
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <iostream>


namespace SMS = CGAL::Surface_mesh_simplification;

typedef CGAL::Simple_cartesian<double>                        Kernel;

typedef Kernel::Point_3                                       Point_3;
typedef CGAL::Surface_mesh<Point_3>                           Surface;

typedef SMS::LindstromTurk_cost<Surface>                      Cost;
typedef SMS::LindstromTurk_placement<Surface>                 Placement;
typedef SMS::Polyhedral_envelope_filter<Kernel,SMS::Bounded_normal_change_filter<> > Filter;

struct Stats
{
  std::size_t collected = 0;
  std::size_t processed = 0;
  std::size_t collapsed = 0;
  std::size_t non_collapsable = 0;
  std::size_t cost_uncomputable = 0;
  std::size_t placement_uncomputable = 0;
};

struct My_visitor : SMS::Edge_collapse_visitor_base<Surface>
{
  My_visitor(Stats* s) : stats(s) {}

  // Called during the collecting phase for each edge collected.
  void OnCollected(const Profile&, const boost::optional<double>&)
  {
    ++(stats->collected);
    std::cerr << "\rEdges collected: " << stats->collected << std::endl;
  }

  // Called during the processing phase for each edge selected.
  // If cost is absent the edge won't be collapsed.
  void OnSelected(const Profile&,
                  boost::optional<double> cost,
                  std::size_t /* initial */,
                  std::size_t /* current */)
  {
    ++(stats->processed);
    if(!cost)
      ++(stats->cost_uncomputable);
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
  Surface input_mesh;
  std::ifstream is(argc > 1 ? argv[1] : "data/helmet.off");
  is >> input_mesh;

  SMS::Count_stop_predicate<Surface> stop(0); // go as far as you can while in the envelope

  Stats stats;
  My_visitor vis(&stats);

  std::cout << "Input has " << num_vertices(input_mesh) << " vertices and " << num_edges(input_mesh) << " edges" << std::endl;
  CGAL::Iso_cuboid_3<Kernel> bbox(CGAL::Polygon_mesh_processing::bbox(input_mesh));

  Point_3 cmin = (bbox.min)();
  Point_3 cmax = (bbox.max)();
  const double diag = CGAL::approximate_sqrt(CGAL::squared_distance(cmin, cmax));

  Surface mesh_cpy = input_mesh; // need a copy to keep the AABB tree valid
  Surface small_envelope_mesh = input_mesh;
  Surface big_envelope_mesh = input_mesh;
  Surface huge_envelope_mesh = input_mesh;

  CGAL::Timer t;
  t.start();

  {
    std::cout << "eps = " << 0.005*diag << std::endl;
    Placement placement;
    Filter filter(0.005*diag);
    SMS::edge_collapse(small_envelope_mesh, stop, CGAL::parameters::get_cost(Cost()).filter(filter).get_placement(placement));
    std::cout << "Output has " << vertices(small_envelope_mesh).size() << " vertices and " << edges(small_envelope_mesh).size() << " edges" << std::endl;
    std::cout << t.time() << "sec\n";
    t.reset();
  }

  {
    std::cout << "eps = " << 0.01*diag << std::endl;
    Placement placement;
    Filter filter(0.01*diag);
    SMS::edge_collapse(big_envelope_mesh, stop, CGAL::parameters::get_cost(Cost()).filter(filter).visitor(vis).get_placement(placement));
    std::cout << "Output has " << vertices(big_envelope_mesh).size() << " vertices and " << edges(big_envelope_mesh).size() << " edges" << std::endl;
    std::cout << t.time() << "sec\n";
    t.reset();
  }

  {
    std::cout << "eps = " << 0.02*diag << std::endl;
    Placement placement;
    Filter filter(0.02*diag);
    SMS::edge_collapse(huge_envelope_mesh, stop, CGAL::parameters::get_cost(Cost()).filter(filter).get_placement(placement));
    std::cout << "Output has " << vertices(huge_envelope_mesh).size() << " vertices and " << edges(huge_envelope_mesh).size() << " edges" << std::endl;
    std::cout << t.time() << "sec\n";
  }



  std::cout << "\nEdges collected: "  << stats.collected
            << "\nEdges proccessed: " << stats.processed
            << "\nEdges collapsed: "  << stats.collapsed
            << std::endl
            << "\nEdges not collapsed due to topological constraints: "  << stats.non_collapsable
            << "\nEdge not collapsed due to cost computation constraints: "  << stats.cost_uncomputable
            << "\nEdge not collapsed due to placement computation constraints: " << stats.placement_uncomputable
            << std::endl;

  assert(vertices(input_mesh).size() > vertices(small_envelope_mesh).size());
  assert(vertices(small_envelope_mesh).size() > vertices(big_envelope_mesh).size());
  assert(vertices(big_envelope_mesh).size() > vertices(huge_envelope_mesh).size());

  return EXIT_SUCCESS;
}
