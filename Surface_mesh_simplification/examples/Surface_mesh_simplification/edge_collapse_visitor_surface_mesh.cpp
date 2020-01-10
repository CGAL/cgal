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
typedef Kernel::Point_2                                                 Point_2;

typedef CGAL::Surface_mesh<Point_3>                                     Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor          halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor            vertex_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;

typedef SMS::Edge_profile<Surface_mesh>                                 Profile;
typedef Surface_mesh::Property_map<vertex_descriptor, Point_2>          UV_pmap;

// The following is a Visitor that keeps track of the simplification process.



struct My_visitor : SMS::Edge_collapse_visitor_base<Surface_mesh>
{
  My_visitor(UV_pmap)
    : uv_pmap(uv_pmap)
  {}

  // Called during the processing phase for each edge being collapsed.
  // If placement is absent the edge is left uncollapsed.
  void OnCollapsing(const Profile& prof,
                    boost::optional<Point> placement)
  {
      if (placement) {
        p0 = prof.p0();
        p1 = prof.p1();
        vertex_descriptor v0 = prof.v0();
        vertex_descriptor v1 = prof.v1();
        p0_2 = get(uv_pmap, v0);
        p1_2 = get(uv_pmap, v1);
        p_2 = CGAL::midpoint(p0_2,p1_2);
      }
  }

  // Called after each edge has been collapsed
  void OnCollapsed(const Profile&, vertex_descriptor vd)
  {
    put(uv_pmap, vd, p_2);
  }

  UV_pmap uv_pmap;
  Point_3 p0, p1;
  Point_2 p0_2, p1_2, p_2;
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

  UV_pmap uv_pmap = surface_mesh.add_property_map<vertex_descriptor, Point_2>("v:uv").first;
  

  My_visitor vis(uv_pmap);


  int r = SMS::edge_collapse(surface_mesh, stop, CGAL::parameters::visitor(vis));

 

  return EXIT_SUCCESS;
}
