#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Midpoint placement policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

//Placement wrapper
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

#include <iostream>
#include <fstream>
#include <map>

typedef CGAL::Simple_cartesian<double>                Kernel;
typedef Kernel::Point_3                               Point_3;
typedef CGAL::Polyhedron_3<Kernel>                    Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

// BGL property map which indicates whether an edge is marked as non-removable
struct Border_is_constrained_edge_map
{
  const Surface_mesh* sm_ptr;
  typedef boost::graph_traits<Surface_mesh>::edge_descriptor key_type;
  typedef bool                                               value_type;
  typedef value_type                                         reference;
  typedef boost::readable_property_map_tag                   category;

  Border_is_constrained_edge_map(const Surface_mesh& sm) : sm_ptr(&sm) {}

  friend bool get(Border_is_constrained_edge_map m, const key_type& edge) {
    return CGAL::is_border(edge, *m.sm_ptr);
  }
};

// Placement class
typedef SMS::Constrained_placement<SMS::Midpoint_placement<Surface_mesh>,
                                   Border_is_constrained_edge_map > Placement;

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const char* filename = (argc > 1) ? argv[1] : "data/mesh_with_border.off";
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

  // map used to check that constrained_edges and the points of its vertices
  // are preserved at the end of the simplification
  std::map<Surface_mesh::Halfedge_handle, std::pair<Point_3, Point_3> >constrained_edges;
  std::size_t nb_border_edges=0;

  for(Surface_mesh::Halfedge_iterator hit=surface_mesh.halfedges_begin(),
                                      hit_end=surface_mesh.halfedges_end();
                                      hit!=hit_end; ++hit)
  {
    if(hit->is_border())
    {
      constrained_edges[hit] = std::make_pair(hit->opposite()->vertex()->point(),
                                              hit->vertex()->point());
      ++nb_border_edges;
    }
  }

  // Contract the surface mesh as much as possible
  SMS::Count_stop_predicate<Surface_mesh> stop(0);
  Border_is_constrained_edge_map bem(surface_mesh);

  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface mesh lack an "id()" field.
  std::cout << "Collapsing as many edges of mesh: " << filename << " as possible..." << std::endl;
  int r = SMS::edge_collapse(surface_mesh,
                             stop,
                             CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,surface_mesh))
                                              .halfedge_index_map(get(CGAL::halfedge_external_index  ,surface_mesh))
                                              .edge_is_constrained_map(bem)
                                              .get_placement(Placement(bem)));

  std::cout << "\nFinished!\n" << r << " edges removed.\n"
            << (surface_mesh.size_of_halfedges()/2) << " final edges.\n";

  std::ofstream os(argc > 2 ? argv[2] : "out.off");
  os.precision(17);
  os << surface_mesh;

  // now check!
  for(Surface_mesh::Halfedge_iterator hit=surface_mesh.halfedges_begin(),
                                      hit_end=surface_mesh.halfedges_end();
                                      hit!=hit_end; ++hit)
  {
    if(hit->is_border())
    {
      --nb_border_edges;
      assert(constrained_edges[hit] == std::make_pair(hit->opposite()->vertex()->point(),
                                                      hit->vertex()->point()));
    }
  }

  assert(nb_border_edges == 0);

  return EXIT_SUCCESS;
}
