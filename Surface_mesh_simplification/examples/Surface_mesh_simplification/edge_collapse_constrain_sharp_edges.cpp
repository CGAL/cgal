#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/property_map.h>

#include <cmath>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Polyhedron_3<Kernel>                              Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor      edge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_iterator        edge_iterator;

namespace SMS = CGAL::Surface_mesh_simplification;

// BGL property map which indicates whether an edge is marked as non-removable
struct Constrained_edge_map
  : public boost::put_get_helper<bool, Constrained_edge_map>
{
  typedef boost::readable_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;

  Constrained_edge_map(const CGAL::Unique_hash_map<key_type,bool>& aConstraints)
    : mConstraints(aConstraints)
  {}

  reference operator[](const key_type& e) const { return  is_constrained(e); }

  bool is_constrained(const key_type& e) const { return mConstraints.is_defined(e); }

private:
  const CGAL::Unique_hash_map<key_type,bool>& mConstraints;
};

bool is_border (edge_descriptor e, const Surface_mesh& sm)
{
  return (face(halfedge(e,sm),sm) == boost::graph_traits<Surface_mesh>::null_face()) ||
         (face(opposite(halfedge(e,sm),sm),sm) == boost::graph_traits<Surface_mesh>::null_face());
}

Point_3 point(vertex_descriptor vd,  const Surface_mesh& sm)
{
  return get(CGAL::vertex_point, sm, vd);
}

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const char* filename = (argc > 1) ? argv[1] : "data/cube-meshed.off";
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

  CGAL::Unique_hash_map<edge_descriptor, bool> constraint_hmap(false);
  Constrained_edge_map constraints_map(constraint_hmap);
  SMS::Constrained_placement<SMS::Midpoint_placement<Surface_mesh>,
                             Constrained_edge_map > placement(constraints_map);

  // map used to check that constrained_edges and the points of its vertices
  // are preserved at the end of the simplification
  // Warning: the computation of the dihedral angle is only an approximation and can
  //          be far from the real value and could influence the detection of sharp
  //          edges after the simplification
  std::map<edge_descriptor,std::pair<Point_3, Point_3> >constrained_edges;
  std::size_t nb_sharp_edges = 0;

  // detect sharp edges
  std::ofstream cst_output("constrained_edges.cgal");
  for(edge_descriptor ed : edges(surface_mesh))
  {
    halfedge_descriptor hd = halfedge(ed, surface_mesh);
    if(is_border(ed, surface_mesh))
    {
      std::cerr << "border" << std::endl;
      ++nb_sharp_edges;
      constraint_hmap[ed] = true;
      constrained_edges[ed] = std::make_pair(point(source(hd, surface_mesh), surface_mesh),
                                             point(target(hd, surface_mesh), surface_mesh));
    }
    else
    {
      double angle = CGAL::approximate_dihedral_angle(point(target(opposite(hd, surface_mesh), surface_mesh), surface_mesh),
                                                      point(target(hd, surface_mesh), surface_mesh),
                                                      point(target(next(hd, surface_mesh), surface_mesh), surface_mesh),
                                                      point(target(next(opposite(hd, surface_mesh), surface_mesh), surface_mesh), surface_mesh));
      if(CGAL::abs(angle) < 100)
      {
        ++nb_sharp_edges;
        constraint_hmap[ed] = true;
        Point_3 p = point(source(hd, surface_mesh), surface_mesh);
        Point_3 q = point(target(hd, surface_mesh), surface_mesh);
        constrained_edges[ed] = std::make_pair(p,q);
        cst_output << "2 " << p << " "  << q << "\n";
      }
    }
  }
  cst_output.close();

  std::cerr << "# sharp edges = " << nb_sharp_edges << std::endl;

  // Contract the surface mesh as much as possible
  SMS::Count_stop_predicate<Surface_mesh> stop(0);

  std::cout << "Collapsing as many non-sharp edges of mesh: " << filename << " as possible..." << std::endl;
  int r = SMS::edge_collapse(surface_mesh, stop,
                             CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, surface_mesh))
                                              .halfedge_index_map(get(CGAL::halfedge_external_index, surface_mesh))
                                              .edge_is_constrained_map(constraints_map)
                                              .get_placement(placement));

  std::cout << "\nFinished!\n" << r << " edges removed.\n"
            << num_edges(surface_mesh) << " final edges.\n";
  std::ofstream os(argc > 2 ? argv[2] : "out.off"); os << surface_mesh;

  std::cout  << "Checking sharped edges were preserved...\n";
  // check sharp edges were preserved
  for(edge_descriptor ed : edges(surface_mesh))
  {
    halfedge_descriptor hd = halfedge(ed, surface_mesh);
    if(is_border(ed, surface_mesh))
    {
      --nb_sharp_edges;
      assert(constrained_edges[ed] == std::make_pair(point(source(hd, surface_mesh), surface_mesh),
                                                     point(target(hd, surface_mesh), surface_mesh)));
    }
    else
    {
      double angle = approximate_dihedral_angle(point(target(opposite(hd, surface_mesh), surface_mesh), surface_mesh),
                                                point(target(hd, surface_mesh), surface_mesh),
                                                point(target(next(hd, surface_mesh), surface_mesh), surface_mesh),
                                                point(target(next(opposite(hd, surface_mesh), surface_mesh), surface_mesh), surface_mesh));
      if(CGAL::abs(angle) < 100)
      {
        --nb_sharp_edges;
        assert(constrained_edges[ed] == std::make_pair(point(source(hd, surface_mesh), surface_mesh),
                                                       point(target(hd, surface_mesh), surface_mesh)));
      }
    }
  }

  std::cout  << "OK\n";

  std::cout << "Check that no removable edge has been forgotten..." << std::endl;
  r = SMS::edge_collapse(surface_mesh,
                         stop,
                         CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, surface_mesh))
                                          .halfedge_index_map(get(CGAL::halfedge_external_index, surface_mesh))
                                          .edge_is_constrained_map(constraints_map)
                                          .get_placement(placement)
  );

  assert(r == 0);

  if(r == 0) {
    std::cout << "OK\n";
  } else {
    std::cout << "ERROR! " << r << " edges removed!\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
