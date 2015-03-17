#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <boost/property_map/property_map.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point;
typedef Kernel::Compare_dihedral_angle_3                    Compare_dihedral_angle_3;

typedef CGAL::Surface_mesh<Point>                           Mesh;


template <typename G>
struct Constraint : public boost::put_get_helper<bool,Constraint<G> >
{
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef boost::readable_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;

  Constraint()
    :g_(NULL)
  {}

  Constraint(G & g, double bound) 
    : g_(&g), bound_(bound)
  {}

  bool operator[](edge_descriptor e) const {
    return compare_((*g_).point(source(e,*g_)),
                   (*g_).point(target(e,*g_)),
                   (*g_).point(target(next(halfedge(e,*g_),*g_),*g_)),
                   (*g_).point(target(next(opposite(halfedge(e,*g_),*g_),*g_),*g_)),
                   bound_) == CGAL::SMALLER;
  }

  G* g_;
  Compare_dihedral_angle_3 compare_;
  double bound_;
};


int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/eight.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  const double bound = std::cos(0.7* CGAL_PI);

  std::vector<face_descriptor> cc;
  face_descriptor fd = *faces(mesh).first;
  CGAL::Polygon_mesh_processing::connected_component(fd,
                                        mesh,
                                        std::back_inserter(cc));

  std::cerr << "connected components without edge constraints" << std::endl;
  std::cerr << cc.size() << " faces in the CC of " << fd << std::endl;


  std::cerr << "\nconnected components with edge constraints (dihedral angle < 3/4 pi)" << std::endl;
  Mesh::Property_map<face_descriptor, std::size_t> fccmap;
  fccmap = mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
  std::size_t num = CGAL::Polygon_mesh_processing::connected_components(mesh,
                                        fccmap,
                                        Constraint<Mesh>(mesh, bound));
  
  std::cerr << "The graph has " << num << " connected components (face connectivity)" << std::endl;
  BOOST_FOREACH(face_descriptor f , faces(mesh)){
    std::cout  << f << " in connected component " << fccmap[f] << std::endl;
  }
 
  std::cerr << "We keep the two largest components" << std::endl; 
  CGAL::Polygon_mesh_processing::keep_largest_connected_components(mesh,
                                        2,
                                        Constraint<Mesh>(mesh, bound));

  return 0;
}
