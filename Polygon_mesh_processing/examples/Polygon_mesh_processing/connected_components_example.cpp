#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef Kernel::Point_3                    Point;
typedef Kernel::Compare_dihedral_angle_3   Compare_dihedral_angle_3;
typedef CGAL::Surface_mesh<Point>          Mesh;


template <typename G>
struct Constraint : public boost::put_get_helper<bool,Constraint<G> >{
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef boost::readable_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;

  Constraint()
    :g(NULL)
  {}

  Constraint(G & g, double bound) 
    : g(&g), bound(bound)
  {}

  bool operator[](edge_descriptor e) const {
    return compare((*g).point(source(e,*g)),
                   (*g).point(target(e,*g)),
                   (*g).point(target(next(halfedge(e,*g),*g),*g)),
                   (*g).point(target(next(opposite(halfedge(e,*g),*g),*g),*g)),
                   bound) == CGAL::SMALLER;
  }

  G* g;
  Compare_dihedral_angle_3 compare;
  double bound;
};


int main(int, char* argv[]) 
{
  typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  const double bound = std::cos(0.7* CGAL_PI);
  Mesh sm;
  std::ifstream in(argv[1]);
  in >> sm;
  
  std::vector<face_descriptor> cc;
  face_descriptor fd = *faces(sm).first;
  CGAL::Polygon_mesh_processing::connected_component(fd,
                                                     sm,
                                                     std::back_inserter(cc));

  std::cerr << "connected components without edge constraints" << std::endl;
  std::cerr << cc.size() << " faces in the CC of " << fd << std::endl;


  std::cerr << "\nconnected components with edge constraints (dihedral angle < 3/4 pi)" << std::endl;
  Mesh::Property_map<face_descriptor,std::size_t> fccmap;
  fccmap = sm.add_property_map<face_descriptor,std::size_t>("f:CC").first; 
  std::size_t num = PMP::connected_components(sm,
                                              fccmap,
                                              Constraint<Mesh>(sm,bound)
                                              );
  
 std::cerr << "The graph has " << num << " connected components (face connectivity)" << std::endl;
 BOOST_FOREACH(face_descriptor f , faces(sm)){
   std::cout  << f << " in connected component " << fccmap[f] << std::endl;
  }
 
 std::cerr << "We keep the two largest components" << std::endl; 
 PMP::keep_largest_connected_components(sm,2,Constraint<Mesh>(sm,bound));

 std::cout << "mesh:\n" << sm << std::endl;
  return 0;
}
