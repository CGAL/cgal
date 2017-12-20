#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/helpers.h>
#include <iostream>
#include <fstream>
#include <cstring>

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


int main(int argc, char* argv[])
{

  typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  const double bound = std::cos(0.7* CGAL_PI);
  const char* filename = (argc > 1) ? argv[1] : "data/blobby_3cc.off";
  Mesh sm;
  std::ifstream in(filename);
  in >> sm;

  std::cout << "VEF " << num_vertices(sm) << " " << num_edges(sm) << " " << num_faces(sm) << "\n";

  std::vector<face_descriptor> cc;
  face_descriptor fd = *faces(sm).first;
  CGAL::Polygon_mesh_processing::connected_component(fd,
                                                     sm,
                                                     std::back_inserter(cc));

  std::cerr << "connected components without edge constraints" << std::endl;
  std::cerr << cc.size() << " faces in the CC of " << fd << std::endl;
  if (strcmp(filename, "data/blobby_3cc.off") == 0)
    assert(cc.size() == 1452);

  std::cerr << "\nconnected components with edge constraints (dihedral angle < 3/4 pi)" << std::endl;
  Mesh::Property_map<face_descriptor,std::size_t> fccmap;
  fccmap = sm.add_property_map<face_descriptor,std::size_t>("f:CC").first;
  std::size_t num =
    PMP::connected_components(sm,
      fccmap
    );

 std::cerr << "The graph has " << num << " connected components (face connectivity)" << std::endl;
 if (strcmp(filename, "data/blobby_3cc.off") == 0)
   assert(num == 3);

 std::vector<face_descriptor> one_face_per_cc(num);
 std::vector<std::size_t> cc_size(num,0);


 BOOST_FOREACH(face_descriptor f , faces(sm)){
 //  std::cout  << f << " in connected component " << fccmap[f] << std::endl;
    std::size_t ccid=fccmap[f];
    if (++cc_size[ccid]==1)
      one_face_per_cc[ccid]=f;
 }

  std::size_t id_of_cc_to_remove =
    std::distance(cc_size.begin(),
                  std::min_element(cc_size.begin(), cc_size.end()));

  Mesh copy1 = sm;
  Mesh copy2 = sm;

  // remove cc from copy1
  std::vector<face_descriptor> ff;
  for (std::size_t i=0;i<num;++i)
    if (i!=id_of_cc_to_remove)
      ff.push_back(one_face_per_cc[i]);
  PMP::keep_connected_components(copy1, ff,
     PMP::parameters::edge_is_constrained_map(Constraint<Mesh>(copy1, bound)));

  // remove cc from copy2
  ff.clear();
  ff.push_back(one_face_per_cc[id_of_cc_to_remove]);
  PMP::remove_connected_components(copy2, ff,
     PMP::parameters::edge_is_constrained_map(Constraint<Mesh>(copy2, bound)));

  std::cerr << "We keep the " << num-1 << " largest components" << std::endl;
  PMP::keep_largest_connected_components(sm, num-1,
    PMP::parameters::edge_is_constrained_map(Constraint<Mesh>(sm, bound)));

  sm.collect_garbage();
  copy1.collect_garbage();
  copy2.collect_garbage();

  assert( num_vertices(sm)==num_vertices(copy1) && num_vertices(copy1)==num_vertices(copy2) );
  assert( num_edges(sm)==num_edges(copy1) && num_edges(copy1)==num_edges(copy2) );
  assert( num_faces(sm)==num_faces(copy1) && num_faces(copy1)==num_faces(copy2) );

  {
    Mesh m;
    Point p(0,0,0), q(1,0,0), r(0,1,0), s(0,0,1);
    CGAL::make_tetrahedron(p,q,r,s,m);
    CGAL::make_triangle(p,q,r,m);
    CGAL::make_tetrahedron(p,q,r,s,m);
    PMP::keep_large_connected_components(m, 4);
    m.collect_garbage();
    assert( num_vertices(m) == 8);
  }

  return 0;
}
