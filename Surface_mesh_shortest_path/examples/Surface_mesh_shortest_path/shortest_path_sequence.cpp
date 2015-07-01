#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Random.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/variant.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron_3> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef Traits::Barycentric_coordinate Barycentric_coordinate;
typedef boost::graph_traits<Polyhedron_3> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;
typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;

// A model of SurfaceMeshShortestPathVisitor storing simplicies
// using boost::variant
struct Sequence_collector
{
  typedef boost::variant< vertex_descriptor,
                         std::pair<halfedge_descriptor,double>,
                         std::pair<face_descriptor, Barycentric_coordinate> > Simplex;
  std::vector< Simplex > sequence;

  void operator()(halfedge_descriptor he, double alpha)
  {

    sequence.push_back( std::make_pair(he, alpha) );
  }

  void operator()(vertex_descriptor v)
  {
    sequence.push_back( v );
  }

  void operator()(face_descriptor f, Barycentric_coordinate alpha)
  {
    sequence.push_back( std::make_pair(f, alpha) );
  }
};

// A visitor to print what a variant contains using boost::apply_visitor
struct Print_visitor : public boost::static_visitor<> {
  int i;
  Polyhedron_3& g;

  Print_visitor(Polyhedron_3& g) :i(-1), g(g) {}

  void operator()(vertex_descriptor v)
  {
    std::cout << "#" << ++i << " : Vertex : " << get(boost::vertex_index, g)[v] << "\n";
  }

  void operator()(const std::pair<halfedge_descriptor,double>& h_a)
  {
    std::cout << "#" << ++i << " : Edge : " << get(CGAL::halfedge_index, g)[h_a.first] << " , ("
                                            << 1.0 - h_a.second << " , "
                                            << h_a.second << ")\n";
  }

  void operator()(const std::pair<face_descriptor, Barycentric_coordinate>& f_bc)
  {
    std::cout << "#" << ++i << " : Face : " << get(CGAL::face_index, g)[f_bc.first] << " , ("
                                            << f_bc.second[0] << " , "
                                            << f_bc.second[1] << " , "
                                            << f_bc.second[2] << ")\n";
  }
};

int main(int argc, char** argv)
{
  // read input polyhedron
  Polyhedron_3 polyhedron;
  std::ifstream input((argc>1)?argv[1]:"data/elephant.off");
  input >> polyhedron;
  input.close();

  // initialize indices of vertices, halfedges and facets
  CGAL::set_halfedgeds_items_id(polyhedron);

  // pick up a random face
  const size_t randSeed = argc > 2 ? std::atoi(argv[2]) : 7915421;
  CGAL::Random rand(randSeed);
  const int target_face_index = rand.get_int(0, num_faces(polyhedron));
  face_iterator face_it = faces(polyhedron).first;
  std::advance(face_it,target_face_index);
  // ... and define a barycentric coordinate inside the face
  Barycentric_coordinate face_location = {{0.25, 0.5, 0.25}};

  // construct a shortest path query object and add a source point
  Surface_mesh_shortest_path shortest_paths(polyhedron);
  shortest_paths.add_source_point(*face_it, face_location);

  // pick a random target point inside a face
  face_it = faces(polyhedron).first;
  std::advance(face_it, rand.get_int(0, num_faces(polyhedron)));

  // collect the sequence of simplicies crossed by the shortest path
  Sequence_collector sequence_collector;
  shortest_paths.shortest_path_sequence_to_source_points(*face_it, face_location, sequence_collector);

  // print the sequence using the visitor pattern
  Print_visitor print_visitor(polyhedron);
  for (size_t i = 0; i < sequence_collector.sequence.size(); ++i)
    boost::apply_visitor(print_visitor, sequence_collector.sequence[i]);

  return 0;
}
