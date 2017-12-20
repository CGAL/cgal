#include <cstdlib>
#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Random.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_shortest_path.h>

#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef Traits::Barycentric_coordinates Barycentric_coordinates;
typedef boost::graph_traits<Triangle_mesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;
typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;

// A model of SurfacemeshShortestPathVisitor storing simplicies
// using boost::variant
struct Sequence_collector
{
  typedef boost::variant< vertex_descriptor,
                         std::pair<halfedge_descriptor,double>,
                         std::pair<face_descriptor, Barycentric_coordinates> > Simplex;
  std::vector< Simplex > sequence;

  void operator()(halfedge_descriptor he, double alpha)
  {

    sequence.push_back( std::make_pair(he, alpha) );
  }

  void operator()(vertex_descriptor v)
  {
    sequence.push_back( v );
  }

  void operator()(face_descriptor f, Barycentric_coordinates alpha)
  {
    sequence.push_back( std::make_pair(f, alpha) );
  }
};

// A visitor to print what a variant contains using boost::apply_visitor
struct Print_visitor : public boost::static_visitor<> {
  int i;
  Triangle_mesh& g;

  Print_visitor(Triangle_mesh& g) :i(-1), g(g) {}

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

  void operator()(const std::pair<face_descriptor, Barycentric_coordinates>& f_bc)
  {
    std::cout << "#" << ++i << " : Face : " << get(CGAL::face_index, g)[f_bc.first] << " , ("
                                            << f_bc.second[0] << " , "
                                            << f_bc.second[1] << " , "
                                            << f_bc.second[2] << ")\n";
  }
};

int main(int argc, char** argv)
{
  Triangle_mesh tmesh;
  std::ifstream input((argc>1)?argv[1]:"data/elephant.off");
  input >> tmesh;
  input.close();

  // pick up a random face
  const unsigned int randSeed = argc > 2 ? boost::lexical_cast<unsigned int>(argv[2]) : 7915421;
  CGAL::Random rand(randSeed);
  const int target_face_index = rand.get_int(0, static_cast<int>(num_faces(tmesh)));
  face_iterator face_it = faces(tmesh).first;
  std::advance(face_it,target_face_index);
  // ... and define a barycentric coordinates inside the face
  Barycentric_coordinates face_location = {{0.25, 0.5, 0.25}};

  // construct a shortest path query object and add a source point
  Surface_mesh_shortest_path shortest_paths(tmesh);
  shortest_paths.add_source_point(*face_it, face_location);

  // pick a random target point inside a face
  face_it = faces(tmesh).first;
  std::advance(face_it, rand.get_int(0, static_cast<int>(num_faces(tmesh))));

  // collect the sequence of simplicies crossed by the shortest path
  Sequence_collector sequence_collector;
  shortest_paths.shortest_path_sequence_to_source_points(*face_it, face_location, sequence_collector);

  // print the sequence using the visitor pattern
  Print_visitor print_visitor(tmesh);
  for (size_t i = 0; i < sequence_collector.sequence.size(); ++i)
    boost::apply_visitor(print_visitor, sequence_collector.sequence[i]);

  return 0;
}
