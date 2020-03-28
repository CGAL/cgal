#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/selection.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;


struct Vector_pmap_wrapper{
  std::vector<bool>& vect;
  Vector_pmap_wrapper(std::vector<bool>& v) : vect(v) {}
  friend bool get(const Vector_pmap_wrapper& m, face_descriptor f)
  {
    return m.vect[f];
  }
  friend void put(const Vector_pmap_wrapper& m, face_descriptor f, bool b)
  {
    m.vect[f]=b;
  }
};

int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/blobby.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/eight.off";
  std::ifstream input(filename1);

  Mesh mesh1, mesh2;
  if (!input || !(input >> mesh1))
  {
    std::cerr << "First mesh is not a valid off file." << std::endl;
    return 1;
  }
  input.close();
  input.open(filename2);
  if (!input || !(input >> mesh2))
  {
    std::cerr << "Second mesh is not a valid off file." << std::endl;
    return 1;
  }

  //create a property on edges to indicate whether they are constrained
  Mesh::Property_map<edge_descriptor,bool> is_constrained_map =
    mesh1.add_property_map<edge_descriptor,bool>("e:is_constrained", false).first;

  // update mesh1 to contain the mesh bounding the difference
  // of the two input volumes.
  bool valid_difference =
    PMP::corefine_and_compute_difference(mesh1,
                                         mesh2,
                                         mesh1,
                                         params::all_default(), // default parameters for mesh1
                                         params::all_default(), // default parameters for mesh2
                                         params::edge_is_constrained_map(is_constrained_map));

  if (valid_difference)
  {
    std::cout << "Difference was successfully computed\n";
    std::ofstream output("difference.off");
    output.precision(17);
    output << mesh1;
  }
  else{
    std::cout << "Difference could not be computed\n";
    return 1;
  }

  // collect faces incident to a constrained edge
  std::vector<face_descriptor> selected_faces;
  std::vector<bool> is_selected(num_faces(mesh1), false);
  for(edge_descriptor e : edges(mesh1))
    if (is_constrained_map[e])
    {
      // insert all faces incident to the target vertex
      for(halfedge_descriptor h :
                    halfedges_around_target(halfedge(e,mesh1),mesh1))
      {
        if (!is_border(h, mesh1) )
        {
          face_descriptor f=face(h, mesh1);
          if ( !is_selected[f] )
          {
            selected_faces.push_back(f);
            is_selected[f]=true;
          }
        }
      }
    }

  // increase the face selection
  CGAL::expand_face_selection(selected_faces, mesh1, 2,
    Vector_pmap_wrapper(is_selected), std::back_inserter(selected_faces));

  std::cout << selected_faces.size()
            << " faces were selected for the remeshing step\n";

  // remesh the region around the intersection polylines
  PMP::isotropic_remeshing(
    selected_faces,
    0.02,
    mesh1,
    params::edge_is_constrained_map(is_constrained_map) );

  std::ofstream output("difference_remeshed.off");
  output.precision(17);
  output << mesh1;

  return 0;
}
