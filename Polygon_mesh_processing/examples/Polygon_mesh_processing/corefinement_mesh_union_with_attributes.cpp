#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <boost/container/flat_map.hpp>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Surface_mesh<K::Point_3>                       Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

struct Visitor :
  public PMP::Corefinement::Default_visitor<Mesh>
{
  typedef Mesh::Face_index face_descriptor;

  boost::container::flat_map<const Mesh*, Mesh::Property_map<Mesh::Face_index, int> > properties;
  int face_id;

  Visitor()
  {
    properties.reserve(3);
    face_id=-1;
  }

// visitor API overloaded
  void before_subface_creations(face_descriptor f_split,Mesh& tm)
  {
    face_id = properties[&tm][f_split];
  }

  void after_subface_created(face_descriptor f_new,Mesh& tm)
  {
    properties[&tm][f_new] = face_id;
  }

  void after_face_copy(face_descriptor f_src, Mesh& tm_src,
                       face_descriptor f_tgt, Mesh& tm_tgt)
  {
    properties[&tm_tgt][f_tgt] = properties[&tm_src][f_src];
  }
};

int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/blobby.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/eight.off";

  Mesh mesh1, mesh2;
  if(!PMP::IO::read_polygon_mesh(filename1, mesh1) || !PMP::IO::read_polygon_mesh(filename2, mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  Mesh out;

  // add a property in each mesh to track the parent mesh for each face of the output
  Mesh::Property_map<Mesh::Face_index, int>
    mesh1_id = mesh1.add_property_map<Mesh::Face_index, int>("f:id", -1).first,
    mesh2_id = mesh2.add_property_map<Mesh::Face_index, int>("f:id", -1).first,
    out_id = out.add_property_map<Mesh::Face_index, int>("f:id", -1).first;

  // init the face ids (for the purpose of the example but choosing 1 (2) as default value of the map would avoid the loop)
  for(Mesh::Face_index f : faces(mesh1))
    mesh1_id[f] = 1;
  for(Mesh::Face_index f : faces(mesh2))
    mesh2_id[f] = 2;

  Visitor visitor;
  visitor.properties[&mesh1] = mesh1_id;
  visitor.properties[&mesh2] = mesh2_id;
  visitor.properties[&out] = out_id;

  bool valid_union = PMP::corefine_and_compute_union(mesh1, mesh2, out, PMP::parameters::visitor(visitor));

  for(Mesh::Face_index f : faces(mesh1))
    assert( mesh1_id[f] == 1 );
  for(Mesh::Face_index f : faces(mesh2))
    assert( mesh2_id[f] == 2 );
  for(Mesh::Face_index f : faces(out))
    assert( out_id[f]==1 || out_id[f]==2);

  if (valid_union)
  {
    CGAL::IO::write_polygon_mesh("union.off", out, CGAL::parameters::stream_precision(17));
    return 0;
  }

  std::cout << "Union could not be computed\n";
  return 1;
}
