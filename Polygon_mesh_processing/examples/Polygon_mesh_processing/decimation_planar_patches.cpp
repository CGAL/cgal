#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <iostream>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;
int main()
{
  Surface_mesh sm;
  std::ifstream in("data/cube_quad.off");
  in >> sm;

  // triangulate faces;
  PMP::triangulate_faces(sm);
  assert(faces(sm).size()==12);

  Surface_mesh::Property_map<Surface_mesh::Edge_index, bool> ecm =
    sm.add_property_map<Surface_mesh::Edge_index, bool>("ecm",false).first;

  // detect sharp edges of the cube
  PMP::detect_sharp_edges(sm, 60, ecm);

  // create a remeshed version of the cube with many elements
  PMP::isotropic_remeshing(faces(sm), 0.1, sm, CGAL::parameters::edge_is_constrained_map(ecm));
  std::ofstream("cube_remeshed.off") << sm;
  assert(faces(sm).size()>100);
  
  // decimate the mesh
  PMP::remesh_planar_patches(sm);
  std::ofstream("cube_decimated.off") << sm;

  // we should be back to 12 faces
  assert(faces(sm).size()==12);

  return 0;
}
