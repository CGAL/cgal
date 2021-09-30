#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/planar_segmentation.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;
int main(int argc, char** argv)
{
  Surface_mesh sm;
  // When using input with double coordinates, better use minimum_cosinus_squared
  // parameter set to 0.999 because the default value 1 may not succeed for RG because
  // it will try to fit regions exactly!
  std::ifstream in(argc > 1 ? argv[1] : "data/cube_quad.off");
  in >> sm;

  // triangulate faces;
  // PMP::triangulate_faces(sm);
  // assert(faces(sm).size()==12);

  // Surface_mesh::Property_map<Surface_mesh::Edge_index, bool> ecm =
  //   sm.add_property_map<Surface_mesh::Edge_index, bool>("ecm",false).first;

  // // detect sharp edges of the cube
  // PMP::detect_sharp_edges(sm, 60, ecm);

  // // create a remeshed version of the cube with many elements
  // PMP::isotropic_remeshing(faces(sm), 0.1, sm, CGAL::parameters::edge_is_constrained_map(ecm));
  // std::ofstream("cube_remeshed.off") << sm;
  // assert(faces(sm).size()>100);

  // decimate the mesh
  // PMP::decimate(sm); // use region growing approach by default
  // PMP::decimate(sm, CGAL::parameters::use_region_growing(false)); // use connected components approach
  // PMP::decimate(sm, CGAL::parameters::use_region_growing(true)); // use region growing approach
  // PMP::decimate(sm, CGAL::parameters::use_region_growing(false).maximum_Frechet_distance(30.0)); // use PCA approach

  std::cout << std::endl;
  // PMP::decimate(sm, CGAL::parameters::use_region_growing(false));
  const bool success = PMP::decimate(sm); // this is a stable version with approximate parameters
  std::ofstream("cube_decimated.off") << sm;
  std::cout << " SUCCESS: " << success << std::endl;
  std::cout << std::endl;

  // we should be back to 12 faces
  // assert(faces(sm).size()==12);

  return 0;
}
