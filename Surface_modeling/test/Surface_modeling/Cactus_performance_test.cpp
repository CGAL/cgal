#include "Surface_modeling_test_commons.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Timer.h>
#include <CGAL/Deform_mesh.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>  Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_descriptor      edge_descriptor;

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; // use id field of vertices
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map;   // use id field of edges

// #define CGAL_USE_EXPERIMENTAL_POLAR
#ifdef CGAL_USE_EXPERIMENTAL_POLAR
#include <CGAL/internal/Surface_modeling/Deformation_Eigen_polar_closest_rotation_traits_3.h>
typedef CGAL::Deformation_Eigen_polar_closest_rotation_traits_3 Closest_rotation_model;
#else
typedef CGAL::Deformation_Eigen_closest_rotation_traits_3 Closest_rotation_model;
#endif

typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map, 
  CGAL::ORIGINAL_ARAP, CGAL::Default, CGAL::Default, Closest_rotation_model> Deform_mesh_arap;

typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map, 
  CGAL::SPOKES_AND_RIMS, CGAL::Default, CGAL::Default, Closest_rotation_model> Deform_mesh_spoke;


template<class DeformMesh>
bool preprocess_and_deform(DeformMesh& deform_mesh, int deformation_iteration) 
{
  boost::optional<typename DeformMesh::Handle_group> active_handle_group = 
    read_rois(deform_mesh, "data/cactus_roi.txt", "data/cactus_handle.txt");
  if(!active_handle_group) { return false; }

  CGAL::Timer timer; timer.start();

  if(!deform_mesh.preprocess()) {
    std::cerr << "Error: preprocess() failed!" << std::endl;
    return false;
  }
  std::cerr << "Preprocess time: " << timer.time() << std::endl;
  timer.reset();

  deform_mesh.translate(*active_handle_group, typename DeformMesh::Vector(-0.55, -0.30, 0.0) );
  deform_mesh.deform(deformation_iteration, 0);
  timer.stop();

  std::cerr << "Deformation time (one handle translated and deform method iterated " <<
    deformation_iteration << " times): " << timer.time() << std::endl;
  return true;
}

int main()
{
  Polyhedron mesh_1;
  if(!read_to_polyhedron("data/cactus.off", mesh_1)) { return EXIT_FAILURE; }
  Polyhedron mesh_2 = mesh_1;
  
  Deform_mesh_arap deform_mesh_arap(mesh_1, Vertex_index_map(), Edge_index_map()); 
  Deform_mesh_spoke deform_mesh_spoke(mesh_2, Vertex_index_map(), Edge_index_map()); 

  const int deformation_iteration = 500;
  const double x = -0.55; const double y = -0.50; const double z = -0.0;

  std::cerr << "ORIGINAL_ARAP performance: " << std::endl;
  bool successed = preprocess_and_deform(deform_mesh_arap,
    "data/cactus_roi.txt",
    "data/cactus_handle.txt",
    Deform_mesh_arap::Vector(x, y, z),
    deformation_iteration);
  if(!successed) { return EXIT_FAILURE; }

  std::cerr << "SPOKES_AND_RIMS performance: " << std::endl;
  successed = preprocess_and_deform(deform_mesh_spoke,
    "data/cactus_roi.txt",
    "data/cactus_handle.txt",
    Deform_mesh_spoke::Vector(x, y, z),
    deformation_iteration);
  if(!successed) { return EXIT_FAILURE; }

  std::cerr << "Save deformed models" << std::endl;
  std::ofstream("data/cactus_deformed_arap.off") << mesh_1;
  std::ofstream("data/cactus_deformed_spokes.off") << mesh_2;
  std::cerr << "All done!" << std::endl;
  return EXIT_SUCCESS;
}

