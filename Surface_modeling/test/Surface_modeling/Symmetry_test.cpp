#include "Surface_modeling_test_commons.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Deform_mesh.h>
#include <CGAL/IO/Polyhedron_iostream.h>
  
typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>  Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_descriptor      edge_descriptor;

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; // use id field of vertices
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map;   // use id field of edges

typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map, CGAL::ORIGINAL_ARAP> Deform_mesh_arap;
typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map, CGAL::SPOKES_AND_RIMS> Deform_mesh_spoke;

int main()
{
  Polyhedron mesh_1;
  if(!read_to_polyhedron("data/square.off", mesh_1)) { return EXIT_FAILURE; }
  Polyhedron mesh_2 = mesh_1;

  Deform_mesh_arap deform_mesh_arap(mesh_1, Vertex_index_map(), Edge_index_map()); 
  Deform_mesh_spoke deform_mesh_spoke(mesh_2, Vertex_index_map(), Edge_index_map()); 

  const int deformation_iteration = 500;
  const double x = -0.45; const double y = -0.65; const double z = -0.0;

  std::cerr << "ORIGINAL_ARAP performance: " << std::endl;
  bool successed = preprocess_and_deform(deform_mesh_arap,
    "data/Symmetry_test_roi.txt",
    "data/Symmetry_test_handle.txt",
    Deform_mesh_arap::Vector(x, y, z),
    deformation_iteration);
  if(!successed) { return EXIT_FAILURE; }

  std::cerr << "SPOKES_AND_RIMS performance: " << std::endl;
  successed = preprocess_and_deform(deform_mesh_spoke,
    "data/Symmetry_test_roi.txt",
    "data/Symmetry_test_handle.txt",
    Deform_mesh_spoke::Vector(x, y, z),
    deformation_iteration);
  if(!successed) { return EXIT_FAILURE; }

  std::cerr << "Save deformed models" << std::endl;
  std::ofstream("data/Symmetry_test_deformed_arap.off") << mesh_1;
  std::ofstream("data/Symmetry_test_deformed_spokes.off") << mesh_2;
  std::cerr << "All done!" << std::endl;
  return EXIT_SUCCESS;
}

