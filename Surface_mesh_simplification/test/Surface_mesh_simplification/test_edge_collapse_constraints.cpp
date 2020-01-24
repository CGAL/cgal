// Copyright (c) 2010  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/property_map.h>

#include <cmath>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                              Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor      edge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_iterator        edge_iterator;

namespace SMS = CGAL::Surface_mesh_simplification;

// BGL property map which indicates whether an edge is marked as non-removable
struct Constrained_edge_map
  : public boost::put_get_helper<bool, Constrained_edge_map>
{
  typedef boost::readable_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;

  Constrained_edge_map(const CGAL::Unique_hash_map<key_type,bool>& aConstraints)
    : mConstraints(aConstraints)
  {}

  reference operator[](const key_type& e) const { return  is_constrained(e); }

  bool is_constrained(const key_type& e) const { return mConstraints.is_defined(e); }

private:
  const CGAL::Unique_hash_map<key_type,bool>& mConstraints;
};

bool is_border (edge_descriptor e, const Surface_mesh& sm)
{
  return (face(halfedge(e,sm),sm) == boost::graph_traits<Surface_mesh>::null_face()) ||
         (face(opposite(halfedge(e,sm),sm),sm) == boost::graph_traits<Surface_mesh>::null_face());
}

Point_3 point(vertex_descriptor vd,  const Surface_mesh& sm)
{
  return get(CGAL::vertex_point, sm, vd);
}

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const char* filename = (argc > 1) ? argv[1] : "data/arrow.off";
  std::ifstream is(filename);
  if(!is || !(is >> surface_mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Unique_hash_map<edge_descriptor, bool> constraint_hmap(false);
  Constrained_edge_map constraints_map(constraint_hmap);
  SMS::Midpoint_placement<Surface_mesh> placement;
  std::size_t nb_sharp_edges = 0;

  // detect border edges
  for(edge_descriptor ed : edges(surface_mesh))
  {
    halfedge_descriptor hd = halfedge(ed, surface_mesh);
    if(is_border(ed, surface_mesh))
    {
      ++nb_sharp_edges;
      constraint_hmap[ed] = true;
    }
  }

  std::cerr << "# border edges = " << nb_sharp_edges << std::endl;

  // Contract the surface mesh as much as possible
  SMS::Count_stop_predicate<Surface_mesh> stop(0);
  
  SMS::edge_collapse(surface_mesh, stop,
                     CGAL::parameters::edge_is_constrained_map(constraints_map)
                     .max_normal_angle_change(5*CGAL_PI/180.0)
                     .get_placement(placement));
  assert(surface_mesh.number_of_edges() == 15 );
  
  SMS::edge_collapse(surface_mesh, stop,
                     CGAL::parameters::edge_is_constrained_map(constraints_map)
                     .max_normal_angle_change(45*CGAL_PI/180.0)
                     .get_placement(placement));
  
  assert(surface_mesh.number_of_edges() == 12 );
  
  
  SMS::edge_collapse(surface_mesh, stop,
                     CGAL::parameters::edge_is_constrained_map(constraints_map)
                     .max_normal_angle_change(CGAL_PI)
                     .get_placement(placement));
  std::ofstream out("out.off");
  out << surface_mesh;
  out.close();
assert(surface_mesh.number_of_edges() == 9 );
  return EXIT_SUCCESS;
}
