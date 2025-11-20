// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_LINEAR_CELL_COMPLEX_OCTREE_GENERATION_H
#define CGAL_LINEAR_CELL_COMPLEX_OCTREE_GENERATION_H 1

//#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
//#include <CGAL/IO/OFF.h>
#include <CGAL/assertions.h>
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Polyhedron_3.h>
//#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <type_traits>
#include <string>
#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

namespace CGAL {

 /** @file Linear_cell_complex_octree_generation.h
  * Octree generation for linear cell complexes:
  * - Preferred overload: from an in-memory FaceGraph (e.g. CGAL::Surface_mesh)
  * - Convenience overloads: from an OFF filename (std::string / const char*)
  * A boolean `regularized` (currently a stub) is reserved for future 1‑irregular balancing
  */

namespace internal {

// Internal AABB intersector for octree generation
/*
template<typename Kernel>
class Simple_AABB_intersector {
private:
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
  typedef AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
  typedef AABB_traits_3<Kernel, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  Tree tree;
  Polyhedron polyhedron;
  bool valid;

public:
  Simple_AABB_intersector() : valid(false) {}

  explicit Simple_AABB_intersector(const std::string& off_filename) : valid(false) {
    std::ifstream off_file(off_filename);
    if (!off_file.good()) {
      std::cerr << "Error: cannot open " << off_filename << std::endl;
      return;
    }

    off_file >> polyhedron;
    CGAL::Polygon_mesh_processing::triangulate_faces(polyhedron);

    // Compute AABB tree
    tree.insert(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.accelerate_distance_queries();

    if (!tree.empty()) {
      valid = true;
    }
  }

  bool empty() const { return !valid || tree.empty(); }
  typename Tree::Bounding_box bbox() const { return tree.bbox(); }

bool is_outside(double x1, double y1, double z1,
               double x2, double y2, double z2) const {
  if (!valid) return true;

// Create the cube
  typename Kernel::Iso_cuboid_3 cube(
    typename Kernel::Point_3(x1, y1, z1),
    typename Kernel::Point_3(x2, y2, z2));

  // Direct intersection test with AABB tree
  return !tree.do_intersect(cube);
}

bool is_intersect(double x1, double y1, double z1,
                 double x2, double y2, double z2) const {
  if (!valid) return false;

  typename Kernel::Iso_cuboid_3 cube(
    typename Kernel::Point_3(x1, y1, z1),
    typename Kernel::Point_3(x2, y2, z2));

  return tree.do_intersect(cube);
}
};
*/

template<typename Kernel, typename FaceGraph>
class Simple_AABB_intersector_fgraph {
  typedef AABB_face_graph_triangle_primitive<FaceGraph> Primitive;
  typedef AABB_traits_3<Kernel, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  const FaceGraph* fg;
  Tree tree;
  bool valid;

public:
  Simple_AABB_intersector_fgraph(const FaceGraph& graph)
    : fg(&graph), valid(false)
  {
    if(faces(graph).first == faces(graph).second) return;
    tree.insert(faces(graph).first, faces(graph).second, graph);
    tree.accelerate_distance_queries();
    valid = !tree.empty();
  }

  bool empty() const { return !valid; }
  typename Tree::Bounding_box bbox() const { return tree.bbox(); }

  bool is_outside(double x1,double y1,double z1,
                  double x2,double y2,double z2) const
  {
    if(!valid) return true;
    typename Kernel::Iso_cuboid_3 cube(
      typename Kernel::Point_3(x1,y1,z1),
      typename Kernel::Point_3(x2,y2,z2));
    return !tree.do_intersect(cube);
  }

  bool is_intersect(double x1,double y1,double z1,
                    double x2,double y2,double z2) const
  {
    if(!valid) return false;
    typename Kernel::Iso_cuboid_3 cube(
      typename Kernel::Point_3(x1,y1,z1),
      typename Kernel::Point_3(x2,y2,z2));
    return tree.do_intersect(cube);
    }
};

// ==== regularization stub ========================================
template<typename LCC>
void regularize_octree(LCC& /*lcc*/)
{
  // Stub: current octree is uniform (single level) nothing to balance yet.
}
// ===========================================================================

template<typename Intersector>
void compute_initial_grid_size(unsigned int init,
                              const Intersector& intersector,
                              int& longestAxis,
                              double& sx, double& sy, double& sz,
                              unsigned int& initX, unsigned int& initY, unsigned int& initZ)
{
  const auto& bbox = intersector.bbox();
  sx = CGAL::to_double(bbox.xmax() - bbox.xmin());
  sy = CGAL::to_double(bbox.ymax() - bbox.ymin());
  sz = CGAL::to_double(bbox.zmax() - bbox.zmin());

  if (sx >= sy && sx >= sz)     longestAxis = 0;
  else if(sy >= sx && sy >= sz) longestAxis = 1;
  else                          longestAxis = 2;

  if(longestAxis == 0) {
    initX = init;
    initY = std::max(1u, static_cast<unsigned int>(std::ceil(init * (sy / sx))));
    initZ = std::max(1u, static_cast<unsigned int>(std::ceil(init * (sz / sx))));
  } else if(longestAxis == 1) {
    initX = std::max(1u, static_cast<unsigned int>(std::ceil(init * (sx / sy))));
    initY = init;
    initZ = std::max(1u, static_cast<unsigned int>(std::ceil(init * (sz / sy))));
  } else {
    initX = std::max(1u, static_cast<unsigned int>(std::ceil(init * (sx / sz))));
    initY = std::max(1u, static_cast<unsigned int>(std::ceil(init * (sy / sz))));
    initZ = init;
  }

  sx /= initX; sy /= initY; sz /= initZ;
}

template<typename LCC, typename Intersector>
void create_initial_hexahedral_grid(LCC& lcc,
                                   const Intersector& intersector,
                                   double sx, double sy, double sz,
                                   unsigned int initX, unsigned int initY, unsigned int initZ,
                                   bool create_all_voxels = false)
{
  CGAL_precondition(!intersector.empty());

  const auto& bbox = intersector.bbox();
  double startx = CGAL::to_double(bbox.xmin());
  double starty = CGAL::to_double(bbox.ymin());
  double startz = CGAL::to_double(bbox.zmin());

  std::size_t created_count = 0;

  if (create_all_voxels) {
    // Create all voxels in bounding box
    for (unsigned int x = 0; x < initX; ++x) {
      for (unsigned int y = 0; y < initY; ++y) {
        for (unsigned int z = 0; z < initZ; ++z) {
          lcc.make_hexahedron(
            typename LCC::Point(startx + x * sx, starty + y * sy, startz + z * sz),
            typename LCC::Point(startx + (x+1) * sx, starty + y * sy, startz + z * sz),
            typename LCC::Point(startx + (x+1) * sx, starty + (y+1) * sy, startz + z * sz),
            typename LCC::Point(startx + x * sx, starty + (y+1) * sy, startz + z * sz),
            typename LCC::Point(startx + x * sx, starty + (y+1) * sy, startz + (z+1) * sz),
            typename LCC::Point(startx + x * sx, starty + y * sy, startz + (z+1) * sz),
            typename LCC::Point(startx + (x+1) * sx, starty + y * sy, startz + (z+1) * sz),
            typename LCC::Point(startx + (x+1) * sx, starty + (y+1) * sy, startz + (z+1) * sz));
          ++created_count;
        }
      }
    }
  } else {
    // Create only voxels that intersect
    for (unsigned int x = 0; x < initX; ++x) {
      for (unsigned int y = 0; y < initY; ++y) {
        for (unsigned int z = 0; z < initZ; ++z) {
          double x1 = startx + x * sx, y1 = starty + y * sy, z1 = startz + z * sz;
          double x2 = startx + (x+1) * sx, y2 = starty + (y+1) * sy, z2 = startz + (z+1) * sz;

          if (!intersector.is_outside(x1, y1, z1, x2, y2, z2)) {
            lcc.make_hexahedron(
              typename LCC::Point(x1, y1, z1), typename LCC::Point(x2, y1, z1),
              typename LCC::Point(x2, y2, z1), typename LCC::Point(x1, y2, z1),
              typename LCC::Point(x1, y2, z2), typename LCC::Point(x1, y1, z2),
              typename LCC::Point(x2, y1, z2), typename LCC::Point(x2, y2, z2));
            ++created_count;
          }
        }
      }
    }
  }

  // Sew adjacent faces
  lcc.sew3_same_facets();

  if (IO::is_pretty(std::cout)) {
    std::cout << "Initial grid created: " << initX << "x" << initY << "x" << initZ
              << " (" << created_count << " hexahedra generated)" << std::endl;
  }
}

} // namespace internal

// /**
//  * Creates an octree approximation of a 3D surface represented by an OFF file.
//  * The octree is built using hexahedral cells in a Linear Cell Complex.
//  *
//  * @tparam LCC a model of `LinearCellComplex` with dimension >= 3
//  * @param lcc the linear cell complex where the octree will be created
//  * @param off_filename path to the OFF file containing the 3D surface
//  * @param initial_grid_size number of initial subdivisions on the longest axis
//  * @param max_subdivision_level maximum octree subdivision level (currently not used)
//  * @param create_all_voxels if true, creates all voxels in bounding box;
//  *                         if false, creates only intersecting voxels
//  * @param no_remove_outside if true, keeps voxels outside the surface (currently not used)
//  * @param regularized if true, requests a 1-irregular (balanced) octree.
//  *                    Currently a stub (octree is uniform).
//  *
//  * @pre `LCC::dimension >= 3`
//  * @pre `LCC::ambient_dimension == 3`
//  */
// template<typename LCC>
// void compute_octree(LCC& lcc,
//                    const std::string& off_filename,
//                    unsigned int initial_grid_size = 2,
//                    unsigned int max_subdivision_level = 5,
//                    bool create_all_voxels = false,
//                    bool no_remove_outside = false,
//                    bool regularized = false)
// {
//   typedef Simple_cartesian<double> Kernel;
//   static_assert(LCC::dimension >= 3, "LCC dimension must be >= 3");
//   static_assert(LCC::ambient_dimension == 3, "LCC ambient dimension must be 3");
//
//   // 1. Create AABB intersector from OFF file
// internal::Simple_AABB_intersector<Kernel> intersector(off_filename);
//
//   if (intersector.empty()) {
//     std::cerr << "Error: cannot create intersector from " << off_filename << std::endl;
//     return;
//   }
//
//   // 2. Compute initial grid size
//   double sx, sy, sz;
//   int longestAxis;
//   unsigned int initX, initY, initZ;
//   internal::compute_initial_grid_size(initial_grid_size, intersector,
//                                      longestAxis, sx, sy, sz,
//                                      initX, initY, initZ);
//
//   // 3. Create initial hexahedral grid
//   internal::create_initial_hexahedral_grid(lcc, intersector,
//                                           sx, sy, sz, initX, initY, initZ,
//                                           create_all_voxels);
//
//   if (regularized) {
//     internal::regularize_octree(lcc);
//   }
//
//   if (IO::is_pretty(std::cout)) {
//     std::cout << "Basic octree generated (level 0/" << max_subdivision_level << ")" << std::endl;
//     std::cout << "Parameters: grid=" << initial_grid_size
//               << ", all_voxels=" << create_all_voxels
//               << ", no_remove_outside=" << no_remove_outside
//               << ", regularized=" << regularized << std::endl;
//   }
//
//   // Unused parameters (for future implementation)
//   (void)max_subdivision_level;
//   (void)no_remove_outside;
// }
// // END OF compute_octree with filename
//
// // Overload for const char*
// template<typename LCC>
// void compute_octree(LCC& lcc,
//                     const char* filename,
//                     unsigned int initial_grid_size = 2,
//                     unsigned int max_subdivision_level = 5,
//                     bool create_all_voxels = false,
//                     bool no_remove_outside = false,
//                     bool regularized = false)
// {
//   compute_octree(lcc, std::string(filename),
//                  initial_grid_size, max_subdivision_level,
//                  create_all_voxels, no_remove_outside, regularized);
// }

/**
 * FaceGraph-based overload: builds the octree from an already loaded mesh.
 * @tparam LCC  model of LinearCellComplex (dim >=3, ambient 3)
 * @tparam FaceGraph model of FaceGraph (triangulated or triangulable)
 * @param regularized If true requests 1-irregular balancing (stub now).
 */
template<typename LCC, typename FaceGraph, class NamedParameters = parameters::Default_named_parameters/*,
         typename std::enable_if<
           std::is_class<typename std::decay<FaceGraph>::type>::value &&
           !std::is_same<typename std::decay<FaceGraph>::type, std::string>::value &&
           !std::is_pointer<typename std::decay<FaceGraph>::type>::value &&
           !std::is_array<FaceGraph>::value
         , int>::type = 0*/>
void compute_octree(LCC& lcc,
                    const FaceGraph& fg,
                    const NamedParameters& np = parameters::default_values()
                    /*unsigned int initial_grid_size = 2,
                    unsigned int max_subdivision_level = 5,
                    bool create_all_voxels = false,
                    bool no_remove_outside = false,
                    bool regularized = false*/)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using Kernel = typename GetGeomTraits<FaceGraph, NamedParameters>::type;

  const int initial_grid_size = parameters::choose_parameter(parameters::get_parameter(np, internal_np::initial_grid_size), 2);
  //const int max_subdivision_level = parameters::choose_parameter(parameters::get_parameter(np, internal_np::max_subdivision_level), 5);

  const bool create_all_voxels = parameters::choose_parameter(parameters::get_parameter(np, internal_np::create_all_voxels), false);
  //const bool no_remove_outside = parameters::choose_parameter(parameters::get_parameter(np, internal_np::no_remove_outside), false);
  const bool regularized = parameters::choose_parameter(parameters::get_parameter(np, internal_np::regularized), false);

  static_assert(LCC::dimension >= 3, "LCC dimension must be >= 3");
  static_assert(LCC::ambient_dimension == 3, "LCC ambient dimension must be 3");

  internal::Simple_AABB_intersector_fgraph<Kernel, FaceGraph> intersector(fg);
  if(intersector.empty()) {
    std::cerr << "Error: empty FaceGraph passed to compute_octree\n";
    return;
  }

  double sx, sy, sz;
  int longestAxis;
  unsigned int initX, initY, initZ;
  internal::compute_initial_grid_size(initial_grid_size, intersector,
                                      longestAxis, sx, sy, sz,
                                      initX, initY, initZ);

  internal::create_initial_hexahedral_grid(lcc, intersector,
                                           sx, sy, sz, initX, initY, initZ,
                                           create_all_voxels);

  if (regularized) internal::regularize_octree(lcc);

  if (IO::is_pretty(std::cout)) {
    std::cout << "Basic octree generated (level 0/"
              //<< max_subdivision_level
              << ")\n"
              << "Parameters: grid=" << initial_grid_size
              << ", all_voxels=" << create_all_voxels
              //<< ", no_remove_outside=" << no_remove_outside
              << ", regularized=" << regularized << std::endl;
  }
}

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_OCTREE_GENERATION_H