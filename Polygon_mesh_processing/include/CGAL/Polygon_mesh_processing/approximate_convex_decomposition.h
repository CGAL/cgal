// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sven Oesau
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_CONVEX_DECOMPOSITION_H
#define CGAL_POLYGON_MESH_PROCESSING_CONVEX_DECOMPOSITION_H

#include <CGAL/license/Polygon_mesh_processing/combinatorial_repair.h>

#include <unordered_set>

#include <CGAL/Iso_cuboid_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/triangle.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/utils_classes.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <atomic>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_priority_queue.h>
#else
#include <unordered_map>
#endif

#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <queue>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

using Vec3_uint = std::array<unsigned int, 3>;

struct Bbox_uint {
  Vec3_uint lower;
  Vec3_uint upper;
  Bbox_uint(const Vec3_uint &lower, const Vec3_uint &upper) : lower(lower), upper(upper) {}
};

enum Grid_cell : int8_t {
  OUTSIDE = -1,
  SURFACE = 0,
  INSIDE = 1
};

void export_grid(const std::string& filename, const Bbox_3& bb, std::vector<int8_t>& grid, const Vec3_uint& grid_size, double voxel_size) {
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };
  std::ofstream stream(filename);

  stream <<
    "ply" << std::endl <<
    "format ascii 1.0" << std::endl <<
    "element vertex " << (grid_size[0] * grid_size[1] * grid_size[2]) << std::endl <<
    "property double x" << std::endl <<
    "property double y" << std::endl <<
    "property double z" << std::endl <<
    "property uchar red" << std::endl <<
    "property uchar green" << std::endl <<
    "property uchar blue" << std::endl <<
    "end_header" << std::endl;

  for (unsigned int x = 0; x < grid_size[0]; x++)
    for (unsigned int y = 0; y < grid_size[1]; y++)
      for (unsigned int z = 0; z < grid_size[2]; z++) {
        stream << (bb.xmin() + (x + 0.5) * voxel_size) << " " << (bb.ymin() + (y + 0.5) * voxel_size) << " " << (bb.zmin() + (z + 0.5) * voxel_size) << " ";
        switch (vox(x, y, z)) {
        case INSIDE:
          stream << "175 175 100" << std::endl;
          break;
        case OUTSIDE:
          stream << "125 125 175" << std::endl;
          break;
        case SURFACE:
          stream << "200 100 100" << std::endl;
          break;
        default:
          stream << "0 0 0" << std::endl;
          break;
        }
      }

  stream.close();
}

template<typename Filter>
void export_grid(const std::string& filename, const Bbox_3& bb, std::vector<int8_t>& grid, const Vec3_uint& grid_size, double voxel_size, Filter &filter) {
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };
  std::ofstream stream(filename);

  std::size_t count = 0;
  for (unsigned int x = 0; x < grid_size[0]; x++)
    for (unsigned int y = 0; y < grid_size[1]; y++)
      for (unsigned int z = 0; z < grid_size[2]; z++)
        if (filter(vox(x, y, z))) count++;

  stream <<
    "ply" << std::endl <<
    "format ascii 1.0" << std::endl <<
    "element vertex " << count << std::endl <<
    "property double x" << std::endl <<
    "property double y" << std::endl <<
    "property double z" << std::endl <<
    "property uchar red" << std::endl <<
    "property uchar green" << std::endl <<
    "property uchar blue" << std::endl <<
    "end_header" << std::endl;

  for (unsigned int x = 0; x < grid_size[0]; x++)
    for (unsigned int y = 0; y < grid_size[1]; y++)
      for (unsigned int z = 0; z < grid_size[2]; z++) {
        if (!filter(vox(x, y, z))) continue;
        stream << (bb.xmin() + (x + 0.5) * voxel_size) << " " << (bb.ymin() + (y + 0.5) * voxel_size) << " " << (bb.zmin() + (z + 0.5) * voxel_size) << " ";
        switch (vox(x, y, z)) {
        case INSIDE:
            stream << "175 175 100" << std::endl;
          break;
        case OUTSIDE:
          stream << "125 125 175" << std::endl;
          break;
        case SURFACE:
          stream << "200 100 100" << std::endl;
          break;
        default:
          stream << "0 0 0" << std::endl;
          break;
        }
      }

  stream.close();
}

void export_voxels(const std::string& filename, const Bbox_3& bb, std::vector<Vec3_uint>& voxels, double voxel_size) {
  std::ofstream stream(filename);

  stream <<
    "ply" << std::endl <<
    "format ascii 1.0" << std::endl <<
    "element vertex " << voxels.size() << std::endl <<
    "property double x" << std::endl <<
    "property double y" << std::endl <<
    "property double z" << std::endl <<
    "end_header" << std::endl;

  for (const Vec3_uint& v : voxels) {
    stream << (bb.xmin() + (v[0] + 0.5) * voxel_size) << " " << (bb.ymin() + (v[1] + 0.5) * voxel_size) << " " << (bb.zmin() + (v[2] + 0.5) * voxel_size) << " ";
  }

  stream.close();
}

template<typename Point_3>
void export_points(const std::string& filename, const Bbox_3& bb, std::vector<Point_3>& points) {
  std::ofstream stream(filename);

  stream <<
    "ply" << std::endl <<
    "format ascii 1.0" << std::endl <<
    "element vertex " << points.size() << std::endl <<
    "property double x" << std::endl <<
    "property double y" << std::endl <<
    "property double z" << std::endl <<
    "end_header" << std::endl;

  for (const Point_3& p : points) {
    stream << (bb.xmin() + p.x()) << " " << (bb.ymin() + p.y()) << " " << (bb.zmin() + p.z()) << " ";
  }

  stream.close();
}

template<typename Range>
void export_points(const std::string& filename, Range& points) {
  std::ofstream stream(filename);
  stream << std::setprecision(18);

  for (const auto& p : points) {
    stream << p.x() << " " << p.y() << " " << p.z() << std::endl;
  }

  stream.close();
}

template<typename Box>
Box box_union(const Box& a, const Box& b) {
  using FT = decltype(a.xmin());
  return Box(
    (std::min<FT>)(a.xmin(), b.xmin()),
    (std::min<FT>)(a.ymin(), b.ymin()),
    (std::min<FT>)(a.zmin(), b.zmin()),
    (std::max<FT>)(a.xmax(), b.xmax()),
    (std::max<FT>)(a.ymax(), b.ymax()),
    (std::max<FT>)(a.zmax(), b.zmax()));
}

template<typename FT>
std::tuple<Vec3_uint, FT> calculate_grid_size(Bbox_3& bb, const std::size_t number_of_voxels) {
  std::size_t max_voxels_axis = std::pow(number_of_voxels, 1.0 / 3.0);
  assert(max_voxels_axis > 3);
  // get longest axis
  FT longest = 0;

  if (bb.x_span() >= bb.y_span() && bb.x_span() >= bb.z_span())
    longest = bb.x_span();
  else if (bb.y_span() >= bb.x_span() && bb.y_span() >= bb.z_span())
    longest = bb.y_span();
  else if (bb.z_span() >= bb.x_span() && bb.z_span() >= bb.y_span())
    longest = bb.z_span();

  const FT voxel_size = longest * FT(1.0 / (max_voxels_axis - 3));

  FT s = 1.5 * voxel_size;

  bb = Bbox_3(to_double(bb.xmin() - s), to_double(bb.ymin() - s), to_double(bb.zmin() - s), to_double(bb.xmax() + s), to_double(bb.ymax() + s), to_double(bb.zmax() + s));

  return { Vec3_uint{static_cast<unsigned int>(to_double(bb.x_span() / voxel_size + 0.5)), static_cast<unsigned int>(to_double(bb.y_span() / voxel_size + 0.5)), static_cast<unsigned int>(to_double(bb.z_span() / voxel_size + 0.5))}, voxel_size};
}

template<typename GeomTraits, typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
const typename GeomTraits::Point_3 &point(typename boost::graph_traits<PolygonMesh>::face_descriptor fd,
  const PolygonMesh& pmesh,
  const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
    vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
      get_const_property_map(CGAL::vertex_point, pmesh));

  return get(vpm, target(halfedge(fd, pmesh), pmesh));
}

template<typename FaceGraph, typename FT>
Bbox_uint grid_bbox_face(const FaceGraph& mesh, const typename boost::graph_traits<FaceGraph>::face_descriptor fd, const Bbox_3& bb, const FT& voxel_size) {
  Bbox_3 face_bb = face_bbox(fd, mesh);
  return Bbox_uint({
    static_cast<unsigned int>((face_bb.xmin() - bb.xmin()) / to_double(voxel_size) - 0.5),
    static_cast<unsigned int>((face_bb.ymin() - bb.ymin()) / to_double(voxel_size) - 0.5),
    static_cast<unsigned int>((face_bb.zmin() - bb.zmin()) / to_double(voxel_size) - 0.5)
    }, {
    static_cast<unsigned int>((face_bb.xmax() - bb.xmin()) / to_double(voxel_size) + 0.5),
    static_cast<unsigned int>((face_bb.ymax() - bb.ymin()) / to_double(voxel_size) + 0.5),
    static_cast<unsigned int>((face_bb.zmax() - bb.zmin()) / to_double(voxel_size) + 0.5)
    });
}

template<typename GeomTraits>
Iso_cuboid_3<GeomTraits> bbox_voxel(const Vec3_uint& voxel, const Bbox_3& bb, const typename GeomTraits::FT& voxel_size) {
  return Bbox_3(
    bb.xmin() + voxel[0] * to_double(voxel_size),
    bb.ymin() + voxel[1] * to_double(voxel_size),
    bb.zmin() + voxel[2] * to_double(voxel_size),
    bb.xmin() + (voxel[0] + 1) * to_double(voxel_size),
    bb.ymin() + (voxel[1] + 1) * to_double(voxel_size),
    bb.zmin() + (voxel[2] + 1) * to_double(voxel_size)
    );
}

void scanline_floodfill(Grid_cell label, std::vector<int8_t>& grid, const Vec3_uint& grid_size, std::deque<Vec3_uint>& todo) {
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };

  while (!todo.empty()) {
    auto [x, y, z0] = todo.front();
    todo.pop_front();
    if (vox(x, y, z0) != INSIDE)
      return;

    bool xneg = false, xpos = false;
    bool yneg = false, ypos = false;
    bool zneg = false, zpos = false;

    // positive direction
    for (unsigned int z = z0; z < grid_size[2]; z++) {
      if (vox(x, y, z) != INSIDE)
        break;

      vox(x, y, z) = label;
      if (x > 0) {
        if (vox(x - 1, y, z) == INSIDE) {
          if (!xneg) {
            xneg = true;
            todo.push_back({ x - 1, y, z });
          }
        }
        else xneg = false;
      }

      if (x < grid_size[0] - 1) {
        if (vox(x + 1, y, z) == INSIDE) {
          if (!xpos) {
            xpos = true;
            todo.push_back({ x + 1, y, z });
          }
        }
        else xpos = false;
      }

      if (y > 0) {
        if (vox(x, y - 1, z) == INSIDE) {
          if (!yneg) {
            yneg = true;
            todo.push_front({ x, y - 1, z });
          }
        }
        else yneg = false;
      }

      if (y < grid_size[1] - 1) {
        if (vox(x, y + 1, z) == INSIDE) {
          if (!ypos) {
            ypos = true;
            todo.push_front({ x, y + 1, z });
          }
        }
        else ypos = false;
      }
    }

    xneg = xpos = yneg = ypos = zneg = zpos = false;

    if (z0 == 0)
      continue;

    for (unsigned int z = z0 - 1; z > 0; z--) {
      if (vox(x, y, z) != INSIDE)
        break;

      vox(x, y, z) = label;
      if (x > 0) {
        if (vox(x - 1, y, z) == INSIDE) {
          if (!xneg) {
            xneg = true;
            todo.push_back({ x - 1, y, z });
          }
        }
        else xneg = false;
      }

      if (x < grid_size[0] - 1) {
        if (vox(x + 1, y, z) == INSIDE) {
          if (!xpos) {
            xpos = true;
            todo.push_back({ x + 1, y, z });
          }
        }
        else xpos = false;
      }

      if (y > 0) {
        if (vox(x, y - 1, z) == INSIDE) {
          if (!yneg) {
            yneg = true;
            todo.push_front({ x, y - 1, z });
          }
        }
        else yneg = false;
      }

      if (y < grid_size[1] - 1) {
        if (vox(x, y + 1, z) == INSIDE) {
          if (!ypos) {
            ypos = true;
            todo.push_front({ x, y + 1, z });
          }
        }
        else ypos = false;
      }
    }
  }
}

// Valid voxel grids separate OUTSIDE from INSIDE via SURFACE
bool check_grid(std::vector<int8_t>& grid, const Vec3_uint& grid_size) {
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };

  std::size_t in = 0, out = 0, surface = 0, other = 0, violations = 0;

  for (unsigned int x = 1; x < grid_size[0] - 1; x++)
    for (unsigned int y = 1; y < grid_size[1] - 1; y++)
      for (unsigned int z = 1; z < grid_size[2] - 1; z++) {
        switch (vox(x, y, z)) {
        case INSIDE:
          if ( vox(x - 1, y, z) == OUTSIDE
            || vox(x + 1, y, z) == OUTSIDE
            || vox(x, y - 1, z) == OUTSIDE
            || vox(x, y + 1, z) == OUTSIDE
            || vox(x, y, z - 1) == OUTSIDE
            || vox(x, y, z + 1) == OUTSIDE) {
            std::cout << "touching I-O: " << x << " " << y << " " << z << std::endl;
            violations++;
          }
          in++;
          break;
        case OUTSIDE:
          if (vox(x - 1, y, z) == INSIDE
            || vox(x + 1, y, z) == INSIDE
            || vox(x, y - 1, z) == INSIDE
            || vox(x, y + 1, z) == INSIDE
            || vox(x, y, z - 1) == INSIDE
            || vox(x, y, z + 1) == INSIDE) {
            std::cout << "touching O-I: " << x << " " << y << " " << z << std::endl;
            violations++;
          }
          out++;
          break;
        case SURFACE:
          surface++;
          break;
        default:
          std::cout << "other " << x << " " << y << " " << z << std::endl;
          other++;
          break;
        }
      }
  std::cout << "i: " << in << " o: " << out << " s: " << surface << " " << other << std::endl;
  std::cout << "violations: " << violations << std::endl;

  return violations == 0;
}

// Only works for closed meshes
void label_floodfill(std::vector<int8_t>& grid, const Vec3_uint& grid_size) {
  // Walk around boundary and start floodfill when voxel label is INSIDE
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };

  std::deque<Vec3_uint> todo;

  // xmin/xmax
  for (unsigned int y = 0; y < grid_size[1]; y++)
    for (unsigned int z = 0; z < grid_size[2]; z++) {
      if (vox(0, y, z) == INSIDE) {
        todo.push_front({0, y, z});
        scanline_floodfill(OUTSIDE, grid, grid_size, todo);
      }

      if (vox(grid_size[0] - 1, y, z) == INSIDE) {
        todo.push_front({ grid_size[0] - 1, y, z });
        scanline_floodfill(OUTSIDE, grid, grid_size, todo);
      }
    }

  // ymin/ymax
  for (unsigned int x = 0; x < grid_size[0]; x++)
    for (unsigned int z = 0; z < grid_size[2]; z++) {
      if (vox(x, 0, z) == INSIDE) {
        todo.push_front({ x, 0, z });
        scanline_floodfill(OUTSIDE, grid, grid_size, todo);
      }

      if (vox(x, grid_size[1] - 1, z) == INSIDE) {
        todo.push_front({ x, grid_size[1] - 1, z });
        scanline_floodfill(OUTSIDE, grid, grid_size, todo);
      }
    }

  // ymin/ymax
  for (unsigned int x = 0; x < grid_size[0]; x++)
    for (unsigned int y = 0; y < grid_size[1]; y++) {
      if (vox(x, y, 0) == INSIDE) {
        todo.push_front({ x, y, 0 });
        scanline_floodfill(OUTSIDE, grid, grid_size, todo);
      }

      if (vox(x, y, grid_size[2] - 1) == INSIDE) {
        todo.push_front({ x, y, grid_size[2] - 1 });
        scanline_floodfill(OUTSIDE, grid, grid_size, todo);
      }
    }
}

void naive_floodfill(std::vector<int8_t>& grid, const Vec3_uint& grid_size) {
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };

  std::queue<Vec3_uint> queue;
  queue.push({0, 0, 0});

  while (!queue.empty()) {
    auto [x, y, z] = queue.front();
    queue.pop();

    if (vox(x, y, z) != INSIDE)
      continue;

    vox(x, y, z) = OUTSIDE;
    if (x > 0)
      if (vox(x - 1, y, z) == INSIDE) {
          queue.push({ x - 1, y, z });
      }

    if (x < grid_size[0] - 1)
      if (vox(x + 1, y, z) == INSIDE) {
        queue.push({ x + 1, y, z });
      }

    if (y > 0)
      if (vox(x, y - 1, z) == INSIDE) {
        queue.push({ x, y - 1, z });
      }

    if (y < grid_size[1] - 1)
      if (vox(x, y + 1, z) == INSIDE) {
        queue.push({ x, y + 1, z });
      }

    if (z > 0)
      if (vox(x, y, z - 1) == INSIDE) {
        queue.push({ x, y, z - 1 });
      }

    if (z < grid_size[2] - 1)
      if (vox(x, y, z + 1) == INSIDE) {
        queue.push({ x, y, z + 1 });
      }
  }
}

template<typename FaceGraph, typename GeomTraits>
void rayshooting_fill(std::vector<int8_t>& grid, const Vec3_uint& grid_size, const Bbox_3& bb, const typename GeomTraits::FT& voxel_size, const FaceGraph& mesh, CGAL::Parallel_tag) {
#ifdef CGAL_LINKED_WITH_TBB
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };
  using face_descriptor = typename boost::graph_traits<FaceGraph>::face_descriptor;

  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  using Ray_3 = typename GeomTraits::Ray_3;

  using Primitive = CGAL::AABB_face_graph_triangle_primitive<FaceGraph>;
  using Traits = CGAL::AABB_traits_3<GeomTraits, Primitive>;
  using Tree = CGAL::AABB_tree<Traits>;
  using Ray_intersection = std::optional<typename Tree::template Intersection_and_primitive_id<Ray_3>::Type>;

  Tree tree(faces(mesh).first, faces(mesh).second, mesh);

  std::array<Vector_3, 6> dirs = { Vector_3{1, 0, 0}, Vector_3{-1, 0, 0}, Vector_3{0, 1, 0}, Vector_3{0,-1, 0}, Vector_3{0, 0, 1}, Vector_3{0, 0,-1} };

  tbb::parallel_for(std::size_t(0), std::size_t(grid_size[0]), [&](const std::size_t x)
  {
    for (std::size_t y = 0; y < grid_size[1]; y++)
      for (std::size_t z = 0; z < grid_size[2]; z++) {
        if (vox(x, y, z) == SURFACE)
          continue;

        Point_3 c(bb.xmin() + (x + 0.5) * voxel_size, bb.ymin() + (y + 0.5) * voxel_size, bb.zmin() + (z + 0.5) * voxel_size);

        unsigned int inside = 0;
        unsigned int outside = 0;

        for (std::size_t i = 0; i < 6; i++) {
          Ray_intersection intersection = tree.first_intersection(Ray_3(c, dirs[i]));
          if (intersection) {
            // A segment intersection is not helpful as it means the triangle normal is orthogonal to the ray
            if (std::get_if<Point_3>(&(intersection->first))) {
              face_descriptor fd = intersection->second;
              Vector_3 n = compute_face_normal(fd, mesh);
              if (dirs[i] * n > 0)
                inside++;
              else
                outside++;
            }
          }
        }

        if (inside >= 3 && outside == 0) {
          vox(x, y, z) = INSIDE;
        }
        else
          vox(x, y, z) = OUTSIDE;
      }
  }
    );
#else
  CGAL_USE(grid);
  CGAL_USE(grid_size);
  CGAL_USE(bb);
  CGAL_USE(voxel_size);
  CGAL_USE(mesh);
#endif
}

template<typename FaceGraph, typename GeomTraits>
void rayshooting_fill(std::vector<int8_t>& grid, const Vec3_uint& grid_size, const Bbox_3& bb, const typename GeomTraits::FT& voxel_size, const FaceGraph& mesh, CGAL::Sequential_tag) {
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };
  using face_descriptor = typename boost::graph_traits<FaceGraph>::face_descriptor;

  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  using Ray_3 = typename GeomTraits::Ray_3;

  using Primitive = CGAL::AABB_face_graph_triangle_primitive<FaceGraph>;
  using Traits = CGAL::AABB_traits_3<GeomTraits, Primitive>;
  using Tree = CGAL::AABB_tree<Traits>;
  using Ray_intersection = std::optional<typename Tree::template Intersection_and_primitive_id<Ray_3>::Type>;

  Tree tree(faces(mesh).first, faces(mesh).second, mesh);

  std::array<Vector_3, 6> dirs = { Vector_3{1, 0, 0}, Vector_3{-1, 0, 0}, Vector_3{0, 1, 0}, Vector_3{0,-1, 0}, Vector_3{0, 0, 1}, Vector_3{0, 0,-1} };

  for (std::size_t x = 0; x < grid_size[0]; x++)
  {
    for (std::size_t y = 0; y < grid_size[1]; y++)
      for (std::size_t z = 0; z < grid_size[2]; z++) {
        if (vox(x, y, z) == SURFACE)
          continue;

        Point_3 c(bb.xmin() + (x + 0.5) * voxel_size, bb.ymin() + (y + 0.5) * voxel_size, bb.zmin() + (z + 0.5) * voxel_size);

        unsigned int inside = 0;
        unsigned int outside = 0;

        for (std::size_t i = 0; i < 6; i++) {
          Ray_intersection intersection = tree.first_intersection(Ray_3(c, dirs[i]));
          if (intersection) {
            // A segment intersection is not helpful as it means the triangle normal is orthogonal to the ray
            if (std::get_if<Point_3>(&(intersection->first))) {
              face_descriptor fd = intersection->second;
              Vector_3 n = compute_face_normal(fd, mesh);
              if (dirs[i] * n > 0)
                inside++;
              else
                outside++;
            }
          }
        }

        if (inside >= 3 && outside == 0) {
          vox(x, y, z) = INSIDE;
        }
        else
          vox(x, y, z) = OUTSIDE;
      }
  }
}

template<typename GeomTraits>
struct Convex_hull {
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;
  Iso_cuboid_3<GeomTraits> bbox;
  FT voxel_volume;
  FT volume;
  FT volume_error;
  std::vector<Point_3> points;
  std::vector<std::array<unsigned int, 3>> indices;

  Convex_hull() noexcept : voxel_volume(0), volume(0), volume_error(0) {}
  Convex_hull(Convex_hull&& o) noexcept {
    bbox = o.bbox;
    voxel_volume = o.voxel_volume;
    volume = o.volume;
    volume_error = o.volume_error;
    points = std::move(o.points);
    indices = std::move(o.indices);
  }

  Convex_hull(const Convex_hull& o) {
    exit(4);
    bbox = o.bbox;
    voxel_volume = o.voxel_volume;
    volume = o.volume;
    volume_error = o.volume_error;
    points = o.points;
    indices = o.indices;
  }

  Convex_hull<GeomTraits>& operator= (const Convex_hull<GeomTraits>& o) {
    exit(3);
    bbox = o.bbox;
    voxel_volume = o.voxel_volume;
    volume = o.volume;
    volume_error = o.volume_error;
    points = o.points;
    indices = o.indices;

    return *this;
  }

  Convex_hull<GeomTraits>& operator= (Convex_hull<GeomTraits>&& o) noexcept {
    bbox = o.bbox;
    voxel_volume = o.voxel_volume;
    volume = o.volume;
    volume_error = o.volume_error;
    points = std::move(o.points);
    indices = std::move(o.indices);

    return *this;
  }
};

template<typename GeomTraits>
struct Candidate {
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;
  std::size_t index;
  std::vector<Vec3_uint> surface;
  std::vector<Vec3_uint> new_surface;
  std::vector<Vec3_uint> inside;
  std::size_t depth;
  Bbox_uint bbox;
  Convex_hull<GeomTraits> ch;

  Candidate() : depth(0), bbox({ 0, 0, 0 }, { 0, 0, 0 }) {index = cidx++;}
  Candidate(std::size_t depth, Bbox_uint &bbox) : depth(depth), bbox(bbox) { index = cidx++; }
private:
  inline static std::size_t cidx = 0;
};

template<typename GeomTraits>
typename GeomTraits::FT volume(const std::vector<typename GeomTraits::Point_3> &pts, const std::vector<std::array<unsigned int, 3>> &indices) {
  ::CGAL::internal::Evaluate<typename GeomTraits::FT> evaluate;

  typename GeomTraits::FT vol = 0;
  typename GeomTraits::Compute_volume_3 cv3;

  typename GeomTraits::Point_3 origin(0, 0, 0);

  for (const std::array<unsigned int, 3> &i : indices) {
    vol += cv3(origin, pts[i[0]], pts[i[1]], pts[i[2]]);
    evaluate(vol);
  }

  return vol;
}


template<typename PointRange, typename Point_3 = typename PointRange::value_type>
void convex_hull(const PointRange& pts, std::vector<Point_3> &hull_points, std::vector<std::array<unsigned int, 3> > &hull_indices) {
  using Mesh = CGAL::Surface_mesh<Point_3>;
  Mesh m;

  convex_hull_3(pts.begin(), pts.end(), m);
//   if (pts.size() < 20) {
//     export_points("convex_hull" + std::to_string(cnt) + ".xyz", pts);
//     std::cout << pts.size() << " " << cnt++ << std::endl;
//     convex_hull_3(pts.begin(), pts.end(), hull_points, hull_indices);// needs bugfix in convex_hull_3 for indexed triangle list
//   }
  hull_points.resize(m.number_of_vertices());
  hull_indices.resize(m.number_of_faces());

  std::size_t idx = 0;
  for (const typename Mesh::vertex_index v : m.vertices())
    hull_points[idx++] = m.point(v);

  idx = 0;
  for (const typename Mesh::face_index f : m.faces()) {
    auto he = m.halfedge(f);
    hull_indices[idx][0] = m.source(he);
    he = m.next(he);
    hull_indices[idx][1] = m.source(he);
    he = m.next(he);
    hull_indices[idx++][2] = m.source(he);
  }
}

void enlarge(Bbox_uint& bbox, const Vec3_uint& v) {
  bbox.lower[0] = (std::min<unsigned int>)(bbox.lower[0], v[0]);
  bbox.lower[1] = (std::min<unsigned int>)(bbox.lower[1], v[1]);
  bbox.lower[2] = (std::min<unsigned int>)(bbox.lower[2], v[2]);
  bbox.upper[0] = (std::max<unsigned int>)(bbox.upper[0], v[0]);
  bbox.upper[1] = (std::max<unsigned int>)(bbox.upper[1], v[1]);
  bbox.upper[2] = (std::max<unsigned int>)(bbox.upper[2], v[2]);
}

template <typename T, typename = std::void_t<>>
struct is_std_hashable : std::false_type {};

template <typename T>
struct is_std_hashable<T, std::void_t<decltype(std::declval<std::hash<T>>()(std::declval<T>()))>> : std::true_type {};

template<typename Point_3, typename H = typename std::conjunction<is_std_hashable<typename Kernel_traits<Point_3>::type::FT>, typename std::is_same<typename Kernel_traits<Point_3>::type::Kernel_tag, CGAL::Cartesian_tag>::type>::type>
struct Set {
};

template<typename Point_3>
struct Set<Point_3, std::true_type> {
  using type = std::unordered_set<Point_3>;
};

template<typename Point_3>
struct Set<Point_3, std::false_type> {
  using type = std::set<Point_3>;
};

template<typename GeomTraits>
void compute_candidate(Candidate<GeomTraits> &c, const Bbox_3& bb, typename GeomTraits::FT voxel_size) {
  using Point_3 = typename GeomTraits::Point_3;
  using FT = typename GeomTraits::FT;

  c.bbox.lower = c.bbox.upper = c.surface[0];

  typename Set<Point_3>::type voxel_points;

  for (const Vec3_uint& v : c.surface) {
    enlarge(c.bbox, v);
    FT xmin = bb.xmin() + FT(v[0]) * voxel_size;
    FT ymin = bb.ymin() + FT(v[1]) * voxel_size;
    FT zmin = bb.zmin() + FT(v[2]) * voxel_size;
    FT xmax = bb.xmin() + FT(v[0] + 1) * voxel_size;
    FT ymax = bb.ymin() + FT(v[1] + 1) * voxel_size;
    FT zmax = bb.zmin() + FT(v[2] + 1) * voxel_size;
    voxel_points.insert(Point_3(xmin, ymin, zmin));
    voxel_points.insert(Point_3(xmin, ymax, zmin));
    voxel_points.insert(Point_3(xmin, ymin, zmax));
    voxel_points.insert(Point_3(xmin, ymax, zmax));
    voxel_points.insert(Point_3(xmax, ymin, zmin));
    voxel_points.insert(Point_3(xmax, ymax, zmin));
    voxel_points.insert(Point_3(xmax, ymin, zmax));
    voxel_points.insert(Point_3(xmax, ymax, zmax));
  }

  for (const Vec3_uint& v : c.new_surface) {
    enlarge(c.bbox, v);
    FT xmin = bb.xmin() + FT(v[0]) * voxel_size;
    FT ymin = bb.ymin() + FT(v[1]) * voxel_size;
    FT zmin = bb.zmin() + FT(v[2]) * voxel_size;
    FT xmax = bb.xmin() + FT(v[0] + 1) * voxel_size;
    FT ymax = bb.ymin() + FT(v[1] + 1) * voxel_size;
    FT zmax = bb.zmin() + FT(v[2] + 1) * voxel_size;
    voxel_points.insert(Point_3(xmin, ymin, zmin));
    voxel_points.insert(Point_3(xmin, ymax, zmin));
    voxel_points.insert(Point_3(xmin, ymin, zmax));
    voxel_points.insert(Point_3(xmin, ymax, zmax));
    voxel_points.insert(Point_3(xmax, ymin, zmin));
    voxel_points.insert(Point_3(xmax, ymax, zmin));
    voxel_points.insert(Point_3(xmax, ymin, zmax));
    voxel_points.insert(Point_3(xmax, ymax, zmax));
  }

  convex_hull(voxel_points, c.ch.points, c.ch.indices);

  c.ch.volume = volume<GeomTraits>(c.ch.points, c.ch.indices);

  CGAL_assertion(c.ch.volume > 0);

  c.ch.voxel_volume = (voxel_size * voxel_size * voxel_size) * FT(double(c.inside.size() + c.surface.size() + c.new_surface.size()));
  c.ch.volume_error = CGAL::abs(c.ch.volume - c.ch.voxel_volume) / c.ch.voxel_volume;
}

template<typename FaceGraph, typename GeomTraits, typename Concurrency_tag>
void fill_grid(Candidate<GeomTraits> &c, std::vector<int8_t> &grid, const FaceGraph &mesh, const Bbox_3& bb, const Vec3_uint& grid_size, const typename GeomTraits::FT& voxel_size, Concurrency_tag tag) {
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };

  for (const typename boost::graph_traits<FaceGraph>::face_descriptor fd : faces(mesh)) {
    Bbox_uint face_bb = grid_bbox_face(mesh, fd, bb, voxel_size);
    CGAL_assertion(face_bb.lower[0] <= face_bb.upper[0]);
    CGAL_assertion(face_bb.lower[1] <= face_bb.upper[1]);
    CGAL_assertion(face_bb.lower[2] <= face_bb.upper[2]);
    CGAL_assertion(face_bb.upper[0] < grid_size[0]);
    CGAL_assertion(face_bb.upper[1] < grid_size[1]);
    CGAL_assertion(face_bb.upper[2] < grid_size[2]);
    for (unsigned int x = face_bb.lower[0]; x <= face_bb.upper[0]; x++)
      for (unsigned int y = face_bb.lower[1]; y <= face_bb.upper[1]; y++)
        for (unsigned int z = face_bb.lower[2]; z <= face_bb.upper[2]; z++) {
          Iso_cuboid_3 box = bbox_voxel<GeomTraits>({ x, y, z }, bb, voxel_size);
          const typename GeomTraits::Point_3 &p = point<GeomTraits>(fd, mesh);
          if (do_intersect(triangle(fd, mesh), box) || box.has_on_bounded_side(p))
            vox(x, y, z) = Grid_cell::SURFACE;
        }
  }

  //export_grid("before.ply", bb, grid, grid_size, voxel_size);
  //auto filterS = [](const int8_t& l) -> bool {return l == SURFACE; };
  //auto filterO = [](const int8_t& l) -> bool {return l == OUTSIDE; };
  //auto filterI = [](const int8_t& l) -> bool {return l == INSIDE; };
  //export_grid("before_surface_ray.ply", bb, grid, grid_size, voxel_size, filterS);
  //export_grid("before_outside.ply", bb, grid, grid_size, voxel_size, filterO);

  if (CGAL::is_closed(mesh)) {
    naive_floodfill(grid, grid_size);
    //check_grid(grid, grid_size);
  }
  else
    rayshooting_fill<FaceGraph, GeomTraits>(grid, grid_size, bb, voxel_size, mesh, tag);

  //export_grid("after_inside_ray.ply", bb, grid, grid_size, voxel_size, filterI);
  //export_grid("after_outside_ray.ply", bb, grid, grid_size, voxel_size, filterO);

  c.bbox.upper = {grid_size[0] - 1, grid_size[1] - 1, grid_size[2] - 1};

  for (unsigned int x = 0; x < grid_size[0]; x++)
    for (unsigned int y = 0; y < grid_size[1]; y++)
      for (unsigned int z = 0; z < grid_size[2]; z++)
        if (vox(x, y, z) == INSIDE)
          c.inside.push_back({x, y, z});
        else if (vox(x, y, z) == SURFACE)
          c.surface.push_back({x, y, z});

  //export_grid("after_inside.ply", bb, grid, grid_size, voxel_size, filterI);
  //export_grid("after_outside.ply", bb, grid, grid_size, voxel_size, filterO);
}

template<typename GeomTraits, typename FaceGraph, typename Concurrency_tag>
void init(Candidate<GeomTraits> &c, const FaceGraph& mesh, std::vector<int8_t>& grid, const Bbox_3& bb, const Vec3_uint& grid_size, const typename GeomTraits::FT& voxel_size, Concurrency_tag tag) {
  internal::fill_grid(c, grid, mesh, bb, grid_size, voxel_size, tag);
  compute_candidate(c, bb, voxel_size);
}

template<typename Candidate_>
void split(std::vector<Candidate_> &candidates, Candidate_& c, unsigned int axis, unsigned int location) {
  //Just split the voxel bbox along 'axis' after voxel index 'location'
  Candidate_ upper(c.depth + 1, c.bbox);
  Candidate_ lower(c.depth + 1, c.bbox);

  CGAL_assertion(c.bbox.lower[axis] < location);
  CGAL_assertion(c.bbox.upper[axis] > location);

  upper.bbox.lower[axis] = location + 1;
  lower.bbox.upper[axis] = location;

  for (const Vec3_uint& v : c.surface) {
    CGAL_assertion(c.bbox.lower[0] <= v[0] && v[0] <= c.bbox.upper[0]);
    CGAL_assertion(c.bbox.lower[1] <= v[1] && v[1] <= c.bbox.upper[1]);
    CGAL_assertion(c.bbox.lower[2] <= v[2] && v[2] <= c.bbox.upper[2]);
    if (location < v[axis])
      upper.surface.push_back(v);
    else
      lower.surface.push_back(v);
  }

  for (const Vec3_uint& v : c.new_surface) {
    CGAL_assertion(c.bbox.lower[0] <= v[0] && v[0] <= c.bbox.upper[0]);
    CGAL_assertion(c.bbox.lower[1] <= v[1] && v[1] <= c.bbox.upper[1]);
    CGAL_assertion(c.bbox.lower[2] <= v[2] && v[2] <= c.bbox.upper[2]);
    if (location < v[axis])
      upper.new_surface.push_back(v);
    else
      lower.new_surface.push_back(v);
  }

  for (const Vec3_uint& v : c.inside) {
    CGAL_assertion(c.bbox.lower[0] <= v[0] && v[0] <= c.bbox.upper[0]);
    CGAL_assertion(c.bbox.lower[1] <= v[1] && v[1] <= c.bbox.upper[1]);
    CGAL_assertion(c.bbox.lower[2] <= v[2] && v[2] <= c.bbox.upper[2]);
    if (location < v[axis]) {
      if ((location + 1) == v[axis])
        upper.new_surface.push_back(v);
      else
        upper.inside.push_back(v);
    }
    else {
      if (location == v[axis])
        lower.new_surface.push_back(v);
      else
        lower.inside.push_back(v);
    }
  }

  if (!upper.surface.empty())
    candidates.emplace_back(std::move(upper));

  if (!lower.surface.empty())
    candidates.emplace_back(std::move(lower));
}

std::size_t concavity(const Vec3_uint& s, const Vec3_uint& e, int axis, std::vector<int8_t>& grid, const Vec3_uint& grid_size) {
  const auto vox = [&grid, &grid_size](const Vec3_uint &v) -> int8_t& {
    return grid[v[2] + (v[1] * grid_size[2]) + (v[0] * grid_size[1] * grid_size[2])];
    };
  std::size_t i;
  for (i = s[axis];i<e[axis];i++) {
    Vec3_uint v = s;
    v[axis] = i;
    if (vox(v) != OUTSIDE)
      break;
  }

  if (i == e[axis])
    return 0;

  std::size_t j;
  for (j = e[axis]; j > s[axis]; j--) {
    Vec3_uint v = s;
    v[axis] = j;
    if (vox(v) != OUTSIDE)
      break;
  }

  std::size_t res = (i - s[axis]) + (e[axis] - j);
  if(res >= grid_size[axis])
    std::cout << "violation!" << std::endl;
  return (i - s[axis]) + (e[axis] - j);
}

void choose_splitting_location_by_concavity(unsigned int& axis, unsigned int& location, const Bbox_uint &bbox, std::vector<int8_t>& grid, const Vec3_uint& grid_size) {
  std::size_t length = bbox.upper[axis] - bbox.lower[axis] + 1;

  CGAL_assertion(length >= 8);
  CGAL_precondition(axis <= 2);

  std::size_t idx0 = 0, idx1 = 1, idx2 = 2;

  switch(axis) {
  case 0:
    idx0 = 1;
    idx1 = 2;
    idx2 = 0;
    break;
  case 1:
    idx0 = 0;
    idx1 = 2;
    idx2 = 1;
    break;
  case 2:
    idx0 = 0;
    idx1 = 1;
    idx2 = 2;
    break;
  }

  std::vector<int> diam(length, 0), diam2(length, 0);

  for (std::size_t i = bbox.lower[idx2]; i <= bbox.upper[idx2]; i++) {
    for (std::size_t j = bbox.lower[idx0]; j <= bbox.upper[idx0]; j++) {
      Vec3_uint s, e;
      s[idx2] = i;
      s[idx1] = bbox.lower[idx1];
      s[idx0] = j;
      e[idx2] = i;
      e[idx1] = bbox.upper[idx1];
      e[idx0] = j;
      diam[i - bbox.lower[idx2]] += concavity(s, e, idx1, grid, grid_size);
    }
  }

  for (std::size_t i = bbox.lower[idx2]; i <= bbox.upper[idx2]; i++) {
    for (std::size_t j = bbox.lower[idx1]; j <= bbox.upper[idx1]; j++) {
      Vec3_uint s, e;
      s[idx0] = bbox.lower[idx0];
      s[idx2] = i;
      s[idx1] = j;
      e[idx0] = bbox.upper[idx0];
      e[idx2] = i;
      e[idx1] = j;
      diam2[i - bbox.lower[idx2]] += concavity(s, e, idx0, grid, grid_size);
    }
  }

  // Skip initial border
  std::size_t border = (length / 10) + 0.5;
  std::size_t pos1, end1 = length;
  int grad = diam[0] - diam[1];
  for (pos1 = 2; pos1 < border; pos1++) {
    int grad1 = diam[pos1 - 1] - diam[pos1];
    // Stop if the gradient flips or flattens significantly
    if (!(grad * grad1 > 0 && grad1 > (grad>>1)))
      break;
    if (grad < grad1)
      grad = grad1;
  }

  grad = diam[length - 1] - diam[length - 2];
  for (end1 = length - 3; end1 > (length - border - 1); end1--) {
    int grad1 = diam[end1 + 1] - diam[end1];
    // Stop if the gradient flips or flattens significantly
    if (!(grad * grad1 > 0 && grad1 > (grad >> 1)))
      break;
    if (grad < grad1)
      grad = grad1;
  }

  std::size_t pos2, end2 = length;
  grad = diam2[0] - diam2[1];
  for (pos2 = 2; pos2 < border; pos2++) {
    int grad2 = diam[pos2 - 1] - diam[pos2];
    // Stop if the gradient flips or flattens significantly
    if (!(grad * grad2 > 0 && grad2 > (grad >> 1)))
      break;
    if (grad < grad2)
      grad = grad2;
  }

  grad = diam2[length - 1] - diam2[length - 2];
  for (end2 = length - 3; end2 > (length - border - 1); end2--) {
    int grad2 = diam[end2 + 1] - diam[end2];
    // Stop if the gradient flips or flattens significantly
    if (!(grad * grad2 > 0 && grad2 > (grad >> 1)))
      break;
    if (grad < grad2)
      grad = grad2;
  }

  std::size_t conc1 = abs(diam[pos1 + 1] - diam[pos1]);
  std::size_t conc2 = abs(diam2[pos2 + 1] - diam2[pos2]);

  for (std::size_t i = pos1;i<end1;i++)
      if (unsigned(abs(diam[i] - diam[i - 1])) > conc1) {
        pos1 = i - 1;
        conc1 = abs(diam[i] - diam[i - 1]);
      }

  for (std::size_t i = pos2; i < end2; i++)
      if (unsigned(abs(diam2[i] - diam2[i - 1])) > conc2) {
        pos2 = i - 1;
        conc2 = abs(diam2[i] - diam2[i - 1]);
      }

  if (conc2 > conc1)
    pos1 = pos2;

  if (pos1 < 2 || (length - 3) < pos1)
    location = (bbox.upper[axis] + bbox.lower[axis]) / 2;
  else
    location = ((conc1 > conc2) ? pos1 : pos2) + bbox.lower[axis];
}

template<typename GeomTraits, typename NamedParameters>
void choose_splitting_plane(Candidate<GeomTraits>& c, unsigned int &axis, unsigned int &location, std::vector<int8_t>& grid, const Vec3_uint& grid_size, const NamedParameters& np) {
  const bool search_concavity = parameters::choose_parameter(parameters::get_parameter(np, internal_np::split_at_concavity), true);
  const std::array<unsigned int, 3> span = {c.bbox.upper[0] - c.bbox.lower[0], c.bbox.upper[1] - c.bbox.lower[1], c.bbox.upper[2] - c.bbox.lower[2]};

  // Split longest axis
  axis = (span[0] >= span[1]) ? 0 : 1;
  axis = (span[axis] >= span[2]) ? axis : 2;

  if (span[axis] >= 8 && search_concavity)
    choose_splitting_location_by_concavity(axis, location, c.bbox, grid, grid_size);
  else
    location = (c.bbox.upper[axis] + c.bbox.lower[axis]) / 2;
}

template<typename GeomTraits, typename NamedParameters>
bool finished(Candidate<GeomTraits> &c, const NamedParameters& np) {
  const typename GeomTraits::FT max_error = parameters::choose_parameter(parameters::get_parameter(np, internal_np::volume_error), 1);

  if (c.ch.volume_error <= max_error)
    return true;

  std::size_t max_span = 0;
  for (std::size_t i = 0;i<3;i++) {
    const std::size_t span = c.bbox.upper[i] - c.bbox.lower[i];
    max_span = (std::max<std::size_t>)(max_span, span);
  }

  if (max_span <= 1)
    return true;

  return false;
}

template<typename GeomTraits, typename NamedParameters>
void recurse(std::vector<Candidate<GeomTraits>>& candidates, std::vector<int8_t>& grid, const Vec3_uint& grid_size, const Bbox_3& bbox, const typename GeomTraits::FT& voxel_size, const NamedParameters& np, CGAL::Parallel_tag) {
#ifdef CGAL_LINKED_WITH_TBB
  const std::size_t max_depth = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_depth), 10);

  std::vector<internal::Candidate<GeomTraits>> final_candidates;

  std::size_t depth = 0;

  while (!candidates.empty() && depth < max_depth) {
    depth++;
    std::vector<Candidate<GeomTraits>> former_candidates = std::move(candidates);
    for (Candidate<GeomTraits>& c : former_candidates) {

      if (finished(c, np)) {
        CGAL::Bbox_3 bbox = CGAL::bbox_3(c.ch.points.begin(), c.ch.points.end());
        bbox.scale(1.1);
        c.ch.bbox = bbox;
        final_candidates.push_back(std::move(c));
        continue;
      }
      unsigned int axis = 0, location = 0;
      choose_splitting_plane(c, axis, location, grid, grid_size, np);
      split(candidates, c, axis, location);
    }
    tbb::parallel_for_each(candidates, [&](Candidate<GeomTraits>& c) {
      compute_candidate(c, bbox, voxel_size);
      });
  }

  std::swap(candidates, final_candidates);
#else
  CGAL_USE(candidates);
  CGAL_USE(grid);
  CGAL_USE(grid_size);
  CGAL_USE(bbox);
  CGAL_USE(voxel_size);
  CGAL_USE(np);
#endif
}

template<typename GeomTraits, typename NamedParameters>
void recurse(std::vector<Candidate<GeomTraits>>& candidates, std::vector<int8_t>& grid, const Vec3_uint& grid_size, const Bbox_3& bbox, const typename GeomTraits::FT& voxel_size, const NamedParameters& np, CGAL::Sequential_tag) {
  const std::size_t max_depth = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_depth), 10);

  std::vector<internal::Candidate<GeomTraits>> final_candidates;

  std::size_t depth = 0;

  while (!candidates.empty() && depth < max_depth) {
    depth++;
    std::vector<Candidate<GeomTraits>> former_candidates = std::move(candidates);
    for (Candidate<GeomTraits>& c : former_candidates) {

      if (finished(c, np)) {
        CGAL::Bbox_3 bbox = CGAL::bbox_3(c.ch.points.begin(), c.ch.points.end());
        bbox.scale(1.1);
        c.ch.bbox = bbox;
        final_candidates.push_back(std::move(c));
        continue;
      }
      unsigned int axis = 0, location = 0;
      choose_splitting_plane(c, axis, location, grid, grid_size, np);
      split(candidates, c, axis, location);
    }

    for (Candidate<GeomTraits>& c : candidates)
      compute_candidate(c, bbox, voxel_size);
  }

  std::swap(candidates, final_candidates);
}

template<typename GeomTraits, typename NamedParameters>
void merge(std::vector<Convex_hull<GeomTraits>>& candidates, const typename GeomTraits::FT& hull_volume, const NamedParameters& np, CGAL::Parallel_tag) {
#ifdef CGAL_LINKED_WITH_TBB
  const std::size_t max_convex_hulls = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_number_of_convex_hulls), 64);
  if (candidates.size() <= max_convex_hulls)
    return;

  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;

  struct Merged_candidate {
    std::size_t ch_a, ch_b;
    int ch;
    FT volume_error;

    bool operator < (const Merged_candidate& other) const {
      return (volume_error > other.volume_error);
    }

    Merged_candidate() : ch_a(-1), ch_b(-1) {}
    Merged_candidate(std::size_t ch_a, std::size_t ch_b) : ch_a(ch_a), ch_b(ch_b) {}
  };

  tbb::concurrent_unordered_map<std::size_t, Convex_hull<GeomTraits>> hulls;
  std::atomic<std::size_t> num_hulls = candidates.size();

  std::unordered_set<std::size_t> keep;

  for (std::size_t i = 0; i < candidates.size(); i++) {
    hulls.emplace(i, std::move(candidates[i]));
    keep.insert(i);
  }

  candidates.clear();
  candidates.reserve(max_convex_hulls);

  std::vector<Merged_candidate> todo;
  tbb::concurrent_priority_queue<Merged_candidate> queue;

  const auto do_merge = [hull_volume, &hulls, &num_hulls](Merged_candidate& m) {
    Convex_hull<GeomTraits>& ci = hulls[m.ch_a];
    Convex_hull<GeomTraits>& cj = hulls[m.ch_b];
    m.ch = num_hulls.fetch_add(1);
    Convex_hull<GeomTraits>& ch = hulls[m.ch];

    ch.bbox = box_union(ci.bbox, cj.bbox);
    std::vector<Point_3> pts(ci.points.begin(), ci.points.end());
    pts.reserve(pts.size() + cj.points.size());
    std::copy(cj.points.begin(), cj.points.end(), std::back_inserter(pts));
    convex_hull(pts, ch.points, ch.indices);

    ch.volume = volume<GeomTraits>(ch.points, ch.indices);

    ch.volume_error = m.volume_error = CGAL::abs(ci.volume + cj.volume - ch.volume) / hull_volume;
    };

  for (std::size_t i : keep) {
    const Convex_hull<GeomTraits>& ci = hulls[i];
    for (std::size_t j : keep) {
      if (j <= i)
        continue;
      const Convex_hull<GeomTraits>& cj = hulls[j];
      if (CGAL::do_intersect(ci.bbox, cj.bbox))
        todo.emplace_back(Merged_candidate(i, j));
      else {
        Merged_candidate m(i, j);
        Bbox_3 bbox = box_union(ci.bbox, cj.bbox).bbox();
        m.ch = -1;
        m.volume_error = CGAL::abs(ci.volume + cj.volume - bbox.x_span() * bbox.y_span() * bbox.z_span()) / hull_volume;
        queue.push(std::move(m));
      }
    }
  }

  // parallel for if available
  tbb::parallel_for_each(todo, do_merge);
  for (Merged_candidate& m : todo)
    queue.push(std::move(m));
  todo.clear();

  while (!queue.empty() && keep.size() > max_convex_hulls) {
    Merged_candidate m;
    while (!queue.try_pop(m) && !queue.empty());

    auto ch_a = hulls.find(m.ch_a);
    if (ch_a == hulls.end())
      continue;

    auto ch_b = hulls.find(m.ch_b);
    if (ch_b == hulls.end())
      continue;

    if (m.ch == -1)
      do_merge(m);

    keep.erase(m.ch_a);
    keep.erase(m.ch_b);

    hulls.unsafe_erase(ch_a);
    hulls.unsafe_erase(ch_b);

    const Convex_hull<GeomTraits>& cj = hulls[m.ch];

    for (std::size_t id : keep) {
      const Convex_hull<GeomTraits>& ci = hulls[id];
      if (CGAL::do_intersect(ci.bbox, cj.bbox))
        todo.emplace_back(Merged_candidate(id, m.ch));
      else {
        Merged_candidate merged(id, m.ch);
        Bbox_3 bbox = box_union(ci.bbox, cj.bbox).bbox();
        merged.ch = -1;
        merged.volume_error = CGAL::abs(ci.volume + cj.volume - bbox.x_span() * bbox.y_span() * bbox.z_span()) / hull_volume;
        queue.push(std::move(merged));
      }
    }

    keep.insert(m.ch);

    tbb::parallel_for_each(todo, do_merge);
    for (Merged_candidate& m : todo)
      queue.push(std::move(m));
    todo.clear();
  }

  num_hulls = 0;

  for (std::size_t i : keep)
    candidates.push_back(std::move(hulls[i]));
#else
  CGAL_USE(candidates);
  CGAL_USE(hull_volume);
  CGAL_USE(np);
  assert(false);
#endif
}

template<typename GeomTraits, typename NamedParameters>
void merge(std::vector<Convex_hull<GeomTraits>>& candidates, const typename GeomTraits::FT& hull_volume, const NamedParameters& np, CGAL::Sequential_tag) {
  const std::size_t max_convex_hulls = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_number_of_convex_hulls), 64);
  if (candidates.size() <= max_convex_hulls)
    return;

  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;

  struct Merged_candidate {
    std::size_t ch_a, ch_b;
    int ch;
    FT volume_error;

    bool operator < (const Merged_candidate& other) const {
      return (volume_error > other.volume_error);
    }

    Merged_candidate() : ch_a(-1), ch_b(-1) {}
    Merged_candidate(std::size_t ch_a, std::size_t ch_b) : ch_a(ch_a), ch_b(ch_b) {}
  };

  std::unordered_map<std::size_t, Convex_hull<GeomTraits>> hulls;
  std::size_t num_hulls = candidates.size();

  std::unordered_set<std::size_t> keep;

  for (std::size_t i = 0; i < candidates.size(); i++) {
    hulls.emplace(i, std::move(candidates[i]));
    keep.insert(i);
  }

  candidates.clear();
  candidates.reserve(max_convex_hulls);

  std::priority_queue<Merged_candidate> queue;

  const auto do_merge = [hull_volume, &hulls, &num_hulls](Merged_candidate& m) {
    Convex_hull<GeomTraits>& ci = hulls[m.ch_a];
    Convex_hull<GeomTraits>& cj = hulls[m.ch_b];

    m.ch = num_hulls++;
    Convex_hull<GeomTraits>& ch = hulls[m.ch];

    ch.bbox = box_union(ci.bbox, cj.bbox);
    std::vector<Point_3> pts(ci.points.begin(), ci.points.end());
    pts.reserve(pts.size() + cj.points.size());
    std::copy(cj.points.begin(), cj.points.end(), std::back_inserter(pts));
    convex_hull(pts, ch.points, ch.indices);

    ch.volume = volume<GeomTraits>(ch.points, ch.indices);

    ch.volume_error = m.volume_error = CGAL::abs(ci.volume + cj.volume - ch.volume) / hull_volume;
    };

  for (std::size_t i : keep) {
    const Convex_hull<GeomTraits>& ci = hulls[i];
    for (std::size_t j : keep) {
      if (j <= i)
        continue;
      const Convex_hull<GeomTraits>& cj = hulls[j];
      if (CGAL::do_intersect(ci.bbox, cj.bbox)) {
        Merged_candidate m(i, j);

        m.ch = num_hulls++;
        Convex_hull<GeomTraits>& ch = hulls[m.ch];
        ch.bbox = box_union(ci.bbox, cj.bbox);
        std::vector<Point_3> pts(ci.points.begin(), ci.points.end());
        pts.reserve(pts.size() + cj.points.size());
        std::copy(cj.points.begin(), cj.points.end(), std::back_inserter(pts));
        convex_hull(pts, ch.points, ch.indices);

        ch.volume = volume<GeomTraits>(ch.points, ch.indices);

        ch.volume_error = m.volume_error = CGAL::abs(ci.volume + cj.volume - ch.volume) / hull_volume;
        queue.push(std::move(m));
      }
      else {
        Merged_candidate m(i, j);
        Bbox_3 bbox = box_union(ci.bbox, cj.bbox).bbox();
        m.ch = -1;
        m.volume_error = CGAL::abs(ci.volume + cj.volume - bbox.x_span() * bbox.y_span() * bbox.z_span()) / hull_volume;
        queue.push(std::move(m));
      }
    }
  }

  while (!queue.empty() && keep.size() > max_convex_hulls) {
    Merged_candidate m = queue.top();
    queue.pop();

    auto ch_a = hulls.find(m.ch_a);
    if (ch_a == hulls.end())
      continue;

    auto ch_b = hulls.find(m.ch_b);
    if (ch_b == hulls.end())
      continue;

    if (m.ch == -1)
      do_merge(m);

    keep.erase(m.ch_a);
    keep.erase(m.ch_b);

    hulls.erase(ch_a);
    hulls.erase(ch_b);

    const Convex_hull<GeomTraits>& cj = hulls[m.ch];

    for (std::size_t id : keep) {
      const Convex_hull<GeomTraits>& ci = hulls[id];
      if (CGAL::do_intersect(ci.bbox, cj.bbox)) {
        Merged_candidate merged(id, m.ch);

        merged.ch = num_hulls++;
        Convex_hull<GeomTraits>& ch = hulls[merged.ch];
        ch.bbox = box_union(ci.bbox, cj.bbox);
        std::vector<Point_3> pts(ci.points.begin(), ci.points.end());
        pts.reserve(pts.size() + cj.points.size());
        std::copy(cj.points.begin(), cj.points.end(), std::back_inserter(pts));
        convex_hull(pts, ch.points, ch.indices);

        ch.volume = volume<GeomTraits>(ch.points, ch.indices);

        ch.volume_error = merged.volume_error = CGAL::abs(ci.volume + cj.volume - ch.volume) / hull_volume;
        queue.push(std::move(merged));
      }
      else {
        Merged_candidate merged(id, m.ch);
        Bbox_3 bbox = box_union(ci.bbox, cj.bbox).bbox();
        merged.ch = -1;
        merged.volume_error = CGAL::abs(ci.volume + cj.volume - bbox.x_span() * bbox.y_span() * bbox.z_span()) / hull_volume;
        queue.push(std::move(merged));
      }
    }

    keep.insert(m.ch);
  }

  num_hulls = 0;

  for (std::size_t i : keep)
    candidates.push_back(std::move(hulls[i]));
}

}
/**
 * \ingroup PMP_convex_decomposition_grp
 *
 * \brief approximates the input mesh by a number of convex hulls. The input mesh is voxelized and the voxel intersecting with the mesh are labeled as surface.
 *        The remaining voxels are labeled as outside or inside by floodfill, in case the input mesh is closed, or by ray shooting along the axis. A voxel is only
 *        labeled as inside if at least 3 faces facing away from the voxel have been hit and no face facing towards the voxel. In a next step, the convex hull of the mesh
 *        is hierarchically split until the `volume_error` threshold is satisfied.
 *        Afterwards, a greedy pair-wise merging combines smaller convex hulls until the given number of convex hulls is met.
 *
 * \tparam FaceGraph a model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
 *
 * \tparam OutputIterator must be an output iterator accepting variables of type `std::pair<std::vector<geom_traits::Point_3>, std::vector<std::array<unsigned int, 3> > >`.
 *
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param mesh the input mesh to approximate by convex hulls
 * \param out_hulls output iterator into which convex hulls are recorded
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{maximum_number_of_voxels}
 *     \cgalParamDescription{gives an upper bound of the number of voxels. The longest bounding box side will have a length of the cubic root of `maximum_number_of_voxels`. Cannot be smaller than 64. }
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{1000000}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{maximum_depth}
 *     \cgalParamDescription{maximum depth of hierarchical splits}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{10}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{maximum_number_of_convex_hulls}
 *     \cgalParamDescription{maximum number of convex hulls produced by the method}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{16}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{volume_error}
 *     \cgalParamDescription{maximum difference in fraction of volume of the local convex hull with the sum of voxel volumes. If surpassed, the convex hull will be split.}
 *     \cgalParamType{double}
 *     \cgalParamDefault{0.01}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{split_at_concavity}
 *     \cgalParamDescription{split the local box of the convex hull in the mid of the longest axis (faster) or search the concavity along the longest axis of the bounding box for splitting.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{true}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `mesh`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `FaceGraph`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{concurrency_tag}
 *     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *     \cgalParamDefault{`CGAL::Parallel_if_available_tag`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type of `FaceGraph`, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 * \sa `CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh()`
 */
template<typename FaceGraph, typename OutputIterator, typename NamedParameters = parameters::Default_named_parameters>
std::size_t approximate_convex_decomposition(const FaceGraph& mesh, OutputIterator out_hulls, const NamedParameters& np = parameters::default_values()) {
  using Geom_traits = typename GetGeomTraits<FaceGraph, NamedParameters>::type;

  using FT = typename Geom_traits::FT;
  const std::size_t num_voxels = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_number_of_voxels), 1000000);
  using Concurrency_tag = typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t, NamedParameters, Parallel_if_available_tag>::type;

#ifndef CGAL_LINKED_WITH_TBB
  if constexpr (std::is_same_v<Concurrency_tag, Parallel_tag> || std::is_same_v<Concurrency_tag, Parallel_if_available_tag>) {
    CGAL_error_msg("CGAL was not compiled with TBB support. Use Sequential_tag instead.");
    return 0;
  }
#endif

  const std::size_t max_convex_hulls = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_number_of_convex_hulls), 64);
  assert(max_convex_hulls > 0);

  if (max_convex_hulls == 1) {
    internal::Convex_hull< Geom_traits> ch;
    using Mesh = Surface_mesh<typename Geom_traits::Point_3>;
    Mesh m;
    convex_hull_3(mesh, m);

    ch.points.resize(m.number_of_vertices());
    ch.indices.resize(m.number_of_faces());

    std::size_t idx = 0;
    for (const typename Mesh::vertex_index v : m.vertices())
      ch.points[idx++] = m.point(v);

    idx = 0;
    for (const typename Mesh::face_index f : m.faces()) {
      auto he = m.halfedge(f);
      ch.indices[idx][0] = m.source(he);
      he = m.next(he);
      ch.indices[idx][1] = m.source(he);
      he = m.next(he);
      ch.indices[idx++][2] = m.source(he);
    }

    *out_hulls = std::make_pair(std::move(ch.points), std::move(ch.indices));
    return 1;
  }

  Bbox_3 bb = bbox(mesh);
  const auto [grid_size, voxel_size] = internal::calculate_grid_size<FT>(bb, num_voxels);

  std::vector<int8_t> grid(grid_size[0] * grid_size[1] * grid_size[2], internal::Grid_cell::INSIDE);

  std::vector<internal::Candidate<Geom_traits>> candidates(1);
  init(candidates[0], mesh, grid, bb, grid_size, voxel_size, Concurrency_tag());

  const FT hull_volume = candidates[0].ch.volume;

  recurse(candidates, grid, grid_size, bb, voxel_size, np, Concurrency_tag());

  std::vector<internal::Convex_hull<Geom_traits>> hulls;
  for (internal::Candidate<Geom_traits> &c : candidates)
    hulls.emplace_back(std::move(c.ch));

  candidates.clear();

  // merge until target number is reached
  merge(hulls, hull_volume, np, Concurrency_tag());

  for (std::size_t i = 0; i < hulls.size(); i++)
    *out_hulls++ = std::make_pair(std::move(hulls[i].points), std::move(hulls[i].indices));

  return hulls.size();
}
}
}

#endif