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

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/triangle.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/utils_classes.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_soup_io.h>

#undef CGAL_LINKED_WITH_TBB

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

struct hash {
  std::size_t operator()(const std::array<unsigned int, 3> &a) const {
    return (((a[0] * 1193) ^ a[1]) * 1212) ^ a[2];
  }
};

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
  UNKNOWN = 0,
  INSIDE = 1,
  SURFACE = 2,
  UNSET = 3
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

void export_voxels(const std::string& filename, const Bbox_3& bb, std::vector<Vec3_uint>& voxels, std::vector<int8_t>& grid, const Vec3_uint& grid_size, double voxel_size) {
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };
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

  stream <<
    "ply" << std::endl <<
    "format ascii 1.0" << std::endl <<
    "element vertex " << points.size() << std::endl <<
    "property double x" << std::endl <<
    "property double y" << std::endl <<
    "property double z" << std::endl <<
    "end_header" << std::endl;

  for (const auto& p : points) {
    stream << p.x() << " " << p.y() << " " << p.z() << std::endl;
  }

  stream.close();
}


template<typename FT>
std::tuple<Vec3_uint, FT> calculate_grid_size(Bbox_3& bb, const std::size_t number_of_voxels) {
  std::size_t max_voxels_axis = std::pow(number_of_voxels, 1.0 / 3.0);
  // get largest axis
  FT largest = 0;

  if (bb.x_span() >= bb.y_span() && bb.x_span() >= bb.z_span())
    largest = bb.x_span();
  else if (bb.y_span() >= bb.x_span() && bb.y_span() >= bb.z_span())
    largest = bb.y_span();
  else if (bb.z_span() >= bb.x_span() && bb.z_span() >= bb.y_span())
    largest = bb.z_span();

  const FT voxel_size = largest / (max_voxels_axis - 3);

  FT s = 1.5 * voxel_size;

  bb = Bbox_3(bb.xmin() - s, bb.ymin() - s, bb.zmin() - s, bb.xmax() + s, bb.ymax() + s, bb.zmax() + s);

  return { Vec3_uint{static_cast<unsigned int>(bb.x_span() / voxel_size + 0.5), static_cast<unsigned int>(bb.y_span() / voxel_size + 0.5), static_cast<unsigned int>(bb.z_span() / voxel_size + 0.5)}, voxel_size};
}


template<typename FaceGraph, typename FT>
Bbox_uint grid_bbox_face(const FaceGraph& mesh, const typename boost::graph_traits<FaceGraph>::face_descriptor fd, const Bbox_3& bb, const FT& voxel_size) {
  Bbox_3 face_bb = face_bbox(fd, mesh);
  return Bbox_uint({
    static_cast<unsigned int>((face_bb.xmin() - bb.xmin()) / voxel_size - 0.5),
    static_cast<unsigned int>((face_bb.ymin() - bb.ymin()) / voxel_size - 0.5),
    static_cast<unsigned int>((face_bb.zmin() - bb.zmin()) / voxel_size - 0.5)
    }, {
    static_cast<unsigned int>((face_bb.xmax() - bb.xmin()) / voxel_size + 0.5),
    static_cast<unsigned int>((face_bb.ymax() - bb.ymin()) / voxel_size + 0.5),
    static_cast<unsigned int>((face_bb.zmax() - bb.zmin()) / voxel_size + 0.5)
    });
}

template<typename FT>
Bbox_3 bbox_voxel(const Vec3_uint& voxel, const Bbox_3& bb, const FT& voxel_size) {
  return Bbox_3(
    bb.xmin() + voxel[0] * voxel_size,
    bb.ymin() + voxel[1] * voxel_size,
    bb.zmin() + voxel[2] * voxel_size,
    bb.xmax() + (voxel[0] + 1) * voxel_size,
    bb.ymax() + (voxel[1] + 1) * voxel_size,
    bb.zmax() + (voxel[2] + 1) * voxel_size
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
      if (x > 0)
        if (vox(x - 1, y, z) == INSIDE) {
          if (!xneg) {
            xneg = true;
            todo.push_back({ x - 1, y, z });
          }
        }
        else xneg = false;

      if (x < grid_size[0] - 1)
        if (vox(x + 1, y, z) == INSIDE) {
          if (!xpos) {
            xpos = true;
            todo.push_back({ x + 1, y, z });
          }
        }
        else xpos = false;

      if (y > 0)
        if (vox(x, y - 1, z) == INSIDE) {
          if (!yneg) {
            yneg = true;
            todo.push_front({ x, y - 1, z });
          }
        }
        else yneg = false;

      if (y < grid_size[1] - 1)
        if (vox(x, y + 1, z) == INSIDE) {
          if (!ypos) {
            ypos = true;
            todo.push_front({ x, y + 1, z });
          }
        }
        else ypos = false;
    }

    xneg = xpos = yneg = ypos = zneg = zpos = false;

    if (z0 == 0)
      continue;

    for (unsigned int z = z0 - 1; z > 0; z--) {
      if (vox(x, y, z) != INSIDE)
        break;

      vox(x, y, z) = label;
      if (x > 0)
        if (vox(x - 1, y, z) == INSIDE) {
          if (!xneg) {
            xneg = true;
            todo.push_back({ x - 1, y, z });
          }
        }
        else xneg = false;

      if (x < grid_size[0] - 1)
        if (vox(x + 1, y, z) == INSIDE) {
          if (!xpos) {
            xpos = true;
            todo.push_back({ x + 1, y, z });
          }
        }
        else xpos = false;

      if (y > 0)
        if (vox(x, y - 1, z) == INSIDE) {
          if (!yneg) {
            yneg = true;
            todo.push_front({ x, y - 1, z });
          }
        }
        else yneg = false;

      if (y < grid_size[1] - 1)
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

template<typename GeomTraits>
struct Convex_hull {
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;
  Bbox_3 bbox;
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
    bbox = o.bbox;
    voxel_volume = o.voxel_volume;
    volume = o.volume;
    volume_error = o.volume_error;
    points = o.points;
    indices = o.indices;
  }

  Convex_hull<GeomTraits>& operator= (const Convex_hull<GeomTraits>& o) {
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
  std::vector<Vec3_uint> surface;
  std::vector<Vec3_uint> new_surface;
  std::vector<Vec3_uint> inside;
  std::size_t depth;
  Bbox_uint bbox;
  Convex_hull<GeomTraits> ch;

  Candidate() : depth(0), bbox({ 0, 0, 0 }, { 0, 0, 0 }) {}
  Candidate(std::size_t depth, Bbox_uint &bbox) : depth(depth), bbox(bbox) {}
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
  //std::vector<Point_3> pts;
  //std::vector<std::vector<std::size_t> > indices;
  //convex_hull_3(voxel_points.begin(), voxel_points.end(), pts, indices);// Why does this not work?
  //convex_hull_3(voxel_points.begin(), voxel_points.end(), c.hull_points, c.hull_indices);// Why does this not work?
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

template<typename GeomTraits>
void compute_candidate(Candidate<GeomTraits> &c, const Bbox_3& bb, typename GeomTraits::FT voxel_size) {
  using Point_3 = typename GeomTraits::Point_3;
  using FT = typename GeomTraits::FT;

  std::unordered_set<Point_3> voxel_points;

  // Is it more efficient than using std::vector with plenty of duplicates?
  for (const Vec3_uint& v : c.surface) {
    FT xmin = bb.xmin() + v[0] * voxel_size;
    FT ymin = bb.ymin() + v[1] * voxel_size;
    FT zmin = bb.zmin() + v[2] * voxel_size;
    FT xmax = bb.xmin() + (v[0] + 1) * voxel_size;
    FT ymax = bb.ymin() + (v[1] + 1) * voxel_size;
    FT zmax = bb.zmin() + (v[2] + 1) * voxel_size;
    voxel_points.insert(Point_3(xmin, ymin, zmin));
    voxel_points.insert(Point_3(xmin, ymax, zmin));
    voxel_points.insert(Point_3(xmin, ymin, zmax));
    voxel_points.insert(Point_3(xmin, ymax, zmax));
    voxel_points.insert(Point_3(xmax, ymin, zmin));
    voxel_points.insert(Point_3(xmax, ymax, zmin));
    voxel_points.insert(Point_3(xmax, ymin, zmax));
    voxel_points.insert(Point_3(xmax, ymax, zmax));
  }

  //export_points("surface_points.ply", voxel_points);

  convex_hull(voxel_points, c.ch.points, c.ch.indices);

  c.ch.volume = volume<GeomTraits>(c.ch.points, c.ch.indices);

  assert(c.ch.volume > 0);

  c.ch.voxel_volume = (voxel_size * voxel_size * voxel_size) * (c.inside.size() + c.surface.size() + c.new_surface.size());
  c.ch.volume_error = CGAL::abs(c.ch.volume - c.ch.voxel_volume) / c.ch.voxel_volume;
}

template<typename FaceGraph, typename GeomTraits>
void fill_grid(Candidate<GeomTraits> &c, std::vector<int8_t> &grid, const FaceGraph &mesh, const Bbox_3& bb, const Vec3_uint& grid_size, const typename GeomTraits::FT& voxel_size) {
  const auto vox = [&grid, &grid_size](unsigned int x, unsigned int y, unsigned int z) -> int8_t& {
    return grid[z + (y * grid_size[2]) + (x * grid_size[1] * grid_size[2])];
    };

  for (const typename boost::graph_traits<FaceGraph>::face_descriptor fd : faces(mesh)) {
    Bbox_uint face_bb = grid_bbox_face(mesh, fd, bb, voxel_size);
    assert(face_bb.lower[0] <= face_bb.upper[0]);
    assert(face_bb.lower[1] <= face_bb.upper[1]);
    assert(face_bb.lower[2] <= face_bb.upper[2]);
    assert(face_bb.upper[0] < grid_size[0]);
    assert(face_bb.upper[1] < grid_size[1]);
    assert(face_bb.upper[2] < grid_size[2]);
    for (unsigned int x = face_bb.lower[0]; x <= face_bb.upper[0]; x++)
      for (unsigned int y = face_bb.lower[1]; y <= face_bb.upper[1]; y++)
        for (unsigned int z = face_bb.lower[2]; z <= face_bb.upper[2]; z++)
          if (do_intersect(triangle(fd, mesh), bbox_voxel({x, y, z}, bb, voxel_size)))
            vox(x, y, z) = Grid_cell::SURFACE;
  }

  //export_grid("before.ply", bb, grid, grid_size, voxel_size);

  // For now, only do floodfill
  //label_floodfill(grid, grid_size);
  naive_floodfill(grid, grid_size);

  c.bbox.upper = grid_size;

  for (unsigned int x = 0; x < grid_size[0]; x++)
    for (unsigned int y = 0; y < grid_size[1]; y++)
      for (unsigned int z = 0; z < grid_size[2]; z++)
        if (vox(x, y, z) == INSIDE)
          c.inside.push_back({x, y, z});
        else if (vox(x, y, z) == SURFACE)
          c.surface.push_back({ x, y, z });

  check_grid(grid, grid_size);

  //export_grid("after_fill.ply", bb, grid, grid_size, voxel_size);
}

template<typename GeomTraits, typename FaceGraph>
void init(Candidate<GeomTraits> &c, const FaceGraph& mesh, std::vector<int8_t>& grid, const Bbox_3& bb, const Vec3_uint& grid_size, const typename GeomTraits::FT& voxel_size) {
  internal::fill_grid(c, grid, mesh, bb, grid_size, voxel_size);
  compute_candidate(c, bb, voxel_size);
}

template<typename GeomTraits, typename NamedParameters>
void split(std::vector<Candidate<GeomTraits> > &candidates, Candidate<GeomTraits>& c, unsigned int axis, unsigned int location, std::vector<int8_t>& grid, const Vec3_uint& grid_size, const Bbox_3& bbox, const typename GeomTraits::FT& voxel_size, const NamedParameters& np) {
  //Just split the voxel bbox along 'axis' after voxel index 'location'
  Candidate<GeomTraits> upper(c.depth + 1, c.bbox);
  Candidate<GeomTraits> lower(c.depth + 1, c.bbox);

  upper.bbox.lower[axis] = location;
  lower.bbox.upper[axis] = location - 1;

  for (const Vec3_uint& v : c.surface) {
    assert(c.bbox.lower[0] <= v[0] && v[0] <= c.bbox.upper[0]);
    assert(c.bbox.lower[1] <= v[1] && v[1] <= c.bbox.upper[1]);
    assert(c.bbox.lower[2] <= v[2] && v[2] <= c.bbox.upper[2]);
    if (location <= v[axis])
      upper.surface.push_back(v);
    else
      lower.surface.push_back(v);
  }

  for (const Vec3_uint& v : c.new_surface) {
    assert(c.bbox.lower[0] <= v[0] && v[0] <= c.bbox.upper[0]);
    assert(c.bbox.lower[1] <= v[1] && v[1] <= c.bbox.upper[1]);
    assert(c.bbox.lower[2] <= v[2] && v[2] <= c.bbox.upper[2]);
    if (location <= v[axis])
      upper.new_surface.push_back(v);
    else
      lower.new_surface.push_back(v);
  }

  for (const Vec3_uint& v : c.inside) {
    assert(c.bbox.lower[0] <= v[0] && v[0] <= c.bbox.upper[0]);
    assert(c.bbox.lower[1] <= v[1] && v[1] <= c.bbox.upper[1]);
    assert(c.bbox.lower[2] <= v[2] && v[2] <= c.bbox.upper[2]);
    if (location <= v[axis]) {
      if (upper.bbox.lower[axis] == v[axis])
        upper.new_surface.push_back(v);
      else
        upper.inside.push_back(v);
    }
    else {
      if (lower.bbox.upper[axis] == v[axis])
        lower.new_surface.push_back(v);
      else
        lower.inside.push_back(v);
      //new_surface is used for convex hull calculation
    }
  }

  if (!upper.surface.empty()) {
    compute_candidate(upper, bbox, voxel_size);
    candidates.emplace_back(std::move(upper));
  }

  if (!lower.surface.empty()) {
    compute_candidate(lower, bbox, voxel_size);
    candidates.emplace_back(std::move(lower));
  }
}

template<typename GeomTraits, typename NamedParameters>
void choose_splitting_plane(Candidate<GeomTraits>& c, unsigned int &axis, unsigned int &location, std::vector<int8_t>& grid, const Vec3_uint& grid_size, const NamedParameters& np) {
  //Just split the voxel bbox along 'axis' after voxel index 'location'
  const std::array<unsigned int, 3> span = {c.bbox.upper[0] - c.bbox.lower[0], c.bbox.upper[1] - c.bbox.lower[1], c.bbox.upper[2] - c.bbox.lower[2]};

  // Split largest axis
  axis = (span[0] >= span[1]) ? 0 : 1;
  axis = (span[axis] >= span[2]) ? axis : 2;

  location = (c.bbox.upper[axis] + c.bbox.lower[axis]) / 2;
}

template<typename GeomTraits, typename NamedParameters>
bool finished(Candidate<GeomTraits> &c, const NamedParameters& np) {
  const typename GeomTraits::FT max_error = parameters::choose_parameter(parameters::get_parameter(np, internal_np::volume_error), 1);
  const std::size_t max_depth = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_depth), 10);

  if (c.ch.volume_error <= max_error)
    return true;

  if (c.depth >= max_depth)
    return true;

  return false;
}

template<typename GeomTraits, typename NamedParameters>
void recurse(std::vector<Candidate<GeomTraits>>& candidates, std::vector<int8_t>& grid, const Vec3_uint& grid_size, const Bbox_3 &bbox, const typename GeomTraits::FT &voxel_size, const NamedParameters& np) {
  using FT = typename GeomTraits::FT;
  const FT max_error = parameters::choose_parameter(parameters::get_parameter(np, internal_np::volume_error), 1);
  const std::size_t max_depth = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_depth), 10);

  std::vector<internal::Candidate<GeomTraits>> final_candidates;

  while (!candidates.empty()) {
    std::vector<Candidate<GeomTraits>> former_candidates = std::move(candidates);
    for (Candidate<GeomTraits>& c : former_candidates) {
      // check loop conditions here?
      if (finished(c, np)) {
        c.ch.bbox = CGAL::bbox_3(c.ch.points.begin(), c.ch.points.end());
        c.ch.bbox.scale(1.1); // Enlarge bounding boxes by a small factor for the following merge
        final_candidates.push_back(std::move(c));
        continue;
      }
      unsigned int axis = 0, location = 0;
      choose_splitting_plane(c, axis, location, grid, grid_size, np);
      split(candidates, c, axis, location, grid, grid_size, bbox, voxel_size, np);
    }
  }

  std::swap(candidates, final_candidates);
}

template<typename GeomTraits, typename NamedParameters>
void merge(std::vector<Convex_hull<GeomTraits>>& candidates, std::vector<int8_t>& grid, const Vec3_uint& grid_size, const Bbox_3& bbox, const typename GeomTraits::FT& voxel_size, const typename GeomTraits::FT &hull_volume, const NamedParameters& np) {
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

  // Consider all non-equal pairs
  std::size_t max_merge_candidates = (candidates.size() * (candidates.size() - 1)) >> 1;

#ifdef CGAL_LINKED_WITH_TBB
  tbb::concurrent_unordered_map<std::size_t, Convex_hull<GeomTraits>> hulls;
  std::atomic<std::size_t> num_hulls = candidates.size();
#else
  std::unordered_map<std::size_t, Convex_hull<GeomTraits>> hulls;
  std::size_t num_hulls = candidates.size();
#endif

  std::unordered_set<std::size_t> keep;

  for (std::size_t i = 0;i<candidates.size();i++) {
    hulls.emplace(i, std::move(candidates[i]));
    keep.insert(i);
  }

  candidates.clear();
  candidates.reserve(max_convex_hulls);

#ifdef CGAL_LINKED_WITH_TBB
  std::vector<Merged_candidate> todo;
  tbb::concurrent_priority_queue<Merged_candidate> queue;
#else
  std::priority_queue<Merged_candidate> queue;
#endif

  const auto do_merge = [hull_volume, &hulls, &num_hulls, &queue](Merged_candidate &m) {
    Convex_hull<GeomTraits>& ci = hulls[m.ch_a];
    Convex_hull<GeomTraits>& cj = hulls[m.ch_b];
#ifdef CGAL_LINKED_WITH_TBB
    m.ch = num_hulls.fetch_add(1);
    Convex_hull<GeomTraits>& ch = hulls[m.ch];
#else
    m.ch = num_hulls++;
    Convex_hull<GeomTraits>& ch = hulls[m.ch];
#endif
    ch.bbox = ci.bbox + cj.bbox;
    std::vector<Point_3> pts(ci.points.begin(), ci.points.end());
    pts.reserve(pts.size() + cj.points.size());
    std::copy(cj.points.begin(), cj.points.end(), std::back_inserter(pts));
    convex_hull(pts, ch.points, ch.indices);

    ch.volume = volume<GeomTraits>(ch.points, ch.indices);

    ch.volume_error = m.volume_error = CGAL::abs(ci.volume + cj.volume - ch.volume) / hull_volume;
  };

  //queue.reserve(max_merge_candidates);

  for (std::size_t i : keep) {
    const Convex_hull<GeomTraits>& ci = hulls[i];
    for (std::size_t j : keep) {
      if (j <= i)
        continue;
      const Convex_hull<GeomTraits>& cj = hulls[j];
      if (CGAL::do_intersect(ci.bbox, cj.bbox)) {
#ifdef CGAL_LINKED_WITH_TBB
        // Move this into a task list for later parallelization?
        todo.emplace_back(Merged_candidate(i, j));
#else
        Merged_candidate m(i, j);

        m.ch = num_hulls++;
        Convex_hull<GeomTraits>& ch = hulls[m.ch];
        ch.bbox = ci.bbox + cj.bbox;
        std::vector<Point_3> pts(ci.points.begin(), ci.points.end());
        pts.reserve(pts.size() + cj.points.size());
        std::copy(cj.points.begin(), cj.points.end(), std::back_inserter(pts));
        convex_hull(pts, ch.points, ch.indices);

        ch.volume = volume<GeomTraits>(ch.points, ch.indices);

        ch.volume_error = m.volume_error = CGAL::abs(ci.volume + cj.volume - ch.volume) / hull_volume;
        queue.push(std::move(m));
#endif
      }
      else {
        Merged_candidate m(i, j);
        Bbox_3 bbox = ci.bbox + cj.bbox;
        m.ch = -1;
        m.volume_error = CGAL::abs(ci.volume + cj.volume - bbox.x_span() * bbox.y_span() * bbox.z_span()) / hull_volume;
        queue.push(std::move(m));
      }
    }
  }

  // parallel for if available
#ifdef CGAL_LINKED_WITH_TBB
  tbb::parallel_for_each(todo, do_merge);
  for (Merged_candidate &m : todo)
    queue.push(std::move(m));
  todo.clear();
#endif

  while (!queue.empty() && keep.size() > max_convex_hulls) {
#ifdef CGAL_LINKED_WITH_TBB
    Merged_candidate m;
    while (!queue.try_pop(m) && !queue.empty());
#else
    Merged_candidate m = queue.top();
    queue.pop();
#endif

    auto ch_a = hulls.find(m.ch_a);
    if (ch_a == hulls.end())
      continue;

    auto ch_b = hulls.find(m.ch_b);
    if (ch_b == hulls.end())
      continue;

    if (m.ch == -1)
      do_merge(m);

    hulls.erase(ch_a);
    keep.erase(m.ch_a);
    hulls.erase(ch_b);
    keep.erase(m.ch_b);

    const Convex_hull<GeomTraits>& cj = hulls[m.ch];

    for (std::size_t id : keep) {
      const Convex_hull<GeomTraits>& ci = hulls[id];
      if (CGAL::do_intersect(ci.bbox, cj.bbox)) {
#ifdef CGAL_LINKED_WITH_TBB
        // Move this into a task list for later parallelization?
        todo.emplace_back(Merged_candidate(id, m.ch));
#else
        Merged_candidate merged(id, m.ch);

        merged.ch = num_hulls++;
        Convex_hull<GeomTraits>& ch = hulls[merged.ch];
        ch.bbox = ci.bbox + cj.bbox;
        std::vector<Point_3> pts(ci.points.begin(), ci.points.end());
        pts.reserve(pts.size() + cj.points.size());
        std::copy(cj.points.begin(), cj.points.end(), std::back_inserter(pts));
        convex_hull(pts, ch.points, ch.indices);

        ch.volume = volume<GeomTraits>(ch.points, ch.indices);

        ch.volume_error = merged.volume_error = CGAL::abs(ci.volume + cj.volume - ch.volume) / hull_volume;
        queue.push(std::move(merged));
#endif
      }
      else {
        Merged_candidate merged(id, m.ch);
        Bbox_3 bbox = ci.bbox + cj.bbox;
        merged.ch = -1;
        merged.volume_error = CGAL::abs(ci.volume + cj.volume - bbox.x_span() * bbox.y_span() * bbox.z_span()) / hull_volume;
        queue.push(std::move(merged));
      }
    }

    keep.insert(m.ch);

    std::cout << keep.size() << std::endl;

    // parallel for
#ifdef CGAL_LINKED_WITH_TBB
    tbb::parallel_for_each(todo, do_merge);
    for (Merged_candidate& m : todo)
      queue.push(std::move(m));
    todo.clear();
#endif
  }

  num_hulls = 0;

  for (std::size_t i : keep)
    candidates.push_back(std::move(hulls[i]));
}

}

template<typename FaceGraph, typename OutputIterator, typename NamedParameters = parameters::Default_named_parameters>
std::size_t approximate_convex_decomposition(const FaceGraph& mesh, std::size_t number_of_convex_hulls, OutputIterator out, const NamedParameters& np = parameters::default_values()) {
  CGAL::Memory_sizer memory;
  std::size_t virt_mem = memory.virtual_size();
  std::size_t res_mem = memory.resident_size();
  CGAL::Timer timer;

  using Geom_traits = typename GetGeomTraits<FaceGraph, NamedParameters>::type;
  using Vertex_point_map = typename GetVertexPointMap<FaceGraph, NamedParameters>::type;

  Vertex_point_map vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, mesh));

  using FT = typename Geom_traits::FT;

  // A single voxel grid should be sufficient
  // Why are they creating a voxel mesh and an AABB tree?
  // - just tracing the voxel grid should be sufficient? especially as the voxel grid
  // - tracing a voxel grid is very memory intensive, especially when it is sparse
  // Do they actually have a grid and fill it or only container of voxels?
  //const std::size_t voxel_size = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_voxels), 50);
  const std::size_t num_voxels = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_voxels), 1000000);
  const std::size_t max_depth = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_depth), 10);

  //std::vector<typename boost::graph_traits<FaceGraph>::face_descriptor> degenerate_faces;
  //CGAL::Polygon_mesh_processing::degenerate_faces(mesh, std::back_inserter(degenerate_faces));

  std::cout << (CGAL::is_closed(mesh) ? "input mesh is closed" : "input mesh is not closed") << std::endl;
  //std::cout << degenerate_faces.size() << " degenerate faces" << std::endl;

  Bbox_3 bb = bbox(mesh);
  //FT voxel_size;
  const auto [grid_size, voxel_size] = internal::calculate_grid_size<FT>(bb, num_voxels);

  std::cout << "grid_size: " << grid_size[0] << " " << grid_size[1] << " " << grid_size[2] << std::endl;

  timer.start();

  // if floodfill take INSIDE, otherwise UNKNOWN
  std::vector<int8_t> grid(grid_size[0] * grid_size[1] * grid_size[2], internal::Grid_cell::INSIDE);

  std::vector<internal::Candidate<Geom_traits>> candidates(1);

  init(candidates[0], mesh, grid, bb, grid_size, voxel_size);

  const FT hull_volume = candidates[0].ch.volume;

  double after_init = timer.time();

  recurse(candidates, grid, grid_size, bb, voxel_size, np);

  double after_recurse = timer.time();

//   for (std::size_t i = 0; i < candidates.size(); i++) {
//     CGAL::IO::write_polygon_soup(std::to_string(max_depth) + "-" + std::to_string(i) + "-d" + std::to_string(candidates[i].depth) + "-e" + std::to_string(candidates[i].ch.volume_error) + ".off", candidates[i].ch.points, candidates[i].ch.indices);
//   }

  std::size_t virt_mem2 = memory.virtual_size();
  std::size_t res_mem2 = memory.resident_size();

  std::vector<internal::Convex_hull<Geom_traits>> hulls;
  for (internal::Candidate<Geom_traits> &c : candidates)
    hulls.emplace_back(std::move(c.ch));

  candidates.clear();

  // merge until target number is reached
  merge(hulls, grid, grid_size, bb, voxel_size, hull_volume, np);

  double after_merge = timer.time();

  std::cout << (virt_mem2 - virt_mem) << " additional virtual memory allocated" << std::endl;
  std::cout << (res_mem2 - res_mem) << " additional resident memory occupied\n" << std::endl;

  std::cout << "timing:" << std::endl;
  std::cout << after_init << " initialization" << std::endl;
  std::cout << (after_recurse - after_init) << " recurse" << std::endl;
  std::cout << (after_merge - after_recurse) << " merge" << std::endl;

  std::cout << memory.virtual_size() << " virtual " << memory.resident_size() << " resident memory allocated in total" << std::endl;

  for (std::size_t i = 0;i<hulls.size();i++) {
    CGAL::IO::write_polygon_soup(std::to_string(i) + "-e" + std::to_string(hulls[i].volume_error) + ".off", hulls[i].points, hulls[i].indices);
  }

  return 0;
}
}
}

#endif