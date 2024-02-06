// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H
#define CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Grid_topology_3.h>

#include <CGAL/assertions.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Image_3.h>

#include <array>
#include <fstream>
#include <type_traits>
#include <vector>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domains_grp
 *
 * \brief stores scalar values and gradients at the vertices of a %Cartesian grid.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 */
template <typename GeomTraits>
class Cartesian_grid_3
{
public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;
  using Iso_cuboid_3 = typename Geom_traits::Iso_cuboid_3;

private:
  Iso_cuboid_3 m_bbox;
  std::array<FT,3> m_spacing;
  std::array<std::size_t, 3> m_sizes;

  std::vector<FT> m_values;
  std::vector<Vector_3> m_gradients;

  Geom_traits m_gt;

public:
  /**
   * \brief creates a grid with `xdim * ydim * zdim` grid vertices.
   *
   * The grid covers the space described by a bounding box.
   *
   * \param xdim the number of grid vertices in the `x` direction
   * \param ydim the number of grid vertices in the `y` direction
   * \param zdim the number of grid vertices in the `z` direction
   * \param bbox the bounding box of the grid
   * \param gt the geometric traits
   *
   * \pre `xdim`, `ydim`, and `zdim` are (strictly) positive.
   */
  Cartesian_grid_3(const std::size_t xdim,
                   const std::size_t ydim,
                   const std::size_t zdim,
                   const Iso_cuboid_3& bbox,
                   const Geom_traits& gt = Geom_traits())
    : m_sizes{xdim, ydim, zdim},
      m_bbox{bbox},
      m_gt{gt}
  {
    CGAL_precondition(xdim > 0);
    CGAL_precondition(ydim > 0);
    CGAL_precondition(zdim > 0);

    auto x_coord = m_gt.compute_x_3_object();
    auto y_coord = m_gt.compute_y_3_object();
    auto z_coord = m_gt.compute_z_3_object();
    auto vertex = m_gt.construct_vertex_3_object();

    // pre-allocate memory
    const std::size_t nv = xdim * ydim * zdim;
    m_values.resize(nv);
    m_gradients.resize(nv);

    // calculate grid spacing
    const Point_3& min_p = vertex(bbox, 0);
    const Point_3& max_p = vertex(bbox, 7);
    const FT x_span = x_coord(max_p) - x_coord(min_p);
    const FT y_span = y_coord(max_p) - y_coord(min_p);
    const FT z_span = z_coord(max_p) - z_coord(min_p);

    const FT d_x = x_span / (xdim - 1);
    const FT d_y = y_span / (ydim - 1);
    const FT d_z = z_span / (zdim - 1);
    m_spacing = make_array(d_x, d_y, d_z);
  }

  Cartesian_grid_3(const std::size_t xdim,
                   const std::size_t ydim,
                   const std::size_t zdim,
                   const Point_3& p, const Point_3& q,
                   const Geom_traits& gt = Geom_traits())
    : Cartesian_grid_3(xdim, ydim, zdim, Iso_cuboid_3(p, q), gt)
  { }

  /**
   * \brief creates a grid from a `CGAL::Image_3`.
   *
   * The dimensions and bounding box are read from the image. The values stored
   * in the image must be of type `Geom_traits::FT` or implicitly convertible to it.
   *
   * \param image the image providing the data
   */
  Cartesian_grid_3(const Image_3& image)
  {
    auto point = m_gt.construct_point_3_object();
    auto iso_cuboid = m_gt.construct_iso_cuboid_3_object();

    // compute bounding box
    const FT max_x = image.tx() + (image.xdim() - 1) * image.vx();
    const FT max_y = image.ty() + (image.ydim() - 1) * image.vy();
    const FT max_z = image.tz() + (image.zdim() - 1) * image.vz();
    m_bbox = iso_cuboid(point(image.tx(), image.ty(), image.tz()),
                        point(max_x, max_y, max_z));

    // get spacing
    m_spacing = make_array(image.vx(), image.vy(), image.vz());

    // get sizes
    m_sizes[0] = image.xdim();
    m_sizes[1] = image.ydim();
    m_sizes[2] = image.zdim();

    // pre-allocate
    const std::size_t nv = m_sizes[0] * m_sizes[1] * m_sizes[2];
    m_values.resize(nv);
    m_gradients.resize(nv);

    // copy values
    for(std::size_t x=0; x<m_sizes[0]; ++x)
      for(std::size_t y=0; y<m_sizes[1]; ++y)
        for(std::size_t z=0; z<m_sizes[2]; ++z)
          value(x, y, z) = image.value(x, y, z);
  }

   /**
  * \brief creates a `CGAL::Image_3` from the %Cartesian grid.
  */
  explicit operator Image_3() const;

public:
  /**
   * \return the geometric traits class
   */
  const Geom_traits& geom_traits() const
  {
    return m_gt;
  }

  /**
   * \return the number of grid vertices in the `x` direction
   */
  std::size_t xdim() const { return m_sizes[0]; }

  /**
   * \return the number of grid vertices in the `y` direction
   */
  std::size_t ydim() const { return m_sizes[1]; }

  /**
   * \return the number of grid vertices in the `z` direction
   */
  std::size_t zdim() const { return m_sizes[2]; }

  /**
   * \return the bounding box of the %Cartesian grid.
   */
  const Iso_cuboid_3& bbox() const { return m_bbox; }

  /**
   * \return the spacing of the %Cartesian grid, that is a vector whose coordinates are
   *         the grid steps in the `x`, `y`, and `z` directions, respectively
   */
  const std::array<FT, 3>& spacing() const { return m_spacing; }

public:
  /**
   * \brief gets the geometric position of the grid vertex described by a set of indices.
   *
   * Positions are not stored but calculated from an offset and grid spacing.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   *
   * \return the stored value
   *
   * \pre `i < xdim()` and `j < ydim()` and `k < zdim()`
   */
  Point_3 point(const std::size_t i,
                const std::size_t j,
                const std::size_t k) const
  {
    auto x_coord = m_gt.compute_x_3_object();
    auto y_coord = m_gt.compute_y_3_object();
    auto z_coord = m_gt.compute_z_3_object();
    auto point = m_gt.construct_point_3_object();
    auto vertex = m_gt.construct_vertex_3_object();

    const Point_3& min_p = vertex(m_bbox, 0);
    return point(x_coord(min_p) + i * m_spacing[0],
                 y_coord(min_p) + j * m_spacing[1],
                 z_coord(min_p) + k * m_spacing[2]);
  }

  /**
   * \brief gets the scalar value stored at the grid vertex described by a set of indices.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   *
   * \return the stored value
   *
   * \pre `i < xdim()` and `j < ydim()` and `k < zdim()`
   */
  FT value(const std::size_t i,
           const std::size_t j,
           const std::size_t k) const
  {
    return m_values[linear_index(i, j, k)];
  }

  /**
   * \brief gets the scalar value stored at the grid vertex described by a set of indices.
   *
   * \note This function can be used to set the value at a grid vertex.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   *
   * \return a reference to the stored value
   *
   * \pre `i < xdim()` and `j < ydim()` and `k < zdim()`
   */
  FT& value(const std::size_t i,
            const std::size_t j,
            const std::size_t k)
  {
    return m_values[linear_index(i, j, k)];
  }

public:
  /**
   * \brief gets the gradient stored at the grid vertex described by a set of indices.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   *
   * \pre `i < xdim()` and `j < ydim()` and `k < zdim()`
   */
  const Vector_3& gradient(const std::size_t i,
                           const std::size_t j,
                           const std::size_t k) const
  {
    return m_gradients[linear_index(i, j, k)];
  }

  /**
   * \brief gets the gradient stored at the grid vertex described by a set of indices.
   *
   * \note This function can be used to set the gradient at a grid vertex.
   *
   * \param i the index in the `x` direction
   * \param j the index in the `y` direction
   * \param k the index in the `z` direction
   *
   * \return a reference to the stored gradient
   *
   * \pre `i < xdim()` and `j < ydim()` and `k < zdim()`
   */
  Vector_3& gradient(const std::size_t i,
                     const std::size_t j,
                     const std::size_t k)
  {
    return m_gradients[linear_index(i, j, k)];
  }

private:
  std::size_t linear_index(const std::size_t i,
                           const std::size_t j,
                           const std::size_t k) const
  {
    CGAL_precondition(i < xdim() && j < ydim() && k < zdim());

    // convert (i, j, k) into a linear index to access the scalar values / gradient vectors
    return (k * ydim() + j) * xdim() + i;
  }
};

template <typename GeomTraits>
Cartesian_grid_3<GeomTraits>::
operator Image_3() const
{
  auto x_coord = m_gt.compute_x_3_object();
  auto y_coord = m_gt.compute_y_3_object();
  auto z_coord = m_gt.compute_z_3_object();
  auto vertex = m_gt.construct_vertex_3_object();

  // select number type
  WORD_KIND wordkind;
  if(std::is_floating_point<FT>::value) // @fixme seems meaningless given that vx vy vz are doubles
    wordkind = WK_FLOAT;
  else
    wordkind = WK_FIXED;

  // select signed or unsigned
  SIGN sign;
  if(std::is_signed<FT>::value)
    sign = SGN_SIGNED;
  else
    sign = SGN_UNSIGNED;

  // get spacing
  const double vx = CGAL::to_double(m_spacing[0]);
  const double vy = CGAL::to_double(m_spacing[1]);
  const double vz = CGAL::to_double(m_spacing[2]);

  // create image
  _image* im = _createImage(xdim(), ydim(), zdim(),
                            1,           // vectorial dimension
                            vx, vy, vz,  // voxel size
                            sizeof(FT),  // image word size in bytes
                            wordkind,    // image word kind WK_FIXED, WK_FLOAT, WK_UNKNOWN
                            sign);       // image word sign

  // error handling
  if(im == nullptr || im->data == nullptr)
    throw std::bad_alloc();  // @todo idk?

  // set min coordinates
  const Point_3& min_p = vertex(m_bbox, 0);
  im->tx = float(CGAL::to_double(x_coord(min_p)));
  im->ty = float(CGAL::to_double(y_coord(min_p)));
  im->tz = float(CGAL::to_double(z_coord(min_p)));

  // copy data
  FT* data = static_cast<FT*>(im->data); // @fixme what compatibility with non trivial FTs?
  for(std::size_t x=0; x<xdim(); ++x) {
    for(std::size_t y=0; y<ydim(); ++y) {
      for(std::size_t z=0; z<zdim(); ++z)
      {
        const std::size_t lid = linear_index(x, y, z);
        data[lid] = m_values[lid];
      }
    }
  }

  return Image_3{ im, Image_3::OWN_THE_DATA };
}

namespace IO {

template <typename GeomTraits,
          typename NamedParameters = parameters::Default_named_parameters>
bool write_OBJ(const std::string& filename,
               const Cartesian_grid_3<GeomTraits>& grid,
               const NamedParameters& np = parameters::default_values())
{
  using Point_3 = typename GeomTraits::Point_3;

  auto x_coord = grid.geom_traits().compute_x_3_object();
  auto y_coord = grid.geom_traits().compute_y_3_object();
  auto z_coord = grid.geom_traits().compute_z_3_object();
  auto vertex = grid.geom_traits().construct_vertex_3_object();

  std::ofstream out(filename);
  set_ascii_mode(out); // obj is ASCII only

  set_stream_precision_from_NP(out, np);

  if(out.fail())
    return false;


  // write vertices
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        const Point_3& p = vertex(grid.bbox(), 0);
        const double x_coord_d = CGAL::to_double(x_coord(p) + x * grid.spacing()[0]);
        const double y_coord_d = CGAL::to_double(y_coord(p) + y * grid.spacing()[1]);
        const double z_coord_d = CGAL::to_double(z_coord(p) + z * grid.spacing()[2]);
        out << "v " << x_coord_d << " " << y_coord_d << " " << z_coord_d << std::endl;
      }
    }
  }

  // write faces
  for(std::size_t x=0; x<grid.xdim()-1; ++x) {
    for(std::size_t y=0; y<grid.ydim()-1; ++y) {
      for(std::size_t z=0; z<grid.zdim()-1; ++z)
      {
        const std::size_t v0 = (z * grid.ydim() + y) * grid.xdim() + x;
        const std::size_t v1 = (z * grid.ydim() + y + 1) * grid.xdim() + x;
        const std::size_t v2 = (z * grid.ydim() + y + 1) * grid.xdim() + x + 1;
        const std::size_t v3 = (z * grid.ydim() + y) * grid.xdim() + x + 1;
        out << "f " << v0+1 << " " << v1+1 << " " << v2+1 << " " << v3+1 << std::endl;
      }
    }
  }

  return out.good();
}

} // namespace IO
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H
