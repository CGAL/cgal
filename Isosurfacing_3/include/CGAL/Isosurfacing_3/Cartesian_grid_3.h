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
#include <CGAL/Image_3.h>

#include <array>
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

private:
  Bbox_3 m_bbox;
  Vector_3 m_spacing;
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
   * \param bbox the bounding box
   * \param gt the geometric traits
   *
   * \pre `xdim`, `ydim`, and `zdim` are (strictly) positive.
   */
  Cartesian_grid_3(const std::size_t xdim,
                   const std::size_t ydim,
                   const std::size_t zdim,
                   const Bbox_3& bbox,
                   const Geom_traits& gt = Geom_traits())
    : m_sizes{xdim, ydim, zdim},
      m_bbox{bbox},
      m_gt{gt}
  {
    CGAL_precondition(xdim > 0);
    CGAL_precondition(ydim > 0);
    CGAL_precondition(zdim > 0);

    auto vector = m_gt.construct_vector_3_object();

    // pre-allocate memory
    const std::size_t nv = xdim * ydim * zdim;
    m_values.resize(nv);
    m_gradients.resize(nv);

    // calculate grid spacing
    const FT d_x = FT{bbox.x_span()} / (xdim - 1);
    const FT d_y = FT{bbox.y_span()} / (ydim - 1);
    const FT d_z = FT{bbox.z_span()} / (zdim - 1);
    m_spacing = vector(d_x, d_y, d_z);
  }

  /**
   * \brief creates a grid from an `Image_3`.
   *
   * The dimensions and bounding box are read from the image. The values stored
   * in the image must be of type `Geom_traits::FT` or implicitly convertible to it.
   *
   * \param image the image providing the data
   */
  Cartesian_grid_3(const Image_3& image)
  {
    from_image(image);
  }

  void from_image(const Image_3& image);
  Image_3 to_image() const;

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
  std::size_t xdim() const
  {
    return m_sizes[0];
  }

  /**
   * \return the number of grid vertices in the `y` direction
   */
  std::size_t ydim() const
  {
    return m_sizes[1];
  }

  /**
   * \return the number of grid vertices in the `z` direction
   */
  std::size_t zdim() const
  {
    return m_sizes[2];
  }

  /**
   * \return the bounding box of the %Cartesian grid.
   */
  const Bbox_3& bbox() const
  {
    return m_bbox;
  }

  /**
   * \return the spacing of the %Cartesian grid, that is a vector whose coordinates are
   *         the grid steps in the `x`, `y`, and `z` directions, respectively
   */
  const Vector_3& spacing() const
  {
    return m_spacing;
  }

  /**
   * \brief gets the geometric position of the grid vertex described by a set of indices.
   *
   * Positions are not stored but calculated from an offset and grid spacing.
   *
   * \param x the index in the `x` direction
   * \param y the index in the `y` direction
   * \param z the index in the `z` direction
   *
   * \return the stored value
   */
  const Point_3& point(const std::size_t x,
                       const std::size_t y,
                       const std::size_t z) const
  {
    typename Geom_traits::Compute_x_3 x_coord = m_gt.compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = m_gt.compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = m_gt.compute_z_3_object();
    typename Geom_traits::Construct_point_3 cp = m_gt.construct_point_3_object();

    return cp(m_bbox.xmin() + x * x_coord(m_spacing),
              m_bbox.ymin() + y * x_coord(m_spacing),
              m_bbox.zmin() + z * x_coord(m_spacing));
  }

  /**
   * \brief gets the scalar value stored at the grid vertex described by a set of indices.
   *
   * \param x the index in the `x` direction
   * \param y the index in the `y` direction
   * \param z the index in the `z` direction
   *
   * \return the stored value
   */
  FT value(const std::size_t x,
           const std::size_t y,
           const std::size_t z) const
  {
    return m_values[linear_index(x, y, z)];
  }

  /**
   * \brief gets the scalar value stored at the grid vertex described by a set of indices.
   *
   * \note This function can be used to set the value at a grid vertex.
   *
   * \param x the index in the `x` direction
   * \param y the index in the `y` direction
   * \param z the index in the `z` direction
   *
   * \return a reference to the stored value
   */
  FT& value(const std::size_t x,
            const std::size_t y,
            const std::size_t z)
  {
    return m_values[linear_index(x, y, z)];
  }

  /**
   * \brief gets the gradient stored at the grid vertex described by a set of indices.
   *
   * \param x the index in the `x` direction
   * \param y the index in the `y` direction
   * \param z the index in the `z` direction
   */
  const Vector_3& gradient(const std::size_t x,
                           const std::size_t y,
                           const std::size_t z) const
  {
    return m_gradients[linear_index(x, y, z)];
  }

  /**
   * \brief gets the gradient stored at the grid vertex described by a set of indices.
   *
   * \note This function can be used to set the gradient at a grid vertex.
   *
   * \param x the index in the `x` direction
   * \param y the index in the `y` direction
   * \param z the index in the `z` direction
   *
   * \return a reference to the stored gradient
   */
  Vector_3& gradient(const std::size_t x,
                     const std::size_t y,
                     const std::size_t z)
  {
    return m_gradients[linear_index(x, y, z)];
  }

private:
  std::size_t linear_index(const std::size_t x,
                           const std::size_t y,
                           const std::size_t z) const
  {
    CGAL_precondition(x < xdim() && y < ydim() && z < zdim());

    // convert (x, y, z) into a linear index to access the scalar value / gradient vectors
    return (z * ydim() + y) * xdim() + x;
  }
};

template <typename GeomTraits>
void
Cartesian_grid_3<GeomTraits>::
from_image(const Image_3& image)
{
  auto vector = m_gt.construct_vector_3_object();

  // compute bounding box
  const FT max_x = image.tx() + (image.xdim() - 1) * image.vx();
  const FT max_y = image.ty() + (image.ydim() - 1) * image.vy();
  const FT max_z = image.tz() + (image.zdim() - 1) * image.vz();
  m_bbox = Bbox_3{image.tx(), image.ty(), image.tz(), max_x, max_y, max_z};

  // get spacing
  m_spacing = vector(image.vx(), image.vy(), image.vz());

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

template <typename GeomTraits>
Image_3
Cartesian_grid_3<GeomTraits>::
to_image() const
{
  // select the number type
  WORD_KIND wordkind;
  if(std::is_floating_point<FT>::value)
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
  const double vx = m_spacing[0];
  const double vy = m_spacing[1];
  const double vz = m_spacing[2];

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
  im->tx = m_bbox.xmin();
  im->ty = m_bbox.ymin();
  im->tz = m_bbox.zmin();

  // copy data
  FT* data = static_cast<FT*>(im->data);
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

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H
