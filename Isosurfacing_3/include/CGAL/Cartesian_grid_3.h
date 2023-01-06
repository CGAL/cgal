// Copyright (c) 2022 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_CARTESIAN_GRID_3_H
#define CGAL_CARTESIAN_GRID_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/Isosurfacing_3/internal/Grid_topology.h>

#include <array>
#include <type_traits>
#include <vector>

namespace CGAL {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Stores scalar values and gradients at points of a cartesian grid.
 *
 * \tparam GeomTraits the traits type
 */
template <typename GeomTraits>
class Cartesian_grid_3
{
public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Vector = typename Geom_traits::Vector_3;

  using VertexDescriptor = Isosurfacing::internal::Grid_topology::Vertex_descriptor;

public:
  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \brief Create a grid with `xdim * ydim * zdim` grid points.
   * The grid covers the space described by a bounding box.
   *
   * \param xdim the number of grid points in x direction
   * \param ydim the number of grid points in y direction
   * \param zdim the number of grid points in z direction
   * \param bbox the bounding box
   */
  Cartesian_grid_3(const std::size_t xdim,
                   const std::size_t ydim,
                   const std::size_t zdim,
                   const Bbox_3& bbox)
    : sizes{xdim, ydim, zdim},
      bbox(bbox)
  {
    // pre-allocate memory
    values.resize(xdim * ydim * zdim);
    gradients.resize(xdim * ydim * zdim);

    // calculate grid spacing
    const FT d_x = bbox.x_span() / (xdim - 1);
    const FT d_y = bbox.y_span() / (ydim - 1);
    const FT d_z = bbox.z_span() / (zdim - 1);
    spacing = Vector(d_x, d_y, d_z);
  }

  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \brief Create a grid from an `Image_3`.
   * The dimensions and bounding box are read from the image.
   * The values stored in the image must be of type `Geom_traits::FT` or implicitly convertible to it.
   *
   * \param image the image providing the data
   */
  Cartesian_grid_3(const Image_3& image)
  {
    from_image(image);
  }

  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \brief Get the value stored in the grid at the grid point described by the vertex.
   *
   * \param v the vertex descriptor of the vertex
   *
   * \return the stored value
   */
  FT operator()(const VertexDescriptor& v) const
  {
    return values[linear_index(v[0], v[1], v[2])];
  }

  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \brief Get the value stored in the grid at the given indices.
   *
   * \param x the index in x direction
   * \param y the index in y direction
   * \param z the index in z direction
   *
   * \return the stored value
   */
  FT value(const std::size_t x,
           const std::size_t y,
           const std::size_t z) const
  {
      return values[linear_index(x, y, z)];
  }

  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \brief Get a reference to the value stored in the grid at the given indices.
   * Can be used to set values of the grid.
   *
   * \param x the index in x direction
   * \param y the index in y direction
   * \param z the index in z direction
   *
   * \return a reference to the stored value
   */
  FT& value(const std::size_t x,
            const std::size_t y,
            const std::size_t z)
  {
    return values[linear_index(x, y, z)];
  }

  Vector gradient(const std::size_t x,
                  const std::size_t y,
                  const std::size_t z) const
  {
    return gradients[linear_index(x, y, z)];
  }

  Vector& gradient(const std::size_t x,
                   const std::size_t y,
                   const std::size_t z)
  {
    return gradients[linear_index(x, y, z)];
  }

  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \return the number of grid points in x direction
   */
  std::size_t xdim() const
  {
    return sizes[0];
  }

  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \return the number of grid points in y direction
   */
  std::size_t ydim() const
  {
    return sizes[1];
  }

  /**
   * \ingroup PkgIsosurfacing3Ref
   *
   * \return the number of grid points in z direction
   */
  std::size_t zdim() const
  {
    return sizes[2];
  }

  const Bbox_3& get_bbox() const
  {
    return bbox;
  }

  const Vector& get_spacing() const
  {
    return spacing;
  }

private:
  std::size_t linear_index(const std::size_t x,
                           const std::size_t y,
                           const std::size_t z) const
  {
    // convert x, y, z into a linear index to access the vector
    return (z * ydim() + y) * xdim() + x;
  }

  void from_image(const Image_3& image);
  Image_3 to_image() const;

private:
  std::vector<FT> values;
  std::vector<Vector> gradients;

  std::array<std::size_t, 3> sizes;

  Bbox_3 bbox;
  Vector spacing;
};

template <typename GeomTraits>
void
Cartesian_grid_3<GeomTraits>::
from_image(const Image_3& image)
{
  // compute bounding box
  const FT max_x = image.tx() + (image.xdim() - 1) * image.vx();
  const FT max_y = image.ty() + (image.ydim() - 1) * image.vy();
  const FT max_z = image.tz() + (image.zdim() - 1) * image.vz();
  bbox = Bbox_3(image.tx(), image.ty(), image.tz(), max_x, max_y, max_z);

  // get spacing
  spacing = Vector(image.vx(), image.vy(), image.vz());

  // get sizes
  sizes[0] = image.xdim();
  sizes[1] = image.ydim();
  sizes[2] = image.zdim();

  // pre-allocate
  values.resize(xdim() * ydim() * zdim());
  gradients.resize(xdim() * ydim() * zdim());

  // copy values
  for(std::size_t x=0; x<sizes[0]; ++x)
    for(std::size_t y=0; y<sizes[1]; ++y)
      for(std::size_t z=0; z<sizes[2]; ++z)
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
  const double vx = spacing()[0];
  const double vy = spacing()[1];
  const double vz = spacing()[2];

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
  im->tx = bbox.xmin();
  im->ty = bbox.ymin();
  im->tz = bbox.zmin();

  // copy data
  FT* data = (FT*)im->data;
  for(std::size_t x=0; x<xdim(); ++x)
    for(std::size_t y=0; y<ydim(); ++y)
      for(std::size_t z=0; z<zdim(); ++z)
        data[(z * ydim() + y) * xdim() + x] = value(x, y, z);

  return Image_3(im, Image_3::OWN_THE_DATA);
}

} // namespace CGAL

#endif // CGAL_CARTESIAN_GRID_3_H
