// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl
//                 Mael Rouxel-Labbé

#ifndef CGAL_ISOSURFACING_3_IO_IMAGE_3_H
#define CGAL_ISOSURFACING_3_IO_IMAGE_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>

#include <CGAL/Image_3.h>

namespace CGAL {
namespace Isosurfacing {
namespace IO {

/**
 * \ingroup IS_IO_functions_grp
 *
 * \brief creates a grid from a `CGAL::Image_3`.
 *
 * The dimensions and bounding box are read from the image. The values stored
 * in the image must be of type `Geom_traits::FT` or implicitly convertible to it.
 *
 * \tparam K must be a model of `IsosurfacingTraits_3`
 *
 * \param image the image providing the data
 * \param k the traits
 */
template <typename K>
std::pair<Cartesian_grid_3<K>,
          Interpolated_discrete_values_3<Cartesian_grid_3<K> > >
read_Image_3(const CGAL::Image_3& image,
             const K& k = K())
{
  using Grid = Cartesian_grid_3<K>;
  using Values = Interpolated_discrete_values_3<Grid>;

  using FT = typename K::FT;
  using Point_3 = typename K::Point_3;
  using Iso_cuboid_3 = typename K::Iso_cuboid_3;

  auto point = k.construct_point_3_object();
  auto iso_cuboid = k.construct_iso_cuboid_3_object();

  Iso_cuboid_3 bbox;

  // compute bounding box
  const FT max_x = image.tx() + (image.xdim() - 1) * image.vx();
  const FT max_y = image.ty() + (image.ydim() - 1) * image.vy();
  const FT max_z = image.tz() + (image.zdim() - 1) * image.vz();
  bbox = iso_cuboid(point(image.tx(), image.ty(), image.tz()),
                    point(max_x, max_y, max_z));

  // get spacing
  // std::array<FT, 3> spacing = make_array(image.vx(), image.vy(), image.vz());

  // get sizes
  std::array<std::size_t, 3> sizes;
  sizes[0] = image.xdim();
  sizes[1] = image.ydim();
  sizes[2] = image.zdim();

  Grid grid { bbox, sizes[0], sizes[1], sizes[2], k };

  Values values { grid };

  // copy values
  for(std::size_t x=0; x<sizes[0]; ++x)
    for(std::size_t y=0; y<sizes[1]; ++y)
      for(std::size_t z=0; z<sizes[2]; ++z)
        values(x, y, z) = image.value(x, y, z);

  return { grid, values };
}


/**
 * \ingroup IS_IO_functions_grp
 *
 * \brief create an `CGAL::Image_3` from a grid and a field of values.
 *
 * \tparam Grid must be `CGAL::Isosurfacing::Cartesian_grid_3<GeomTraits>` with `GeomTraits`
 *              a model of `IsosurfacingTraits_3`
 *
 * \param grid the space partitioning data structure
 * \param values the field of values
 */
template <typename Grid, typename Values>
CGAL::Image_3 write_Image_3(const Grid& grid,
                            const Values& values)
{
  using Geom_traits = typename Grid::Geom_traits;

  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;

  const Geom_traits& gt = grid.geom_traits();
  typename Geom_traits::Compute_x_3 x_coord = gt.compute_x_3_object();
  typename Geom_traits::Compute_y_3 y_coord = gt.compute_y_3_object();
  typename Geom_traits::Compute_z_3 z_coord = gt.compute_z_3_object();
  typename Geom_traits::Construct_vertex_3 vertex = gt.construct_vertex_3_object();

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
  const double vx = CGAL::to_double(grid.spacing()[0]);
  const double vy = CGAL::to_double(grid.spacing()[1]);
  const double vz = CGAL::to_double(grid.spacing()[2]);

  // create image
  _image* im = _createImage(grid.xdim(), grid.ydim(), grid.zdim(),
                            1,           // vectorial dimension
                            vx, vy, vz,  // voxel size
                            sizeof(FT),  // image word size in bytes
                            wordkind,    // image word kind WK_FIXED, WK_FLOAT, WK_UNKNOWN
                            sign);       // image word sign

  // error handling
  if(im == nullptr || im->data == nullptr)
    throw std::bad_alloc();  // @todo idk?

  // set min coordinates
  const Point_3& min_p = vertex(grid.bbox(), 0);
  im->tx = float(CGAL::to_double(x_coord(min_p)));
  im->ty = float(CGAL::to_double(y_coord(min_p)));
  im->tz = float(CGAL::to_double(z_coord(min_p)));

  // copy data
  FT* data = static_cast<FT*>(im->data); // @fixme what compatibility with non trivial FTs?
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        const std::size_t lid = grid.linear_index(x, y, z);
        data[lid] = values(grid.point(lid));
      }
    }
  }

  return Image_3 { im, Image_3::OWN_THE_DATA };
}

} // namespace IO
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_IO_IMAGE_3_H
