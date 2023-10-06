// Copyright (c) 20XX,20XX GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ange Clement

#ifndef CGAL_MESH_3_CONSTRUCT_INITIAL_POINTS_FROM_LABELED_IMAGE_H
#define CGAL_MESH_3_CONSTRUCT_INITIAL_POINTS_FROM_LABELED_IMAGE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/search_for_connected_components_in_labeled_image.h>

#include <CGAL/iterator.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Image_3.h>

template <typename Point>
struct Get_point
{
  const double vx, vy, vz;
  const double tx, ty, tz;
  const std::size_t xdim, ydim, zdim;
  Get_point(const CGAL::Image_3* image)
      : vx(image->vx())
      , vy(image->vy())
      , vz(image->vz())
      , tx(image->tx())
      , ty(image->ty())
      , tz(image->tz())
      , xdim(image->xdim())
      , ydim(image->ydim())
      , zdim(image->zdim())
  {}

  Point operator()(const std::size_t i,
                   const std::size_t j,
                   const std::size_t k) const
  {
    double x = double(i) * vx + tx;
    double y = double(j) * vy + ty;
    double z = double(k) * vz + tz;

    if (i == 0)              x += 1. / 6. * vx;
    else if (i == xdim - 1)  x -= 1. / 6. * vx;
    if (j == 0)              y += 1. / 6. * vy;
    else if (j == ydim - 1)  y -= 1. / 6. * vy;
    if (k == 0)              z += 1. / 6. * vz;
    else if (k == zdim - 1)  z -= 1. / 6. * vz;

    return Point(x, y, z);
  }
};

template <class MeshDomain>
struct Construct_initial_points_labeled_image
{
  const CGAL::Image_3 & image;

  Construct_initial_points_labeled_image(const CGAL::Image_3 & image_)
      : image(image_)
  { }

  template<class OutputIterator>
  OutputIterator operator()(OutputIterator pts, const int n = 20) const
  {
    MeshDomain domain = MeshDomain::create_labeled_image_mesh_domain(image);

    using ValueType = typename OutputIterator::container_type::value_type;

    using Point_3 = typename ValueType::first_type;
    using Index = typename ValueType::second_type;

    using R = typename Point_3::R;
    using Vector_3 = typename R::Vector_3;
    using Weighted_point = typename R::Weighted_point_3;
    using Segment_3 = typename R::Segment_3;
    using Image_word_type = unsigned char;

    const double max_v = (std::max)((std::max)(image.vx(),
                                               image.vy()),
                                               image.vz());

    struct Seed {
      std::size_t i, j, k;
      std::size_t radius;
    };
    using Seeds = std::vector<Seed>;

    Seeds seeds;
    Get_point<Point_3> get_point(&image);
    std::cout << "Searching for connected components..." << std::endl;
    search_for_connected_components_in_labeled_image(image,
                                                     std::back_inserter(seeds),
                                                     CGAL::Emptyset_iterator(),
                                                     CGAL::Identity<Image_word_type>(),
                                                     Image_word_type());
    std::cout << "  " << seeds.size() << " components were found." << std::endl;
    std::cout << "Construct initial points..." << std::endl;
    for(const Seed seed : seeds)
    {
      const Point_3 seed_point = get_point(seed.i, seed.j, seed.k);

      const double radius = double(seed.radius + 1)* max_v;
      CGAL::Random_points_on_sphere_3<Point_3> points_on_sphere_3(radius);
      typename MeshDomain::Construct_intersection construct_intersection =
          domain.construct_intersection_object();

      std::vector<Vector_3> directions;
      if(seed.radius < 2) {
        // shoot in six directions
        directions.push_back(Vector_3(-radius, 0, 0));
        directions.push_back(Vector_3(+radius, 0, 0));
        directions.push_back(Vector_3(0, -radius, 0));
        directions.push_back(Vector_3(0, +radius, 0));
        directions.push_back(Vector_3(0, 0, -radius));
        directions.push_back(Vector_3(0, 0, +radius));
      } else {
        for(int i = 0; i < n; ++i)
        {
          // shoot n random directions
          directions.push_back(*points_on_sphere_3++ - CGAL::ORIGIN);
        }
      }

      for(const Vector_3& v : directions)
      {
        const Point_3 test = seed_point + v;
        const Segment_3 test_segment = Segment_3(seed_point, test);

        const typename MeshDomain::Intersection intersect =
            construct_intersection(test_segment);
        if (std::get<2>(intersect) != 0)
        {
          const Point_3& bpi = std::get<0>(intersect);
          const Index index = std::get<1>(intersect);
          Weighted_point pi = Weighted_point(bpi);

          *pts++ = std::make_pair(bpi, index);
        }
      }
    }
    return pts;
  }
};

#endif // CGAL_MESH_3_CONSTRUCT_INITIAL_POINTS_FROM_LABELED_IMAGE_H
