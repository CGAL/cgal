// Copyright (c) 2015,2016 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau and Ange Clement

#ifndef CGAL_MESH_3_CONSTRUCT_INITIAL_POINTS_LABELED_IMAGE_H
#define CGAL_MESH_3_CONSTRUCT_INITIAL_POINTS_LABELED_IMAGE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/search_for_connected_components_in_labeled_image.h>

#include <CGAL/iterator.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Image_3.h>

namespace CGAL
{

namespace Mesh_3
{
namespace internal
{

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

} // end namespace internal
} // end namespace Mesh_3

/*!
 * \ingroup PkgMesh3Initializers
 *
 * Functor for generating initial points in labeled images.
 * This functor is a model of the `InitialPointsGenerator` concept,
 * and can be passed as a parameter to `CGAL::make_mesh_3` using the
 * `CGAL::parameters::initial_points_generator()` parameter function.
 *
 * On images that contain multiple non-connected components,
 * this functor scans the complete image and
 * outputs points to initialize each component.
 *
 * @tparam C3t3 model of `MeshComplex_3InTriangulation_3`
 * @tparam MeshDomain model of `MeshDomain_3`
 *
 * \cgalModels{InitialPointsGenerator}
 *
 * \sa `CGAL::parameters::initial_points_generator()`
 * \sa `CGAL::make_mesh_3()`
 * \sa `CGAL::Construct_initial_points_gray_image`
 */
template<typename C3t3, typename MeshDomain>
struct Construct_initial_points_labeled_image
{
  const CGAL::Image_3& image_;
  const MeshDomain& domain_;

  Construct_initial_points_labeled_image(const CGAL::Image_3& image,
                                         const MeshDomain& domain)
    : image_(image)
    , domain_(domain)
  { }

  /*!
  * \brief constructs at least `n` initial points,
  * by scanning the image and
  * collecting points in each object in the image,
  * even if they are non-connected components.
  * This ensures that each component will be initialized.
  *
  * @tparam OutputIterator model of `OutputIterator` for
  * tuple-like objects containing
  * - a `Weighted_point_3` for the point
  * - an `int` for the minimal dimension of the subcomplexes on which the point lies
  * - a `MeshDomain::Index` for the corresponding subcomplex index
  */
  template <typename OutputIterator>
  OutputIterator operator()(OutputIterator pts, int n = 20) const
  {
    CGAL_IMAGE_IO_CASE(image_.image(), operator()(pts, CGAL::Identity<Word>(), n));
    return pts;
  }

  /*!
   * \brief Same as above, but a `TransformOperator` that transforms values of the image is used.
   *
   * @tparam OutputIterator model of `OutputIterator` for
   * tuple-like objects containing
   * - a `Weighted_point_3` for the point
   * - an `int` for the minimal dimension of the subcomplexes on which the point lies
   * - a `MeshDomain::Index` for the corresponding subcomplex index
   * @tparam MeshDomain model of `MeshDomain_3`
   * @tparam TransformOperator functor that transforms values of the image.
   *   It must provide the following type:<br>
   *   `result_type` a type that supports the '==' operator<br>
   *   and the following operator:<br>
   *   `result_type operator()(Word v)`
   *   with `Word` the type of the image values.
   * @tparam C3t3 model of `MeshComplex_3InTriangulation_3`
   */
  template <typename OutputIterator, typename TransformOperator>
  OutputIterator operator()(OutputIterator pts, TransformOperator transform, int n = 20) const
  {
    typedef typename MeshDomain::Subdomain     Subdomain;
    typedef typename MeshDomain::Point_3       Point_3;
    typedef typename MeshDomain::Index         Index;

    typedef typename C3t3::Triangulation       Tr;
    typedef typename Tr::Geom_traits           GT;
    typedef typename GT::FT                    FT;
    typedef typename Tr::Weighted_point        Weighted_point;
    typedef typename Tr::Segment               Segment_3;
    typedef typename Tr::Vertex_handle         Vertex_handle;
    typedef typename Tr::Cell_handle           Cell_handle;
    typedef typename GT::Vector_3              Vector_3;

    C3t3 c3t3;
    Tr& tr = c3t3.triangulation();

    typename GT::Compare_weighted_squared_radius_3 cwsr =
      tr.geom_traits().compare_weighted_squared_radius_3_object();
    typename GT::Construct_point_3 cp =
      tr.geom_traits().construct_point_3_object();
    typename GT::Construct_weighted_point_3 cwp =
      tr.geom_traits().construct_weighted_point_3_object();

    const double max_v = (std::max)((std::max)(image_.vx(),
                                               image_.vy()),
                                               image_.vz());

    struct Seed {
      std::size_t i, j, k;
      std::size_t radius;
    };
    using Seeds = std::vector<Seed>;

    Seeds seeds;
    Mesh_3::internal::Get_point<Point_3> get_point(&image_);
    std::cout << "Searching for connected components..." << std::endl;
    CGAL_IMAGE_IO_CASE(image_.image(), search_for_connected_components_in_labeled_image(image_,
                                                     std::back_inserter(seeds),
                                                     CGAL::Emptyset_iterator(),
                                                     transform,
                                                     Word()));
    std::cout << "  " << seeds.size() << " components were found." << std::endl;
    std::cout << "Construct initial points..." << std::endl;
    for(const Seed seed : seeds)
    {
      const Point_3 seed_point = get_point(seed.i, seed.j, seed.k);
      Cell_handle seed_cell = tr.locate(cwp(seed_point));

      const Subdomain seed_label
        = domain_.is_in_domain_object()(seed_point);
      const Subdomain seed_cell_label
        = (   tr.dimension() < 3
           || seed_cell == Cell_handle()
           || tr.is_infinite(seed_cell))
          ? Subdomain()  //seed_point is OUTSIDE_AFFINE_HULL
          : domain_.is_in_domain_object()(
              seed_cell->weighted_circumcenter(tr.geom_traits()));

      if ( seed_label != std::nullopt
        && seed_cell_label != std::nullopt
        && *seed_label == *seed_cell_label)
          continue; //this means the connected component has already been initialized

      const double radius = double(seed.radius + 1)* max_v;
      CGAL::Random_points_on_sphere_3<Point_3> points_on_sphere_3(radius);
      // [construct intersection]
      typename MeshDomain::Construct_intersection construct_intersection =
          domain_.construct_intersection_object();
      // [construct intersection]

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
        const Segment_3 test_segment(seed_point, test);

        // [use construct intersection]
        const typename MeshDomain::Intersection intersect =
            construct_intersection(test_segment);
        // [use construct intersection]
        if (std::get<2>(intersect) != 0)
        {
          // [get construct intersection]
          const Point_3& intersect_point = std::get<0>(intersect);
          const Index& intersect_index = std::get<1>(intersect);
          // [get construct intersection]
          Weighted_point pi(intersect_point);

          // This would cause trouble to optimizers
          // check pi will not be hidden
          typename Tr::Locate_type lt;
          int li, lj;
          Cell_handle pi_cell = tr.locate(pi, lt, li, lj);
          if(lt != Tr::OUTSIDE_AFFINE_HULL) {
            switch (tr.dimension())
            { //skip dimension 0
            case 1:
              if (tr.side_of_power_segment(pi_cell, pi, true) != CGAL::ON_BOUNDED_SIDE)
                continue;
              break;
            case 2:
              if (tr.side_of_power_circle(pi_cell, 3, pi, true) != CGAL::ON_BOUNDED_SIDE)
                continue;
              break;
            case 3:
              if (tr.side_of_power_sphere(pi_cell, pi, true) != CGAL::ON_BOUNDED_SIDE)
                continue;
            }
          }

          //check pi is not inside a protecting ball
          std::vector<Vertex_handle> conflict_vertices;
          if (tr.dimension() == 3)
          {
            tr.vertices_on_conflict_zone_boundary(pi, pi_cell
              , std::back_inserter(conflict_vertices));
          }
          else
          {
            for (typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
                 vit != tr.finite_vertices_end(); ++vit)
            {
              const Weighted_point& wp = tr.point(vit);
              if (cwsr(wp, FT(0)) == CGAL::SMALLER) // 0 < wp's weight
                conflict_vertices.push_back(vit);
            }
          }

          bool pi_inside_protecting_sphere = false;
          for(Vertex_handle cv : conflict_vertices)
          {
            if(tr.is_infinite(cv))
              continue;

            const Weighted_point& cv_wp = tr.point(cv);
            if (cwsr(cv_wp, FT(0)) == CGAL::EQUAL) // 0 == wp's weight
              continue;

            // if the (squared) distance between intersect_point and cv is smaller or equal than cv's weight
            if (cwsr(cv_wp, - tr.min_squared_distance(intersect_point, cp(cv_wp))) != CGAL::LARGER)
            {
              pi_inside_protecting_sphere = true;
              break;
            }
          }
          if (pi_inside_protecting_sphere)
            continue;

          // insert point in temporary triangulation

          /// The following lines show how to insert initial points in the
          /// `c3t3` object. [insert initial points]
          Vertex_handle v = tr.insert(pi);

          // `v` could be null if `pi` is hidden by other vertices of `tr`.
          CGAL_assertion(v != Vertex_handle());

          c3t3.set_dimension(v, 2); // by construction, points are on surface
          c3t3.set_index(v, intersect_index);
          /// [insert initial points]

          *pts++ = std::make_pair(pi, intersect_index); // dimension 2 by construction, points are on surface
        }
      }
    }
    return pts;
  }
};

} // end namespace CGAL

#endif // CGAL_MESH_3_CONSTRUCT_INITIAL_POINTS_LABELED_IMAGE_H
