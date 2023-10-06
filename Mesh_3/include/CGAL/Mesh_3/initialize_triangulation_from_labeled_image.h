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
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_INITIALIZE_TRIANGULATION_FROM_LABELED_IMAGE_H
#define CGAL_MESH_3_INITIALIZE_TRIANGULATION_FROM_LABELED_IMAGE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/make_mesh_3.h>

#include <CGAL/Mesh_3/Construct_initial_points_from_labeled_image.h>

template<class C3T3, class MeshDomain, class MeshCriteria>
void init_tr_from_labeled_image_call_init_features(C3T3&,
                                                   const MeshDomain&,
                                                   const MeshCriteria&,
                                                   CGAL::Tag_false)
{
}
template<class C3T3, class MeshDomain, class MeshCriteria>
void init_tr_from_labeled_image_call_init_features(C3T3& c3t3,
                                                   const MeshDomain& domain,
                                                   const MeshCriteria& criteria,
                                                   CGAL::Tag_true)
{
  CGAL::Mesh_3::internal::init_c3t3_with_features(c3t3,
                                                  domain,
                                                  criteria);
  std::cout << c3t3.triangulation().number_of_vertices()
            << " initial points on 1D-features" << std::endl;
}

template<class C3T3, class MeshDomain, class MeshCriteria,
         typename Image_word_type,
         typename TransformOperator = CGAL::Identity<Image_word_type> >
void initialize_triangulation_from_labeled_image(C3T3& c3t3,
      const MeshDomain&   domain,
      const CGAL::Image_3& image,
      const MeshCriteria& criteria,
      Image_word_type,
      bool protect_features = false,
      TransformOperator transform = CGAL::Identity<Image_word_type>())
{
  typedef typename MeshDomain::Point_3       Point_3;
  typedef typename MeshDomain::Index         Index;
  typedef typename MeshDomain::Subdomain     Subdomain;

  typedef typename std::pair<Point_3, Index> ConstructedPoint;

  typedef typename C3T3::Triangulation       Tr;
  typedef typename Tr::Vertex_handle         Vertex_handle;
  typedef typename Tr::Cell_handle           Cell_handle;
  typedef typename Tr::Weighted_point        Weighted_point;

  Tr& tr = c3t3.triangulation();

  typedef typename Tr::Geom_traits           GT;
  typedef typename GT::FT                    FT;
  typename GT::Compare_weighted_squared_radius_3 cwsr =
    tr.geom_traits().compare_weighted_squared_radius_3_object();
  typename GT::Construct_point_3 cp =
    tr.geom_traits().construct_point_3_object();
  typename GT::Construct_weighted_point_3 cwp =
    tr.geom_traits().construct_weighted_point_3_object();

  if(protect_features) {
    init_tr_from_labeled_image_call_init_features
      (c3t3, domain, criteria,
       CGAL::internal::Has_features<MeshDomain>());
  }

  std::vector<ConstructedPoint> constructedPoints;

  Construct_initial_points_labeled_image<MeshDomain> construct(image);
  construct(std::back_inserter(constructedPoints));

  std::cout << "  " << constructedPoints.size() << " constructed points." << std::endl;

  for (const ConstructedPoint & constructedPoint : constructedPoints)
  {
    const Point_3& point = constructedPoint.first;
    const Index&   index = constructedPoint.second;

    Weighted_point pi = cwp(point);

    Cell_handle seed_cell = tr.locate(cwp(point));

    const Subdomain seed_label
        = domain.is_in_domain_object()(point);
    const Subdomain seed_cell_label
      = (   tr.dimension() < 3
         || seed_cell == Cell_handle()
         || tr.is_infinite(seed_cell))
        ? Subdomain()  //seed_point is OUTSIDE_AFFINE_HULL
        : domain.is_in_domain_object()(
            seed_cell->weighted_circumcenter(tr.geom_traits()));

    if ( seed_label != std::nullopt
      && seed_cell_label != std::nullopt
      && *seed_label == *seed_cell_label)
        continue; //this means the connected component has already been initialized

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

      // if the (squared) distance between point and cv is smaller or equal than cv's weight
      if (cwsr(cv_wp, - tr.min_squared_distance(point, cp(cv_wp))) != CGAL::LARGER)
      {
        pi_inside_protecting_sphere = true;
        break;
      }
    }
    if (pi_inside_protecting_sphere)
      continue;

    /// The following lines show how to insert initial points in the
    /// `c3t3` object. [insert initial points]
    Vertex_handle v = tr.insert(pi);

    // `v` could be null if `pi` is hidden by other vertices of `tr`.
    CGAL_assertion(v != Vertex_handle());

    c3t3.set_dimension(v, 2); // by construction, points are on surface
    c3t3.set_index(v, index);
    /// [insert initial points]
  }


  if ( tr.dimension() != 3 )
  {
    std::cout << "  not enough points: triangulation.dimension() == "
              << tr.dimension() << std::endl;
    CGAL::Mesh_3::internal::init_c3t3(c3t3, domain, criteria, 20);
    std::cout << "  -> " << tr.number_of_vertices() << " initial points." << std::endl;
  }
}

#endif // CGAL_MESH_3_INITIALIZE_TRIANGULATION_FROM_LABELED_IMAGE_H
