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
  typedef typename C3T3::Triangulation       Tr;
  typedef typename Tr::Geom_traits           GT;
  typedef typename Tr::Weighted_point        Weighted_point;
  typedef typename Tr::Vertex_handle         Vertex_handle;
  typedef typename MeshDomain::Point_3       Point_3;
  typedef typename MeshDomain::Index         Index;

  typedef typename std::pair<Point_3, Index> ConstructedPoint;

  typedef MeshDomain                         Mesh_domain;

  Tr& tr = c3t3.triangulation();

  typename GT::Construct_weighted_point_3 cwp =
    tr.geom_traits().construct_weighted_point_3_object();

  if(protect_features) {
    init_tr_from_labeled_image_call_init_features
      (c3t3, domain, criteria,
       CGAL::internal::Has_features<Mesh_domain>());
  }

  std::vector<ConstructedPoint> constructedPoints;

  Construct_initial_points_labeled_image construct(image);
  construct(std::back_inserter(constructedPoints), domain, c3t3);

  std::cout << "  " << constructedPoints.size() << " constructed points" << std::endl;

  for (const ConstructedPoint & constructedPoint : constructedPoints)
  {
    const Point_3& point = constructedPoint.first;
    const Index&   index = constructedPoint.second;

    Weighted_point pi = cwp(point);

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
  std::cout << "  " << tr.number_of_vertices() << " initial points." << std::endl;
  if ( c3t3.triangulation().dimension() != 3 )
  {
    std::cout << "  not enough points: triangulation.dimension() == "
              << c3t3.triangulation().dimension() << std::endl;
    CGAL::Mesh_3::internal::init_c3t3(c3t3, domain, criteria, 20);
    std::cout << "  -> " << tr.number_of_vertices() << " initial points." << std::endl;
  }
}

#endif // CGAL_MESH_3_INITIALIZE_TRIANGULATION_FROM_LABELED_IMAGE_H
