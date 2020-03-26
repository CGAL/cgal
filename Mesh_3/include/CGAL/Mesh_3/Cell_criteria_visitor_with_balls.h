// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************
#ifndef CGAL_MESH_3_CELL_CRITERIA_VISITOR_WITH_BALLS_H
#define CGAL_MESH_3_CELL_CRITERIA_VISITOR_WITH_BALLS_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/mesh_standard_criteria.h>
#include <CGAL/Mesh_3/mesh_standard_cell_criteria.h>

namespace CGAL {
namespace Mesh_3 {


template <typename Tr>
class Cell_criteria_visitor_with_balls
  : public Criterion_visitor<Tr, typename Tr::Cell_handle>
{
  typedef Criterion_visitor<Tr, typename Tr::Cell_handle> Base;
  typedef Cell_criteria_visitor_with_balls<Tr> Self;

public:
  typedef Abstract_criterion<Tr, Self> Criterion;
  typedef Mesh_3::Cell_radius_criterion<Tr, Self> Cell_radius_criterion;
  typedef Mesh_3::Cell_radius_edge_criterion<Tr, Self> Cell_radius_edge_criterion;

  typedef typename Base::Quality Cell_quality;
  typedef typename Base::Is_bad  Is_cell_bad;
  typedef typename Base::Handle  Handle;
  typedef Handle Cell_handle;

  typedef typename Tr::Bare_point      Bare_point;
  typedef typename Tr::Weighted_point  Weighted_point;
  typedef typename Tr::Geom_traits     Gt;
  typedef typename Gt::FT              FT;

  int nb_weighted_points;
  std::vector<Weighted_point> points;
  double radius_ortho_shpere;

  typedef typename Tr::Cell::Surface_patch_index Surface_patch_index;
  typedef typename Tr::Vertex_handle Vertex_handle;

  // Constructor
  Cell_criteria_visitor_with_balls(const Tr& tr, const Cell_handle& ch)
    : Base(tr, ch)
  {
    typename Gt::Compare_weighted_squared_radius_3 compare_sq_radius =
      tr.geom_traits().compare_weighted_squared_radius_3_object();
    typename Gt::Squared_radius_orthogonal_sphere sq_radius_ortho_sphere =
      tr.geom_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object();

    const Weighted_point& p = tr.point(ch, 0);
    const Weighted_point& q = tr.point(ch, 1);
    const Weighted_point& r = tr.point(ch, 2);
    const Weighted_point& s = tr.point(ch, 3);

    if(compare_sq_radius(p, FT(0)) == CGAL::SMALLER) points.push_back(p);
    if(compare_sq_radius(q, FT(0)) == CGAL::SMALLER) points.push_back(q);
    if(compare_sq_radius(r, FT(0)) == CGAL::SMALLER) points.push_back(r);
    if(compare_sq_radius(s, FT(0)) == CGAL::SMALLER) points.push_back(s);

    nb_weighted_points = (int)points.size();

    if(nb_weighted_points ==4)
      radius_ortho_shpere = sq_radius_ortho_sphere(points[0], points[1], points[2], points[3]);
    else if(nb_weighted_points ==3)
      radius_ortho_shpere = sq_radius_ortho_sphere(points[0], points[1], points[2]);
    else if(nb_weighted_points ==2)
      radius_ortho_shpere = sq_radius_ortho_sphere(points[0], points[1]);
  }

  // Destructor
  ~Cell_criteria_visitor_with_balls() { }

  //      void visit(const Criterion& criterion)
  //      {
  //              Base::do_visit(criterion);
  //      }

  void visit(const Cell_radius_criterion& criterion)
  {
    if ( nb_weighted_points == 4 && radius_ortho_shpere <= 0.)
      Base::increment_counter();
    else if(nb_weighted_points == 3 && radius_ortho_shpere <= 0. )
      Base::increment_counter();
    else if(nb_weighted_points == 2 && radius_ortho_shpere <= 0. )
      Base::increment_counter();
    else if(nb_weighted_points == 1)
      Base::do_visit(criterion);
    else
      Base::do_visit(criterion);

  }

  void visit(const Cell_radius_edge_criterion& criterion)
  {
    if ( nb_weighted_points == 4 && radius_ortho_shpere <= 0.)
      Base::increment_counter();
    else if(nb_weighted_points == 3 && radius_ortho_shpere <= 0. )
      Base::increment_counter();
    else if(nb_weighted_points == 2 && radius_ortho_shpere <= 0. )
      Base::increment_counter();
    else if(nb_weighted_points == 1)
      Base::increment_counter();
    else
      Base::do_visit(criterion);
  }

};  // end class Cell_criterion_visitor

} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_CELL_CRITERIA_VISITOR_WITH_BALLS_H
