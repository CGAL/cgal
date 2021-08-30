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

#ifndef CGAL_MESH_3_FACET_CRITERIA_VISITOR_WITH_BALLS_H
#define CGAL_MESH_3_FACET_CRITERIA_VISITOR_WITH_BALLS_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/mesh_standard_criteria.h>
#include <CGAL/Mesh_3/mesh_standard_facet_criteria.h>
#include <CGAL/Mesh_3/Facet_on_same_surface_criterion.h>

namespace CGAL {
namespace Mesh_3 {

template <typename Tr>
class Facet_criterion_visitor_with_balls
  : public Mesh_3::Criterion_visitor<Tr, typename Tr::Facet>
{
  typedef Mesh_3::Criterion_visitor<Tr, typename Tr::Facet> Base;
  typedef Facet_criterion_visitor_with_balls<Tr> Self;

public:
  typedef Mesh_3::Abstract_criterion<Tr, Self> Criterion;
  typedef Mesh_3::Curvature_size_criterion<Tr, Self> Curvature_size_criterion;
  typedef Mesh_3::Aspect_ratio_criterion<Tr, Self> Aspect_ratio_criterion;
  typedef Mesh_3::Facet_on_surface_criterion<Tr, Self> Facet_on_surface_criterion;
  typedef Mesh_3::Uniform_size_criterion<Tr, Self> Uniform_size_criterion;
  typedef Mesh_3::Facet_on_same_surface_criterion<Tr, Self> Facet_on_same_surface_criterion;


  typedef typename Base::Quality Facet_quality;
  typedef typename Base::Is_bad  Is_facet_bad;
  typedef typename Base::Handle  Handle;
  typedef Handle Facet;

  typedef typename Tr::Bare_point      Bare_point;
  typedef typename Tr::Weighted_point  Weighted_point;
  typedef typename Tr::Geom_traits     Gt;
  typedef typename Gt::FT              FT;

  int wp_nb_;
  double radius_ortho_shpere;
  double ratio;

  //typedef typename Tr::Cell::Surface_patch_index Surface_patch_index;
  //typedef typename Tr::Vertex_handle Vertex_handle;


  // Constructor
  Facet_criterion_visitor_with_balls(const Tr& tr, const Facet& fh)
    : Base(tr, fh)
    , wp_nb_(0)
    , radius_ortho_shpere(0.)
    , ratio(0.)
  {
    typename Gt::Compare_weighted_squared_radius_3 compare_sq_radius =
      tr.geom_traits().compare_weighted_squared_radius_3_object();
    typename Gt::Compute_weight_3 cw =
      tr.geom_traits().compute_weight_3_object();
    typename Gt::Construct_point_3 cp =
      tr.geom_traits().construct_point_3_object();
    typename Gt::Squared_radius_orthogonal_sphere sq_radius_ortho_sphere =
      tr.geom_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object();

    Weighted_point wp1 = tr.point(fh.first, (fh.second+1)&3);
    Weighted_point wp2 = tr.point(fh.first, (fh.second+2)&3);
    Weighted_point wp3 = tr.point(fh.first, (fh.second+3)&3);

    if(compare_sq_radius(wp1, FT(0)) == CGAL::SMALLER) { ++wp_nb_; }
    if(compare_sq_radius(wp2, FT(0)) == CGAL::SMALLER) { ++wp_nb_; }
    if(compare_sq_radius(wp3, FT(0)) == CGAL::SMALLER) { ++wp_nb_; }

    switch ( wp_nb_ )
    {
    case 3:
      {
        radius_ortho_shpere = sq_radius_ortho_sphere(wp1, wp2, wp3);
      }
      break;

    case 2:
      {
        if(compare_sq_radius(wp3, FT(0)) == CGAL::SMALLER)
        {
          if(compare_sq_radius(wp1, FT(0)) == CGAL::SMALLER) { std::swap(wp2, wp3); }
          else { std::swap(wp1, wp3); }
        }

        radius_ortho_shpere = sq_radius_ortho_sphere(wp1, wp2);

        double f_size1 = CGAL::squared_distance(cp(wp1), cp(wp3));
        double f_size2 = CGAL::squared_distance(cp(wp2), cp(wp3));

        ratio = (f_size1 < f_size2) ? f_size1 / cw(wp1)
                                    : f_size2 / cw(wp2);
      }
      break;

    case 1:
      {
        if(compare_sq_radius(wp2, FT(0)) == CGAL::SMALLER) { std::swap(wp1, wp2); }
        else if(compare_sq_radius(wp3, FT(0)) == CGAL::SMALLER) { std::swap(wp1, wp3); }

        double f_size = (std::min)(CGAL::squared_distance(cp(wp1), cp(wp2)),
                                   CGAL::squared_distance(cp(wp1), cp(wp3)));

        ratio = f_size / cw(wp1);
      }
      break;

    default: break;
    }

    //std::cerr << "radius: " << radius_ortho_shpere << "  ratio: " << ratio << "  wp_nb:" << wp_nb_ << "\n";
  }

  // Destructor
  ~Facet_criterion_visitor_with_balls() { }

  void visit(const Criterion& criterion)
  {
    if ( wp_nb_ == 3 && radius_ortho_shpere <= 0.)
      Base::increment_counter();
    else
      Base::do_visit(criterion);
  }

  void visit(const Curvature_size_criterion& criterion)
  {
    if ( wp_nb_ >= 2 && radius_ortho_shpere <= 0.)
      Base::increment_counter();
    else if ( wp_nb_ == 1)
    {
      if ( ratio > 1.21 )
      {
        Base::do_visit(criterion);
        return;
      }
      else
        Base::increment_counter();
    }
    else
      Base::do_visit(criterion);

  }

  void visit(const Aspect_ratio_criterion& criterion)
  {
    if ( wp_nb_ >=2  && radius_ortho_shpere <= 0.)
    {
      if ( ratio > 4 )
      {
        Base::do_visit(criterion);
        return;
      }
      else
        Base::increment_counter();
    }
    else if ( wp_nb_ == 1)
      Base::increment_counter();
    else
      Base::do_visit(criterion);

  }

  //    void visit(const Facet_on_surface_criterion& criterion)
  //    {
  //    if ( wp_nb_ == 3 && radius_ortho_shpere <= 0.)
  //      Base::increment_counter();
  //    else
  //      Base::do_visit(criterion);
  //    }
  //
  //    void visit(const Uniform_size_criterion& criterion)
  //    {
  //        if ( wp_nb_ == 3 && radius_ortho_shpere <= 0.)
  //      Base::increment_counter();
  //    else
  //      Base::do_visit(criterion);
  //    }
  //
  //  void visit(const Facet_on_same_surface_criterion& criterion)
  //    {
  //    if ( wp_nb_ == 3 && radius_ortho_shpere <= 0.)
  //      Base::increment_counter();
  //    else
  //      Base::do_visit(criterion);
  //    }

};  // end class Facet_criterion_visitor


} //end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_FACET_CRITERIA_VISITOR_WITH_BALLS_H
