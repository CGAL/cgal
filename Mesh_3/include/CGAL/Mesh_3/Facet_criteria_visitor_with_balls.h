// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
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
  typedef typename Base::Badness Facet_badness;
  typedef typename Base::Handle Handle;
  typedef Handle Facet;

  typedef typename Tr::Weighted_point Weighted_point;
  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::Compute_squared_radius_smallest_orthogonal_sphere_3 
  Squared_radius_orthogonal_sphere;
  int wp_nb_;
  double radius_ortho_shpere;
  double ratio;
  
  //typedef typename Tr::Cell::Surface_patch_index Surface_patch_index;
  //typedef typename Tr::Vertex_handle Vertex_handle;

   
  // Constructor
  Facet_criterion_visitor_with_balls(const Facet& fh)
    : Base(fh)
    , wp_nb_(0)
    , radius_ortho_shpere(0.)
    , ratio(0.)
  {
    Squared_radius_orthogonal_sphere sq_radius_ortho_sphere = 
      Gt().compute_squared_radius_smallest_orthogonal_sphere_3_object();
      
    Weighted_point p1 = fh.first->vertex ((fh.second+1)&3)->point();
    Weighted_point p2 = fh.first->vertex ((fh.second+2)&3)->point();
    Weighted_point p3 = fh.first->vertex ((fh.second+3)&3)->point();

    if(p1.weight() > 0) { ++wp_nb_; }
    if(p2.weight() > 0) { ++wp_nb_; }
    if(p3.weight() > 0) { ++wp_nb_; }
      
    switch ( wp_nb_ )
    {
    case 3:
      {
        radius_ortho_shpere = sq_radius_ortho_sphere(p1,p2,p3);
      }
      break;
          
    case 2:
      {
        if ( p3.weight() > 0 )
        { 
          if ( p1.weight() > 0 ) { std::swap(p2,p3); }
          else { std::swap(p1,p3); }  
        }
          
        radius_ortho_shpere = sq_radius_ortho_sphere(p1,p2);
          
        double f_size1 = CGAL::squared_distance(p1,p3);
        double f_size2 = CGAL::squared_distance(p2,p3);
          
        ratio = (f_size1 < f_size2) ? f_size1 / p1.weight()
          : f_size2 / p2.weight();
      }
      break;
          
    case 1:
      {
        if ( p2.weight() > 0 ) { std::swap(p1,p2); }
        else if ( p3.weight() > 0 ) { std::swap(p1,p3); }
          
        double f_size = (std::min)(CGAL::squared_distance(p1,p2),
                                   CGAL::squared_distance(p1,p3));
          
        ratio = f_size / p1.weight();
      }
      break;
          
    default: break;
    }
      
    //std::cerr << "radius: " << radius_ortho_shpere << "  ratio: " << ratio << "  wp_nb:" << wp_nb_ << "\n";
  }
   
  // Destructor
  ~Facet_criterion_visitor_with_balls() { };
   
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
