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
  typedef typename Base::Badness Cell_badness;
  typedef typename Base::Handle Handle;
  typedef Handle Cell_handle;

  typedef typename Tr::Point Point_3;
  typedef typename Tr::Geom_traits::Compute_squared_radius_smallest_orthogonal_sphere_3 
  Squared_radius_orthogonal_sphere;
  int nb_weighted_points;
  std::vector<Point_3> points;
  double radius_ortho_shpere;
        
  typedef typename Tr::Cell::Surface_patch_index Surface_patch_index;
  typedef typename Tr::Vertex_handle Vertex_handle;

  // Constructor
  Cell_criteria_visitor_with_balls(const Cell_handle& ch)
    : Base(ch) 
  {

    const Point_3& p = ch->vertex(0)->point();
    const Point_3& q = ch->vertex(1)->point();
    const Point_3& r = ch->vertex(2)->point();
    const Point_3& s = ch->vertex(3)->point();

    if(p.weight() > 0) points.push_back(p);
    if(q.weight() > 0) points.push_back(q);
    if(r.weight() > 0) points.push_back(r);
    if(s.weight() > 0) points.push_back(s);

    nb_weighted_points = (int)points.size();
                
    if(nb_weighted_points ==4)
      radius_ortho_shpere = Squared_radius_orthogonal_sphere()( points[0], points[1], points[2], points[3]);
    else if(nb_weighted_points ==3)
      radius_ortho_shpere = Squared_radius_orthogonal_sphere()( points[0], points[1], points[2]);
    else if(nb_weighted_points ==2)
      radius_ortho_shpere = Squared_radius_orthogonal_sphere()( points[0], points[1]);
                
  }

  // Destructor
  ~Cell_criteria_visitor_with_balls() { };
        
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
