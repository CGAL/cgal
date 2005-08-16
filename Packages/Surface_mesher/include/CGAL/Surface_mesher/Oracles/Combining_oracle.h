// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source: 
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_SURFACE_MESHER_COMBINING_ORACLE_H
#define CGAL_SURFACE_MESHER_COMBINING_ORACLE_H

#include <list>
#include <algorithm>

namespace CGAL {

  namespace Surface_mesher {

  /** Combine two oracles
      This class is a model of Oracle that combines two models of Oracle,
      so that the corresponding surface is the union of the two surfaces.
      \param Oracle_a First oracle type.
      \param Oracle_b Second oracle type.
  */
  template <class Oracle_a, class Oracle_b>
  class Combining_oracle
  {
  public:
    // Public types
    
    typedef typename Oracle_a::Point Point;
    typedef typename Oracle_a::Segment Segment;
    typedef typename Oracle_a::Ray Ray;
    typedef typename Oracle_a::Line Line;
    typedef typename Oracle_a::Triangle Triangle;

    typedef typename std::list<Point> Points;
  private:
    Oracle_a& oracle_a;
    Oracle_b& oracle_b;

  public:
    // Constructors
    Combining_oracle (Oracle_a& a, Oracle_b& b) : oracle_a(a), oracle_b(b) {}

    // Predicates and Constructions
    bool is_in_volume(const Point& p)
    {
      return( oracle_a.is_in_volume(p) || oracle_b.is_in_volume(p) );
    }

    Object intersect_segment_surface(Segment s) {      
      Object obj = oracle_a.intersect_segment_surface(s);
      if( obj.is_empty() )
	obj = oracle_b.intersect_segment_surface(s);
      return obj;
    }

    Object intersect_ray_surface(Ray r) {      
      Object obj = oracle_a.intersect_ray_surface(r);
      if( obj.is_empty() )
	obj = oracle_b.intersect_ray_surface(r);
      return obj;
    }
    
    Object intersect_line_surface(Line l) {      
      Object obj = oracle_a.intersect_line_surface(l);
      if( obj.is_empty() )
	obj = oracle_b.intersect_line_surface(l);
      return obj;
    }
    
    // Random points
    Points random_points (int n) {
      CGAL_precondition (n > 0);

      typename Oracle_a::Points points_a = oracle_a.random_points(n);
      typename Oracle_b::Points points_b = oracle_b.random_points(n);
      Points points;
      std::copy(points_a.begin(),
                points_a.end(),
                std::back_inserter(points));
       std::copy(points_b.begin(),
                 points_b.end(),
                 std::back_inserter(points));
      return points;
    }

  };  // end Combining_oracle


  }  // namespace Surface_mesher

} // namespace CGAL


#endif  // CGAL_SURFACE_MESHER_COMBINING_ORACLE_H
