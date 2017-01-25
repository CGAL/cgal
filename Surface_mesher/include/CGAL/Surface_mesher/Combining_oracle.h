// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent Rineau


#ifndef CGAL_SURFACE_MESHER_COMBINING_ORACLE_H
#define CGAL_SURFACE_MESHER_COMBINING_ORACLE_H

#include <CGAL/license/Surface_mesher.h>


#include <list>
#include <algorithm>

#include <CGAL/Object.h>
#include <CGAL/Multi_surface_3.h>

#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>

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
    Oracle_a& oracle_a;
    Oracle_b& oracle_b;
  public:
    // Public types
    typedef typename Oracle_a::Point_3 Point_3;
    typedef typename Oracle_a::Segment_3 Segment_3;
    typedef typename Oracle_a::Ray_3 Ray_3;
    typedef typename Oracle_a::Line_3 Line_3;

    typedef typename Oracle_a::Intersection_point Intersection_point;

    CGAL_static_assertion((::boost::is_same<
                         Intersection_point,
                         typename Oracle_b::Intersection_point>::value));
                        

    typedef ::CGAL::Multi_surface_3<typename Oracle_a::Surface_3,
      typename Oracle_b::Surface_3> Surface_3;

    typedef typename Oracle_a::Intersect_3 Intersect_a;
    typedef typename Oracle_b::Intersect_3 Intersect_b;

  public:
    Combining_oracle(Oracle_a& oracle_a, Oracle_b& oracle_b)
      : oracle_a(oracle_a), oracle_b(oracle_b)
    {
    }

    class Intersect_3 
    {
      Oracle_a& oracle_a;
      Oracle_b& oracle_b;
    public:
      Intersect_3(Oracle_a& oracle_a, Oracle_b& oracle_b)
        : oracle_a(oracle_a), oracle_b(oracle_b)
      {
      }
      
      Object operator()(const Surface_3& surface, Segment_3 s) const
      {
        const Object obj = oracle_a.intersect_3_object()(surface.surface_a(), s);
        if( obj.is_empty() )
          return oracle_b.intersect_3_object()(surface.surface_b(), s);
        return obj;
      }

      Object operator()(const Surface_3& surface, const Ray_3& r) const {
        const Object obj = oracle_a.intersect_3_object()(surface.surface_a(), r);
        if( obj.is_empty() )
          return oracle_b.intersect_3_object()(surface.surface_b(), r);
        return obj;  
      }
      
      Object operator()(const Surface_3& surface, const Line_3& l) const {
        const Object obj = oracle_a.intersect_3_object()(surface.surface_a(), l);
        if( obj.is_empty() )
          return oracle_b.intersect_3_object()(surface.surface_b(), l);
        return obj;  
      }
    }; // end nested class Intersect_3

    class Construct_initial_points
    {
      Oracle_a& oracle_a;
      Oracle_b& oracle_b;
    public:
      Construct_initial_points(Oracle_a& oracle_a, Oracle_b& oracle_b)
        : oracle_a(oracle_a), oracle_b(oracle_b)
      {
      }

      // Random points
      template <typename OutputIteratorPoints>
      OutputIteratorPoints operator() (const Surface_3& surface, 
                                       OutputIteratorPoints out, 
                                       int n = 20) // WARNING: why 20?
      {
        OutputIteratorPoints out2 = 
          oracle_a.construct_initial_points_object()(surface.surface_a(),
                                                     out,
                                                     n);
        return oracle_b.construct_initial_points_object()(surface.surface_b(),
                                                          out2,
                                                          n);
      }
    }; // end nested class Construct_initial_points
     
    Intersect_3 intersect_3_object() const
    {
      return Intersect_3(oracle_a, oracle_b);
    }

    Construct_initial_points construct_initial_points_object() const
    {
      return Construct_initial_points(oracle_a, oracle_b);
    }

    bool is_in_volume(const Surface_3& surface, const Point_3& p) const
    {
      return( oracle_a.is_in_volume(surface.surface_a(), p) || 
              oracle_b.is_in_volume(surface.surface_b(), p) );
    }
  };  // end Combining_oracle

  }  // namespace Surface_mesher

} // namespace CGAL


#endif  // CGAL_SURFACE_MESHER_COMBINING_ORACLE_H
