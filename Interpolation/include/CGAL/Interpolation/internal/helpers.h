// Copyright (c) 2003,2017  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Julia Floetotto

#ifndef CGAL_INTERPOLATION_INTERNAL_HELPERS_H
#define CGAL_INTERPOLATION_INTERNAL_HELPERS_H

#include <CGAL/license/Interpolation.h>


namespace CGAL {
namespace Interpolation {
namespace internal {
  
  template < class InterpolationTraits>
  struct V2P
  {    
    typedef typename InterpolationTraits::Point_d Point;
    typedef typename InterpolationTraits::Weighted_point_d Weighted_point;

    V2P(const InterpolationTraits& traits)
      : traits(traits)
    {}
    
    template <typename VH>
    const Point& operator()(const VH& vh) const
    {
      return traits.construct_point_d_object()(vh->point());
    }
   
    const Point& operator()(const Point& p) const
    {
      return p;
    }
    
    Point operator()(const Weighted_point& wp) const
    {
      return traits.construct_point_d_object()(wp);
    }

  private:
    InterpolationTraits traits;
  };


  template < typename Dt, typename T2>
  struct Vertex2Point {
    typedef typename Dt::Vertex_handle Vertex_handle;
    typedef typename Dt::Point Point;
    
    typedef std::pair<Vertex_handle, T2> argument_type;
    typedef std::pair<Point, T2> result_type;
    
    result_type operator()(const argument_type& vp) const
    {
      return std::make_pair(vp.first->point(), vp.second);
    }
  };

  
  template <typename Dt, typename Map>
  struct Vertex2Vertex {
    typedef typename Dt::Vertex_handle Vertex_handle;
    typedef typename Dt::Geom_traits::FT FT;
    typedef std::pair<Vertex_handle,FT> argument_type;
    typedef std::pair<Vertex_handle,FT> result_type;
    
    const Map& map;
    const Dt& dt;

    Vertex2Vertex(const Map& map, const Dt& dt)
      : map(map), dt(dt)
    {}

    result_type operator()(const argument_type& vp) const
    {
      typename Map::const_iterator it = map.find(vp.first);
      CGAL_assertion(it != map.end());
      CGAL_assertion(dt.tds().is_vertex(it->second));
      return std::make_pair(it->second, vp.second);
    }
  };

  
  template <typename Tr, typename InterpolationTraits, typename FunctorArgType, typename T2>
  struct Output_iterator_functor_selector
  {
    typedef CGAL::Identity<std::pair<FunctorArgType, T2> > type;
  };

  
  template <typename Tr, typename InterpolationTraits, typename T2>
  struct Output_iterator_functor_selector<Tr, InterpolationTraits,
                                          typename InterpolationTraits::Point_d, T2>
  {
    typedef Interpolation::internal::Vertex2Point<Tr, T2> type;
  };

  
  // the struct "Project_vertex_output_iterator"
  // is used in the (next two) functions
  // as well as in regular_neighbor_coordinates_2 and
  // in surface_neighbor_coordinates_3
  //
  //projection of iterator over std::pair <Vertex_handle, T>
  //to iterator over std::pair< Point, T>
  template < class OutputIterator, class Fct = void>
  struct Project_vertex_output_iterator
  {
    // this class wraps the OutputIterator with value type
    // std::pair<Vertex_handle,T>
    // into an output iterator with value type std::pair<Point, T>
    // Conditions: OutputIterator has value type std::pair<Vertex_handle, T>
    //             and Vertex_handle has a function ->point()
    //             with return type const Point&

    OutputIterator _base;
    Fct fct;

    //creation:
    Project_vertex_output_iterator(OutputIterator o, Fct fct)
      : _base(o), fct(fct)
    {}

    OutputIterator base() {return _base;}

    Project_vertex_output_iterator& operator++(){_base++; return *this;}
    Project_vertex_output_iterator& operator++(int){_base++; return *this;}
    Project_vertex_output_iterator& operator*(){return *this;}

    template<class Vertex_pair>
    Project_vertex_output_iterator&
    operator=(const Vertex_pair& vp){
      *_base = fct(vp);
      return *this;
    }
  };

  
} // namespace internal
} // namespace Interpolation
} // namespace CGAL

#endif // CGAL_INTERPOLATION_INTERNAL_HELPERS_H
