// Copyright (c) 2003, 2017 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julia Floetotto

#ifndef CGAL_INTERPOLATION_INTERNAL_HELPERS_H
#define CGAL_INTERPOLATION_INTERNAL_HELPERS_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/assertions.h>

#include <utility>

namespace CGAL {

namespace Interpolation {

namespace internal {

// Extracts the bare point from a weighted point or a vertex handle (and identity
// if the argument is a bare point)
template < typename InterpolationTraits_ >
struct Extract_bare_point
{
  typedef InterpolationTraits_                       Traits;
  typedef typename Traits::Point_d                   Point;
  typedef typename Traits::Weighted_point_d          Weighted_point;

  Extract_bare_point(const Traits& traits = Traits()) : traits(traits) {}

  const Point& operator()(const Point& p) const { return p; }

  Point operator()(const Weighted_point& wp) const {
    return traits.construct_point_d_object()(wp);
  }

  template <typename VH>
  Point operator()(const VH& vh) const {
    CGAL_precondition(vh != VH());
    return traits.construct_point_d_object()(vh->point());
  }

private:
  Traits traits;
};


// Converts a pair<Triangulation_::Vertex_handle, T2_>
// into     a pair<Triangulation_::Point, T2_>
template < typename Triangulation_, typename T2_ >
struct Extract_point_in_pair
{
  typedef Triangulation_                              Tr;
  typedef T2_                                         T2;
  typedef typename Tr::Vertex_handle                  Vertex_handle;
  typedef typename Tr::Point                          Point; // possibly weighted

  typedef std::pair<Vertex_handle, T2>                argument_type;
  typedef std::pair<Point, T2>                        result_type;

  result_type operator()(const argument_type& vp) const
  {
    CGAL_precondition(vp.first != Vertex_handle());
    return std::make_pair(vp.first->point(), vp.second);
  }
};


// Given a map `m` of type `Map`, transforms an object of type `pair<Map::key_type, T2>`
// into an object of type `pair<Map::mapped_type, T2>`, using the `m`.
//
template < typename Map_, typename T2_ >
struct Pair_mapper
{
  typedef Map_                                       Map;
  typedef T2_                                        T2;

  typedef typename Map::key_type                     Map_key_type;
  typedef typename Map::mapped_type                  Map_mapped_type;

  typedef std::pair<Map_key_type, T2>                argument_type;
  typedef std::pair<Map_mapped_type, T2>             result_type;

  const Map& map;

  Pair_mapper(const Map& map) : map(map) {}

  result_type operator()(const argument_type& vp) const
  {
    typename Map::const_iterator it = map.find(vp.first);
    CGAL_assertion(it != map.end());

    return std::make_pair(it->second, vp.second);
  }
};


// The struct "Project_vertex_output_iterator" is used in natural_neighbor_coordinates_2,
// as well as in regular_neighbor_coordinates_2, and in surface_neighbor_coordinates_3.
//
// It wraps the OutputIterator with value type `std::pair<Vertex_handle, T2>`
// into an output iterator with value type `OutputFunctor::result_type`.
//
// \pre OutputIterator has value type std::pair<Vertex_handle, T>
// \pre Vertex_handle has a function ->point() with return type const Point&
// \pre OutputIterator::value_type == OutputFunctor::argument_type
template < class OutputIterator, class OutputFunctor >
struct Project_vertex_output_iterator
{
  OutputIterator out;
  OutputFunctor fct;

  //creation:
  Project_vertex_output_iterator(OutputIterator o, OutputFunctor f)
    : out(o), fct(f)
  {}

  OutputIterator base() {return out;}

  Project_vertex_output_iterator& operator++(){out++; return *this;}
  Project_vertex_output_iterator& operator++(int){out++; return *this;}
  Project_vertex_output_iterator& operator*(){return *this;}

  template<class Vertex_pair>
  Project_vertex_output_iterator& operator=(const Vertex_pair& vp)
  {
    *out = fct(vp);
    return *this;
  }
};

} // namespace internal
} // namespace Interpolation
} // namespace CGAL

#endif // CGAL_INTERPOLATION_INTERNAL_HELPERS_H
