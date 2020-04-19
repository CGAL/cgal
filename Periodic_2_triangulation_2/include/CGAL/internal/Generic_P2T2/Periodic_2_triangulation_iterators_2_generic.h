// Copyright (c) 2019-2020 XXX
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Georg Osang

#ifndef CGAL_PERIODIC_2_TRIANGULATION_ITERATORS_2_GENERIC_H
#define CGAL_PERIODIC_2_TRIANGULATION_ITERATORS_2_GENERIC_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/iterator.h>

namespace CGAL {

template < class T >
class Periodic_2_triangulation_point_iterator_2_generic
{
public:
  typedef typename T::Periodic_point                          Periodic_point;

  // typedef Periodic_2_triangulation_point_iterator_2<T>    Periodic_point_iterator;
  // typedef typename T::Periodic_point                      value_type;
  typedef const typename T::Periodic_point *                  pointer;
  typedef const typename T::Periodic_point &                  reference;

  // typedef typename T::Vertex_handle Vertex_handle;
  typedef typename T::DT2                                     DT2;
  typedef typename DT2::Vertex_iterator                       DT2_vertex_iterator;
  // typedef typename Canonicity_tester<T, DT2_vertex_iterator> > Tester;
  // typedef Filter_iterator<DT2_vertex_iterator, Tester> Filtered_DT2_vertex_iterator;

  typedef typename T::P2DT2                                   P2DT2;
  typedef typename P2DT2::Periodic_point_iterator             P2DT2_periodic_point_iterator;

  typedef Periodic_2_triangulation_point_iterator_2_generic   Self;

  /// Type determining how to iterate over the stored simplices in the triangulation
  typedef typename T::Iterator_type                           Iterator_type;

  // Periodic_2_triangulation_point_iterator_2_generic(Iterator_type it = T::STORED)
  //   : _t(nullptr), _it(it) {}

  // used to initialize the begin iterator
  Periodic_2_triangulation_point_iterator_2_generic(const T* t, Iterator_type it)
    : _t(t), _it(it)
  {
    switch(_it)
    {
      case Iterator_type::IT_P2T2:
        p2t2pos = _t->p2dt2.periodic_points_begin();
        break;
      case Iterator_type::IT_DT2:
        dt2pos = _t->dt2.vertices_begin();
        while(dt2pos != _t->dt2.vertices_end() && !t->is_canonical(dt2pos))
        {
          ++dt2pos;
        }
        break;
      default:
        CGAL_triangulation_assertion(false);
    };
  }

  // used to initialize the past-the-end iterator
  Periodic_2_triangulation_point_iterator_2_generic(const T* t, int, Iterator_type it)
    : _t(t), _it(it)
  {
    switch(_it)
    {
      case Iterator_type::IT_P2T2:
        p2t2pos = _t->p2dt2.periodic_points_end();
        break;
      case Iterator_type::IT_DT2:
        dt2pos = _t->dt2.vertices_end();
        break;
      default:
        CGAL_triangulation_assertion(false);
    };
  }

  Self& operator++()
  {
    switch(_it)
    {
      case Iterator_type::IT_P2T2:
        ++p2t2pos;
        break;
      case Iterator_type::IT_DT2:
        do
        {
          ++dt2pos;
        }
        while(dt2pos != _t->dt2.vertices_end() && !_t->is_canonical(dt2pos));
        break;
      default:
        CGAL_triangulation_assertion(false);
    };
    return *this;
  }

  Self& operator--()
  {
    switch(_it)
    {
      case Iterator_type::IT_P2T2:
        --p2t2pos;
        break;
      case Iterator_type::IT_DT2:
        do
        {
          --dt2pos;
        }
        while(dt2pos != _t->dt2.vertices_begin() && !_t->is_canonical(dt2pos));
        break;
      default:
        CGAL_triangulation_assertion(false);
    };
    return *this;
  }

  Self operator++(int)
  {
    Self tmp(*this);
    ++(*this);
    return tmp;
  }

  Self operator--(int)
  {
    Self tmp(*this);
    --(*this);
    return tmp;
  }

  bool operator==(const Self& pi) const
  {
    // We are only allowed to compare iterators of the same type.
    CGAL_triangulation_assertion(_it == pi._it);
    return _t == pi._t && (_it == Iterator_type::IT_DT2 ? dt2pos == pi.dt2pos : p2t2pos == pi.p2t2pos);
  }

  bool operator!=(const Self& pi) const
  {
    return !(*this == pi);
  }

  reference operator*() const
  {
    switch(_it)
    {
      case Iterator_type::IT_DT2:
        periodic_point = construct_periodic_point_dt2();
        break;
      case Iterator_type::IT_P2T2:
        periodic_point = *p2t2pos;
        break;
      default:
        CGAL_triangulation_assertion(false);
    };
    return periodic_point;
  }

  pointer operator->() const
  {
    switch(_it)
    {
      case Iterator_type::IT_DT2:
        periodic_point = construct_periodic_point_dt2();
        break;
      case Iterator_type::IT_P2T2:
        periodic_point = *p2t2pos;
        break;
      default:
        CGAL_triangulation_assertion(false);
    };
    return &periodic_point;
  }

private:
  Periodic_point construct_periodic_point_dt2() const
  {
    // CGAL_triangulation_assertion(dt2pos != typename T::DT2::Vertex_handle());
    return std::make_pair(dt2pos->point(), dt2pos->offset());
  }

  const T*  _t;
  DT2_vertex_iterator dt2pos; // current vertex.
  P2DT2_periodic_point_iterator p2t2pos; // current vertex.
  Iterator_type _it;
  mutable Periodic_point periodic_point; // current point.
};

} // namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_ITERATORS_2_GENERIC_H
