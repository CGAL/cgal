// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
// 
//
// Author(s)     : Wieger Wesselink, Geert-Jan Giezeman <geert@cs.uu.nl>

#ifndef CGAL_POLYGON_2_EDGE_ITERATOR_H
#define CGAL_POLYGON_2_EDGE_ITERATOR_H

#include <iterator>
#include <CGAL/circulator.h>

namespace CGAL {
#ifndef DOXYGEN_RUNNING //to avoid conflicts
template <class Traits_, class Container_> class Polygon_2;
#endif
template <class Segment_>
class Polygon_2__Segment_ptr
{
public:
    typedef Segment_ Segment;
    Polygon_2__Segment_ptr(Segment const &seg) :m_seg(seg){}
    Segment* operator->() {return &m_seg;}
private:
    Segment m_seg;
};

template <class Traits_, class Container_>
class Polygon_2_edge_iterator {
  public:
    typedef typename std::iterator_traits<typename Container_::iterator>::iterator_category iterator_category;
    typedef typename Traits_::Segment_2 Segment_2;
    typedef typename Traits_::Segment_2 value_type;
    typedef Container_ Container;
    typedef typename Container_::const_iterator const_iterator;
    typedef typename Container_::difference_type difference_type;
    typedef Segment_2*           pointer;
    typedef Segment_2&           reference;
  private:
    const Container_* container;   // needed for dereferencing the last edge
    const_iterator first_vertex;   // points to the first vertex of the edge
  public:
    Polygon_2_edge_iterator() {}
    Polygon_2_edge_iterator(const Container_* c, const_iterator f)
      : container(c), first_vertex(f) {}

    bool operator==(
      const Polygon_2_edge_iterator<Traits_, Container_>& x) const
    {
      return first_vertex == x.first_vertex;
    }
    
    bool operator!=(
      const Polygon_2_edge_iterator<Traits_, Container_>& x) const
    {
      return !(first_vertex == x.first_vertex);
    }

    Segment_2 operator*() const {
      const_iterator second_vertex = first_vertex;
      ++second_vertex;
      if (second_vertex == container->end())
        second_vertex = container->begin();
      typename Traits_::Construct_segment_2 construct_segment_2 = 
            Traits_().construct_segment_2_object();
      return construct_segment_2(*first_vertex, *second_vertex);
    }
    
    Polygon_2__Segment_ptr<Segment_2> operator->() const
        {return Polygon_2__Segment_ptr<Segment_2>(operator*());}

    Polygon_2_edge_iterator<Traits_, Container_>& operator++() {
      ++first_vertex;
      return *this;
    }

    Polygon_2_edge_iterator<Traits_, Container_> operator++(int) {
      Polygon_2_edge_iterator<Traits_, Container_> tmp = *this;
      ++*this;
      return tmp;
    }

    Polygon_2_edge_iterator<Traits_, Container_>& operator--() {
      --first_vertex;
      return *this;
    }

    Polygon_2_edge_iterator<Traits_, Container_> operator--(int) {
      Polygon_2_edge_iterator<Traits_, Container_> tmp = *this;
      --*this;
      return tmp;
    }

// random access iterator requirements
    Polygon_2_edge_iterator<Traits_, Container_>&
    operator+=(difference_type n) {
      first_vertex += n;
      return *this;
    }

    Polygon_2_edge_iterator<Traits_, Container_>
    operator+(difference_type n) const {
      return Polygon_2_edge_iterator<Traits_, Container_>(
        container, first_vertex + n);
    }

    Polygon_2_edge_iterator<Traits_, Container_>&
    operator-=(difference_type n) {
      return (*this) -= n;
    }

    Polygon_2_edge_iterator<Traits_, Container_>
    operator-(difference_type n) const {
      return Polygon_2_edge_iterator<Traits_, Container_>(
        container, first_vertex - n);
    }

    difference_type
    operator-(const Polygon_2_edge_iterator<Traits_, Container_>& a) const {
      return first_vertex - a.first_vertex;
    }

    Segment_2 operator[](int n) const {
      return *Polygon_2_edge_iterator<Traits_, Container_>(
        container, first_vertex+n);
    }

    bool operator<(const Polygon_2_edge_iterator<Traits_, Container_>& a) const
    {
      return first_vertex < a.first_vertex;
    }

    bool operator>(const Polygon_2_edge_iterator<Traits_, Container_>& a) const
    {
      return first_vertex > a.first_vertex;
    }

    bool operator<=(const Polygon_2_edge_iterator<Traits_, Container_>& a) const
    {
      return first_vertex <= a.first_vertex;
    }

    bool operator>=(const Polygon_2_edge_iterator<Traits_, Container_>& a) const
    {
      return first_vertex >= a.first_vertex;
    }

};


template <class Traits_,  class Container_>
typename Container_::difference_type
distance_type(const Polygon_2_edge_iterator<Traits_,Container_>&)
{ return Container_::difference_type(); }

template <class Traits_,  class Container_>
typename Traits_::Segment_2*
value_type(const Polygon_2_edge_iterator<Traits_,Container_>&)
{ return (typename Traits_::Segment_2 *)(0); }


//-----------------------------------------------------------------------//
//                          implementation
//-----------------------------------------------------------------------//

//--------------------------------------------------------------------//
// I don't know how to implement the following function:
//
// template <class Traits_, class Container_>
// inline
// Polygon_2_edge_iterator<Traits_, Container_>
// operator+(Container_::difference_type n,
//           Polygon_2_edge_iterator<Traits_, Container_>& a)
// { return a+n; }
//--------------------------------------------------------------------//

} //namespace CGAL

#endif
