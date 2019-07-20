// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_CGAL_POLYGON_2_CURVE_ITERATOR_H
#define CGAL_CGAL_POLYGON_2_CURVE_ITERATOR_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <iterator>

namespace CGAL {



template <class Xcurve>
class Polygon_2_curve_ptr
{
public:
    typedef Xcurve X_monotnoe_curve_2;
    Polygon_2_curve_ptr(X_monotnoe_curve_2 const &c) :m_curve(c){}
    X_monotnoe_curve_2* operator->() {return &m_curve;}

private:
    X_monotnoe_curve_2 m_curve;
};

template <class X_monotone_curve_2_, class Polygon_>
class Polygon_2_curve_iterator 
{
public:
    typedef Polygon_2_curve_iterator< X_monotone_curve_2_, Polygon_ >    Self;
 
    
    typedef X_monotone_curve_2_                           X_monotone_curve_2;
    typedef typename Polygon_::Container                  Container;
    typedef typename Container::iterator                  Container_iterator;
    typedef typename std::iterator_traits<Container_iterator>::iterator_category 
                                                          iterator_category;
    typedef X_monotone_curve_2                            value_type;
    typedef X_monotone_curve_2*                           pointer;
    typedef X_monotone_curve_2&                           reference;


    typedef Polygon_                                      Polygon;
    typedef typename Polygon::Edge_const_iterator         Edge_const_iterator;
    typedef typename Edge_const_iterator::difference_type difference_type;
 
  private:
    const Polygon*      m_pgn;   // needed for dereferencing the last edge
    Edge_const_iterator m_curr_edge;   // points to the current edge iterator
  
public:
    Polygon_2_curve_iterator< X_monotone_curve_2_, Polygon_ >(){}

    Polygon_2_curve_iterator< X_monotone_curve_2_, Polygon_ >
      (const Polygon* pgn, Edge_const_iterator ci) : m_pgn(pgn),
                                                     m_curr_edge(ci) {}

    bool operator==(const Self& x) const
    {
      return (m_curr_edge == x.m_curr_edge);
    }
    
    bool operator!=(const Self& x) const
    {
      return !(m_curr_edge == x.m_curr_edge);
    }

    X_monotone_curve_2 operator*() 
    {
      return X_monotone_curve_2(*m_curr_edge);
    }
    
    Polygon_2_curve_ptr<X_monotone_curve_2> operator->()
    {
      return Polygon_2_curve_ptr<X_monotone_curve_2>(operator*());
    }

    Self& operator++()
    {
      ++m_curr_edge;
      return *this;
    }

    Self operator++(int)
    {
      Self tmp = *this;
      ++*this;
      return tmp;
    }

    Self& operator--() 
    {
      --m_curr_edge;
      return *this;
    }

    Self operator--(int) {
      Self tmp = *this;
      --*this;
      return tmp;
    }

    // random access iterator requirements
    Self& operator+=(difference_type n) 
    {
      m_curr_edge += n;
      return *this;
    }

    Self operator+(difference_type n)
    {
      return Self(m_pgn, m_curr_edge + n);
    }

    Self& operator-=(difference_type n) 
    {
      return (*this) -= n;
    }

    Self operator-(difference_type n) 
    {
      return Self(m_pgn, m_curr_edge - n);
    }

    difference_type operator-(const Self& a) const
    {
      return (const_cast<Edge_const_iterator&>(m_curr_edge) - a.m_curr_edge);
    }

    X_monotone_curve_2 operator[](int n) 
    {
      return *Self(m_pgn, m_curr_edge+n);
    }

    bool operator<(const Self& a)
    {
      return m_curr_edge < a.m_curr_edge;
    }

    bool operator>(const Self& a)
    {
      return m_curr_edge > a.m_curr_edge;
    }

    bool operator<=(const Self& a)
    {
      return m_curr_edge <= a.m_curr_edge;
    }

    bool operator>=(const Self& a)
    {
      return m_curr_edge >= a.m_curr_edge;
    }

};


template <class X_monotone_curve_2_,  class Polygon_>
typename Polygon_2_curve_iterator<X_monotone_curve_2_,Polygon_>::difference_type
distance_type(const Polygon_2_curve_iterator<X_monotone_curve_2_,Polygon_>&)
{ 
  return Polygon_2_curve_iterator<X_monotone_curve_2_,Polygon_>::difference_type(); 
}

template <class X_monotone_curve_2_,  class Polygon_>
X_monotone_curve_2_*
value_type(const Polygon_2_curve_iterator<X_monotone_curve_2_,Polygon_>&)
{
  return (X_monotone_curve_2_*)(0); 
}



} //namespace CGAL

#endif
