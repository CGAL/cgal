// Copyright (c) 2005 Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s): Baruch Zukerman         <baruchzu@post.tau.ac.il>
#ifndef EDGE_ITERATOR_ADAPTOR
#define EDGE_ITERATOR_ADAPTOR

#include <iterator>

CGAL_BEGIN_NAMESPACE

template <class Base_X_monotone_curve_2,
          class X_monotone_curve_2,
          class Arr_A,
          class Arr_B> 
class Edge_iterator_adaptor
{
public:

  typedef Edge_iterator_adaptor<Base_X_monotone_curve_2,
                                X_monotone_curve_2,
                                Arr_A,
                                Arr_B>           Self;


  typedef std::forward_iterator_tag              iterator_category;
  typedef X_monotone_curve_2                     value_type;
  typedef X_monotone_curve_2*                    pointer;
  typedef X_monotone_curve_2&                    reference;
 

protected:

  typedef typename Arr_A::Edge_const_iterator      Edge_iterator_A;
  typedef typename Arr_B::Edge_const_iterator      Edge_iterator_B;

  typedef typename Arr_A::Halfedge_const_iterator  Halfedge_iterator_A;
  typedef typename Arr_B::Halfedge_const_iterator  Halfedge_iterator_B;

  Edge_iterator_A    m_a_iter;
  Edge_iterator_B    m_b_iter;

  bool               m_is_b;

  Edge_iterator_A    m_a_end;
  Edge_iterator_B    m_b_begin;

  

public:
  typedef typename Edge_iterator_A::difference_type    difference_type;

  Edge_iterator_adaptor<Base_X_monotone_curve_2,
                        X_monotone_curve_2,
                        Arr_A,
                        Arr_B>() {}

  Edge_iterator_adaptor<Base_X_monotone_curve_2,
                        X_monotone_curve_2,
                        Arr_A,
                        Arr_B>(Edge_iterator_A curr,
                               Edge_iterator_A a_end,
                               Edge_iterator_B b_begin):
    m_b_iter(b_begin),
    m_a_end(a_end),
    m_b_begin(b_begin)
  {
    if(curr != a_end)
    {
      m_a_iter = curr;
      m_is_b = false;
    }
    else // first range is empty
    {
      m_b_iter = b_begin;
      m_is_b = true;
    }
  }

  Edge_iterator_adaptor<Base_X_monotone_curve_2,
                        X_monotone_curve_2,
                        Arr_A,
                        Arr_B>(Edge_iterator_B curr):
    m_b_iter(curr),
    m_is_b(true)
  {}

   bool operator==(const Self& x) const
    {
      if(m_is_b != x.m_is_b)
        return false;

      if(!m_is_b)
        return (m_a_iter == x.m_a_iter);

      return (m_b_iter == x.m_b_iter);
    }
    
    bool operator!=(const Self& x) const
    {
      return !(this->operator==(x));
    }

    X_monotone_curve_2 operator*() 
    {
      if(!m_is_b)
      {
        Halfedge_iterator_A he = m_a_iter;
        if (he->direction() == SMALLER)
          he = he->twin();
        
        return (X_monotone_curve_2 (he->curve(), he, Halfedge_iterator_B()));
      }
      else
      {
        Halfedge_iterator_B he = m_b_iter;
        if (he->direction() == SMALLER)
          he = he->twin();
        
        return (X_monotone_curve_2 (he->curve(), Halfedge_iterator_A(), he));
      }
     
    }
    
    Self& operator++()
    {
      if(!m_is_b)
      {
        ++m_a_iter;
        if(m_a_iter == m_a_end)
          m_is_b = true;
      }
      else
        ++m_b_iter;
      
      return (*this);
    }


    Self operator++(int)
    {
      Self tmp = *this;
      ++*this;
      return tmp;
    }

};

CGAL_END_NAMESPACE

#endif
