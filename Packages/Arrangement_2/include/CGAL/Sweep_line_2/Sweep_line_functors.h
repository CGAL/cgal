// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>,
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_SWEEP_LINE_FUNCTORS_H
#define CGAL_SWEEP_LINE_FUNCTORS_H

#include <CGAL/assertions.h>
#include <CGAL/Sweep_line_2/Sweep_line_traits.h>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

template <class SweepLineTraits_2>
class Point_less_functor 
{
public:
  typedef SweepLineTraits_2           Traits;
  typedef typename Traits::Point_2    Point_2;
  
  Point_less_functor(Traits * t) : m_traits(t)
  {}
  
  bool operator()(const Point_2& p1,const Point_2& p2) const  
  { 
    return (m_traits->compare_xy_2_object()(p1,p2) == SMALLER);
  }

private:

  /*! a pointer to a traits object */
  Traits * m_traits;
};




template <class SweepLineTraits_2, class Subcurve> 
class Status_line_curve_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  

  Status_line_curve_less_functor(Traits *t) : m_traits(t) 
  {}

  Comparison_result operator()(const Subcurve * c1, const Subcurve * c2) const 
  {
    return m_traits->compare_y_at_x_2_object() (c1->get_left_end(),
                                           c2->get_last_curve());
  }

  Comparison_result operator()(const Point_2& pt, const Subcurve * c2) const 
  {
    return m_traits->compare_y_at_x_2_object()(pt,c2->get_last_curve());
  }

     

private:

  /*! a pointer to a traits object */
  Traits * m_traits;
};


template <class _Subcurve>
class Curves_pair
{
  public:
  Curves_pair(){}

  Curves_pair(_Subcurve* sc1, _Subcurve* sc2)
  {
    //the smallest pointer will be the first 
    if(sc1 < sc2)
      m_curves_pair = std::make_pair(sc1,sc2);
    else
      m_curves_pair = std::make_pair(sc2,sc1);
  }

  _Subcurve* first() const {return m_curves_pair.first;}
  _Subcurve* second() const {return m_curves_pair.second;}

 private:
  std::pair<_Subcurve*, _Subcurve*> m_curves_pair;
};

template <class _Subcurve>
class Curves_pair_less_functor
{
  typedef Curves_pair<_Subcurve> CurvesPair;

public:

  bool operator()(const CurvesPair& pair1, const CurvesPair& pair2)
  {
    if(pair1.first() < pair2.first())
      return true;
    if(pair1.first() > pair2.first())
      return false;
    if(pair1.second() < pair2.second())
      return true;
    return false;
  }
};


template <class Container>
class random_access_input_iterator
{
  public:
  typedef typename Container::value_type  value_type;

  random_access_input_iterator(Container& _container, unsigned int _index = 0):
      m_container(&_container),
      m_index(_index)
  {}

  value_type& operator*()
  {
    if(m_index >= m_container->capacity())
    {
      m_container->reserve(2 * m_index + 1);
      m_container->resize(m_index+1);
    }
    else
      if(m_index >= m_container->size())
         m_container->resize(m_index+1);
    return (*m_container)[m_index];
  }

  random_access_input_iterator& operator++()
  {
    ++m_index;
    return *this;
  }

  random_access_input_iterator operator++(int)
  {
    random_access_input_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const random_access_input_iterator& raii)
  {
    CGAL_precondition(m_container == raii.m_container);
    return (m_index == raii.m_index);
  }

  bool operator!=(const random_access_input_iterator& raii)
  {
    CGAL_precondition(m_container == raii.m_container);
    return !(*this == raii);
  }

  unsigned int operator-(const random_access_input_iterator& raii)
  {
    CGAL_precondition(m_container == raii.m_container);
    return (m_index - raii.m_index);
  }


private:
  Container* m_container;
  unsigned int m_index;
};





CGAL_END_NAMESPACE

#endif
