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
// $URL$
// $Id$
// 
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>,
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_SWEEP_LINE_FUNCTORS_H
#define CGAL_SWEEP_LINE_FUNCTORS_H

#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
//#include <CGAL/Sweep_line_2/Sweep_line_event.h>

CGAL_BEGIN_NAMESPACE

template <class SweepLineTraits_2, class Event>
class Event_less_functor 
{
public:
  typedef SweepLineTraits_2                      Traits;
  typedef typename Traits::Point_2               Point_2;
  typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Traits::Has_boundary_category Has_boundary_category;
  
  Event_less_functor(Traits * t) : m_traits(t)
  {}
  
  //this operator is called by the multiset assertions only in
  //debug mode (to verify that event was inserted at the right place).
  Comparison_result operator()(const Event* e1, const Event* e2) const
  {
    
    if(e1->is_finite() && e2->is_finite())
      return (m_traits->compare_xy_2_object()(e1->point(),
                                              e2->point()));
    if(e1->is_finite())
      return ( this->operator()(e1->point(), e2) );
    
    if(e2->is_finite())
      return (CGAL::opposite(this->operator()(e2->point(), e1)));
    
    
    
    return (compare_unbounded_curve_and_event(e1->get_unbounded_curve(),
                                              e2,
                                              e1->infinity_at_x(),
                                              e1->infinity_at_y(),
                                              curve_end(e1)));
  }

  Comparison_result operator()(const Point_2& pt, const Event* e2) const
  {
    // e1 is a normal event (given as 'pt').

    if(e2->is_finite())
      // e2 is also a normal event, use compare_xy.
      return (m_traits->compare_xy_2_object()(pt, e2->point()));

    Boundary_type e2_x = e2->infinity_at_x();
    switch (e2_x)
    {
    case MINUS_INFINITY :
      // e2 x-coord is at minus infinty, so e1 is larger then e2.
      return LARGER;

    case PLUS_INFINITY :
      // e2 x-coord is at plus infinty, so e1 is smaller then e2.
      return SMALLER;

    default:
      break;
    }

    Curve_end index = curve_end(e2);
   
    Comparison_result res = 
      m_traits->compare_x_2_object()(pt, e2->get_unbounded_curve(), index);
    if(res != EQUAL)
      return res;

    Boundary_type e2_y = e2->infinity_at_y();
    switch (e2_y)
    {
    case MINUS_INFINITY :
      return LARGER;
    
    case PLUS_INFINITY :
      return SMALLER;
    
    default:
      // doesnt suppose to reach here at all
      CGAL_assertion(e2_y != NO_BOUNDARY);
      return SMALLER;
    }
  }

  Comparison_result operator()(const X_monotone_curve_2& cv, const Event* e2) const
  {
    return compare_unbounded_curve_and_event(cv,
                                             e2,
                                             m_boundary_in_x,
                                             m_boundary_in_y,
                                             m_index);
  }


  Comparison_result 
    compare_unbounded_curve_and_event(const X_monotone_curve_2& cv,
                                      const Event* e2,
                                      Boundary_type boundary_in_x,
                                      Boundary_type boundary_in_y,
                                      Curve_end index) const
  {
    switch(boundary_in_x)
    {
    case MINUS_INFINITY:
      if(e2->infinity_at_x() == MINUS_INFINITY)
        return (m_traits->compare_y_at_x_2_object()(cv,
                                                    e2->get_unbounded_curve(),
                                                    MIN_END));
      return SMALLER;
      
    case PLUS_INFINITY:
      if(e2->infinity_at_x() == PLUS_INFINITY)
        return (m_traits->compare_y_at_x_2_object()(cv,
                                                    e2->get_unbounded_curve(),
                                                    MAX_END));
      return LARGER;

    case NO_BOUNDARY:
    default:
      break;
    }

    //if we have reached here, then e1 is event at y=+-oo
    CGAL_assertion(boundary_in_y == MINUS_INFINITY ||
                   boundary_in_y == PLUS_INFINITY);
    switch(e2->infinity_at_x())
    {
    case MINUS_INFINITY:
      return LARGER;
    case PLUS_INFINITY:
      return SMALLER;
    case NO_BOUNDARY:
    default:
      break;
    }

    // e2 is a normal event or event at y=+-00
    Comparison_result res;
    switch(e2->infinity_at_y())
    {
    case MINUS_INFINITY:
      res =  m_traits->compare_x_2_object()(cv,
                                            index,
                                            e2->get_unbounded_curve(),
                                            curve_end(e2));
      if(res != EQUAL)
        return res;
      
      if(boundary_in_y == MINUS_INFINITY)
        return EQUAL;

      return LARGER;

    case PLUS_INFINITY:
      res =  m_traits->compare_x_2_object()(cv,
                                            index,
                                            e2->get_unbounded_curve(),
                                            curve_end(e2));
       if(res != EQUAL)
         return res;
       
       if(boundary_in_y == PLUS_INFINITY)
        return EQUAL;

      return SMALLER;

    case NO_BOUNDARY:
      // e2 is a normal event
      res = m_traits->compare_x_2_object()(e2->point(), cv, index);
      if(res != EQUAL)
        return CGAL::opposite(res);

      if(boundary_in_y == MINUS_INFINITY)
        return SMALLER;

      CGAL_assertion(boundary_in_y == PLUS_INFINITY);
      return LARGER;
     default:
        CGAL_assertion(false); 
	break;

    }

    //doesnt suppose to reach here (all options have been handles by now)
    CGAL_assertion(false);
    return SMALLER;
  }

  void set_boundary_in_x(Boundary_type s)
  {
    m_boundary_in_x = s;
  }

  void set_boundary_in_y(Boundary_type s)
  {
    m_boundary_in_y = s;
  }

  void set_index(Curve_end index)
  {
    m_index = index;
  }

  Curve_end curve_end(const Event* e) const
  {
    CGAL_assertion(!e->is_finite());
    return ((e->is_left_end()) ? MIN_END : MAX_END);
  }


private:

  /*! a pointer to a traits object */
  Traits*        m_traits;
  
  // when comparing using operator()(cv, e) we store here information 
  Boundary_type     m_boundary_in_x;
  Boundary_type     m_boundary_in_y;
  Curve_end         m_index;
};


// forward decleration for Sweep_line_event
template <class Traits_, class Subcurve_>
class Sweep_line_event;


template <class SweepLineTraits_2, class Subcurve> 
class Status_line_curve_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef Arr_traits_basic_adaptor_2<Traits>            Traits_adaptor;

  typedef typename Traits_adaptor::Point_2 Point_2;
  typedef typename Traits_adaptor::X_monotone_curve_2 X_monotone_curve_2;
  typedef Sweep_line_event<Traits, Subcurve>  Event;
  

  template <class Sweep_event>
  Status_line_curve_less_functor(Traits_adaptor *t, Sweep_event** e_ptr) : m_traits(t),
                                                            m_curr_event(reinterpret_cast<Event**>(e_ptr))
  {}

  // this comprator should never be called in release mode (only debug mode)
  Comparison_result operator()(const Subcurve * c1, const Subcurve * c2) const
  {
    // in case to two curves are right curves at the same event, compare
    // to the right of the event point.
    if(std::find((*m_curr_event)->right_curves_begin(), 
                 (*m_curr_event)->right_curves_end(),
                 c1) != (*m_curr_event)->right_curves_end() &&
       std::find((*m_curr_event)->right_curves_begin(), 
                 (*m_curr_event)->right_curves_end(),
                 c2) != (*m_curr_event)->right_curves_end())

     return m_traits->compare_y_at_x_right_2_object()
       (c1->get_last_curve(), c2->get_last_curve(), (*m_curr_event)->point());

    Boundary_type x_inf_c1 =
      m_traits->boundary_in_x_2_object()(c1->get_last_curve(), MIN_END);
    Boundary_type y_inf_c1 =
      m_traits->boundary_in_y_2_object()(c1->get_last_curve(), MIN_END);

    if(x_inf_c1 == NO_BOUNDARY && y_inf_c1 == NO_BOUNDARY)
      return m_traits->compare_y_at_x_2_object()
        (m_traits->construct_min_vertex_2_object()(c1->get_last_curve()),
                                                   c2->get_last_curve());
    
    CGAL_assertion(x_inf_c1 != PLUS_INFINITY);
    if(x_inf_c1 == MINUS_INFINITY)
      return LARGER;

    if(y_inf_c1 == MINUS_INFINITY) 
       return SMALLER;

    CGAL_assertion(y_inf_c1 == PLUS_INFINITY);
    return LARGER;
  }

  Comparison_result operator()(const Point_2& pt, const Subcurve * c2) const
  {
    return (m_traits->compare_y_at_x_2_object()(pt,c2->get_last_curve()));
  }

  /*! a pointer to a traits object */
  Traits_adaptor *  m_traits;

  Event**   m_curr_event;  
};


template <class _Subcurve>
class Curves_pair
{
  public:
  Curves_pair(){}

  Curves_pair (_Subcurve* sc1, _Subcurve* sc2)
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

  bool operator()(const CurvesPair& pair1, const CurvesPair& pair2) const
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

template <class _Subcurve>
class Curves_pair_hash_functor
{
  typedef Curves_pair<_Subcurve> CurvesPair;

public:

  size_t operator() (const CurvesPair& pair) const
  {
    const size_t  half_n_bits = sizeof(size_t) * 8 / 2;
    const size_t  val1 = reinterpret_cast<size_t> (pair.first());
    const size_t  val2 = reinterpret_cast<size_t> (pair.second());

    return (((val1 << half_n_bits) | (val1 >> half_n_bits)) ^ val2);
  }
};

template <class _Subcurve>
class Curves_pair_equal_functor
{
  typedef Curves_pair<_Subcurve> CurvesPair;

public:

  bool operator() (const CurvesPair& pair1, const CurvesPair& pair2) const
  {
    return (pair1.first() == pair2.first() &&
            pair1.second() == pair2.second());
  }
};

template <class Container>
class random_access_input_iterator
{
  public:
  typedef typename Container::value_type  value_type;

  random_access_input_iterator()
  {}

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

  random_access_input_iterator& operator--()
  {
    --m_index;
    return *this;
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
