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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef ARR_BATCHED_POINT_LOCATION_EVENT_H
#define ARR_BATCHED_POINT_LOCATION_EVENT_H


#include <CGAL/Sweep_line_2/Sweep_line_event.h>

CGAL_BEGIN_NAMESPACE

template<class SweepLineTraits_2, class CurveWrap>
class Arr_batched_point_location_event :
  public Sweep_line_event<SweepLineTraits_2, CurveWrap>
{
public:

  typedef SweepLineTraits_2 Traits;
  typedef CurveWrap SubCurve;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Point_2 Point_2;

  typedef Sweep_line_event<Traits, SubCurve> Base;
  typedef Arr_batched_point_location_event<Traits, SubCurve> Self;


  /*! constructor*/
  Arr_batched_point_location_event(): m_is_query(false)
  {}

  /*! destructor */
  ~Arr_batched_point_location_event(){}

  void init(const Point_2 &point)
  {
    Base::init(point);
  }


  bool is_query() const
  {
    return m_is_query;
  }
  
  void set_query() 
  {
    m_is_query = true;
  }


  bool has_curves() const
  {
    return (Base::has_left_curves() || Base::has_right_curves());
  }



protected:

  bool m_is_query;
};


CGAL_END_NAMESPACE

#endif

