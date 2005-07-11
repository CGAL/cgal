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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_CONSTRUCTION_CURVE_H
#define CGAL_ARR_CONSTRUCTION_CURVE_H

#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>

CGAL_BEGIN_NAMESPACE

/*! @class Arr_construction_curve 
 *  
 * a class that holds information about a curve that is added to 
 * the arrangement.
 * In addition to the information that is contained in Sweep_line_subcurve,
 * when an arrangement is constructed, a pointer to the last handled event  
 * on the curve is stored. This information is used to retrieve
 * hints when a subcurve of this curve is inserted into the planar map.
 *
 * Inherits from Sweep_line_subcurve
 * \sa Sweep_line_subcurve
 */

template<class SweepLineTraits_2>
class Arr_construction_curve : public Sweep_line_subcurve<SweepLineTraits_2>
{
public:
  typedef SweepLineTraits_2                        Traits;
  typedef typename Traits::Point_2                 Point_2;
  typedef typename Traits::X_monotone_curve_2      X_monotone_curve_2;

  typedef Sweep_line_subcurve<SweepLineTraits_2>   Base;
  typedef Arr_construction_curve<Traits>           Self;

  typedef typename Base::StatusLineIter            StatusLineIter;
  typedef void*                                    Event_ptr;


  Arr_construction_curve(): Base(),
                            m_lastEvent(0)                         
  {}

  Arr_construction_curve(X_monotone_curve_2 &curve): Base( curve),
                                                     m_lastEvent(0)
  {}

  
  template<class SweepEvent>
  void init(const X_monotone_curve_2 &curve,
            SweepEvent* left,
            SweepEvent* right)
  {
    Base::init(curve, left, right);
    m_lastEvent = left;
    
  }

  void set_last_event(Event_ptr e) {
    m_lastEvent = e;
  }

  Event_ptr get_last_event() const {
    return m_lastEvent;
  }


protected:

  /*! the last event that was handled on the curve */
  Event_ptr  m_lastEvent;
  
};


CGAL_END_NAMESPACE

#endif 
