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

#ifndef CGAL_PMWX_SWEEP_LINE_CURVE_H
#define CGAL_PMWX_SWEEP_LINE_CURVE_H

#include <set>

#include <CGAL/Sweep_line_2_old/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2_old/Pmwx_insert_info.h>
#include <CGAL/Sweep_line_2_old/Pmwx_sweep_line_event.h>
#include <CGAL/Sweep_line_2_old/Sweep_line_rb_tree.h>

CGAL_BEGIN_NAMESPACE

/*! @class Pmwx_sweep_line_curve 
 *  
 * a class that holds information about a curve that is added to 
 * the planar map.
 * In addition to the information that is contained in Sweep_line_subcurve,
 * when a planar map is constructed, a reference to an event that was 
 * handled last on the curve is stored. This information is used to retrieve
 * hints when a subcurve of this curve is inserted into the planar map.
 *
 * Inherits from Sweep_line_subcurve
 * \sa Sweep_line_subcurve
 */

template<class SweepLineTraits_2, class HalfedgeHandle>
class Pmwx_sweep_line_curve : public Sweep_line_subcurve<SweepLineTraits_2>
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef Sweep_line_subcurve<SweepLineTraits_2> Base;
  typedef Pmwx_sweep_line_curve<Traits, HalfedgeHandle> Self;

  typedef Status_line_curve_less_functor<Traits, Self> StatusLineCurveLess;
  typedef Red_black_tree<Self*, StatusLineCurveLess, CGAL_ALLOCATOR(int)> StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;

  typedef Pmwx_insert_info<HalfedgeHandle> PmwxInsertInfo;
  typedef Pmwx_sweep_line_event<Traits, Self> Event;


  Pmwx_sweep_line_curve(): Base(),
                           m_lastEvent(0)                         
  {}

  Pmwx_sweep_line_curve( X_monotone_curve_2 &curve): Base( curve),
                                                     m_lastEvent(0)
  {}

  
  void init(const X_monotone_curve_2 &curve)
  {
    Base::init(curve);
  }


  void set_hint(StatusLineIter hint) 
  {
    m_hint1 = hint;
  }

  StatusLineIter get_hint() const 
  {
    return m_hint1;
  }

  void set_last_event(Event *e) {
    m_lastEvent = e;
  }

  Event *get_last_event() const {
    return m_lastEvent;
  }


private:

  /*! the last event that was handled on the curve */
  Event *m_lastEvent;
  
  /*! */
  StatusLineIter m_hint1;
};


CGAL_END_NAMESPACE

#endif // CGAL_PMWX_SWEEP_LINE_CURVE_H

