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

#ifndef SWEEP_LINE_DO_CURVES_X_VISITOR
#define SWEEP_LINE_DO_CURVES_X_VISITOR


#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_points_visitor.h>

CGAL_BEGIN_NAMESPACE

template <class Traits>
class Sweep_line_do_curves_x_visitor
{
  typedef Sweep_line_do_curves_x_visitor<Traits>            Self;
  typedef Sweep_line_subcurve<Traits>                  Subcurve;
  typedef Sweep_line_event<Traits, Subcurve>          Event;
  typedef typename Traits::X_monotone_curve_2               X_monotone_curve_2;

   typedef Sweep_line_2_impl<Traits,
                            Event,
                            Subcurve,
                            Self,
                            CGAL_ALLOCATOR(int)>            Sweep_line;

   typedef Sweep_line_points_visitor<Traits,int>            PointsVisitor;

  public:

    Sweep_line_do_curves_x_visitor(): m_found_x(false) {}

  void attach(Sweep_line *sl)
  {
    m_sweep_line = sl;
  }

    void before_handle_event(Event* event){}
    bool after_handle_event(Event* event)
    {
      //if(event->is_internal_intersection_point())
      if(PointsVisitor::is_internal_intersection_point(event))
        m_found_x = true;
      return true;
    }

    void add_subcurve(const X_monotone_curve_2& cv,Subcurve* sc){}

    void init_subcurve(Subcurve* sc){}

    bool found_x()
    {
      return m_found_x;
    }



protected:
  bool m_found_x;
  Sweep_line* m_sweep_line;
};

CGAL_END_NAMESPACE

#endif

