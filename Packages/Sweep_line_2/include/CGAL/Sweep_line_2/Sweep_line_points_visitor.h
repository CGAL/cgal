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

#ifndef SWEEP_LINE_POINTS_NOTIFICATION_H
#define SWEEP_LINE_POINTS_NOTIFICATION_H

#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>

CGAL_BEGIN_NAMESPACE

template <class Traits,class OutputIerator>
class Sweep_line_points_visitor 
{
  typedef Sweep_line_points_visitor<Traits,OutputIerator>     Self;
  typedef Sweep_line_subcurve<Traits>                    Subcurve;
  typedef Sweep_line_event<Traits, Subcurve>            Event;

  typedef Sweep_line_2_impl<Traits,
                            Event,
                            Subcurve,
                            Self,
                            CGAL_ALLOCATOR(int)>              Sweep_line;

   
  typedef typename Traits::X_monotone_curve_2                 X_monotone_curve_2;



public:


  Sweep_line_points_visitor(OutputIerator out,
                            bool endpoints): m_out(out),
                                             m_includeEndPoints(endpoints)
  {}

  void attach(Sweep_line *sl)
  {
    m_sweep_line = sl;
  }

       
  void before_handle_event(Event* event){}

  bool after_handle_event(Event* event)
  {
    if(m_includeEndPoints || event->is_internal_intersection_point())
      *m_out++ = event->get_point();
    return true;
  }

 void add_subcurve(X_monotone_curve_2 cv,Subcurve* sc){}

  OutputIerator get_output_iterator()
  {
    return m_out;
  }

  void init_subcurve(Subcurve* sc){}

  


  protected:

    OutputIerator m_out;
    bool m_includeEndPoints;
    Sweep_line* m_sweep_line;

  };

CGAL_END_NAMESPACE

#endif
