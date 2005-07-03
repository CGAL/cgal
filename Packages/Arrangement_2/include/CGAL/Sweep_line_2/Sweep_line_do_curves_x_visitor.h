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
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>
#include <CGAL/Sweep_line_2_empty_visitor.h>
#include <vector>
#include <iterator>

CGAL_BEGIN_NAMESPACE

template <class Traits>
class Sweep_line_do_curves_x_visitor : public Empty_visitor<Traits>
{
  typedef Sweep_line_do_curves_x_visitor<Traits>       Self;
  typedef typename Traits::X_monotone_curve_2          X_monotone_curve_2;
  typedef typename Traits::Point_2                     Point_2;

  typedef Empty_visitor<Traits>                        Base;
  typedef typename Base::Event                         Event;
  typedef typename Base::Subcurve                      Subcurve;
  typedef typename Base::SL_iterator                   SL_iterator;

  typedef Basic_sweep_line_2<Traits, Self>             Sweep_line;

  public:

    Sweep_line_do_curves_x_visitor(Traits *tr): m_traits(tr),
                                                m_found_x(false)
    {}



  template <class CurveIterator>
  void sweep(CurveIterator begin, CurveIterator end)
  {
    std::vector<X_monotone_curve_2> curves_vec;
    std::vector<Point_2> points_vec;
    curves_vec.reserve(std::distance(begin,end));
    make_x_monotone(begin,
                    end,
                    std::back_inserter(curves_vec),
                    std::back_inserter(points_vec),
                    m_traits);
   
    //Perform the sweep
    reinterpret_cast<Sweep_line*>(m_sweep_line) -> sweep(curves_vec.begin(),
                                                         curves_vec.end(),
                                                         points_vec.begin(),
                                                         points_vec.end());
  }

  bool after_handle_event(Event* event,SL_iterator iter, bool flag)
  {
    if(event->is_intersection() || event->is_weak_intersection())
    {
      m_found_x = true;
      reinterpret_cast<Sweep_line*>(m_sweep_line) -> stop_sweep();
    }
    return true;
  }


  bool found_x()
  {
    return m_found_x;
  }



protected:

  Traits* m_traits;

  bool m_found_x;

};

CGAL_END_NAMESPACE

#endif

