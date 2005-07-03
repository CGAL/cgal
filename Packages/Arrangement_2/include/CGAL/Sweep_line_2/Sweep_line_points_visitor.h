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
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>
#include <CGAL/Sweep_line_2_empty_visitor.h>
#include <vector>
#include <iterator>

CGAL_BEGIN_NAMESPACE


template <class Traits,class OutputIerator>
class Sweep_line_points_visitor : 
  public Empty_visitor<Traits>

{
  typedef Sweep_line_points_visitor<Traits,OutputIerator>     Self;
  typedef Empty_visitor<Traits>                        Base;
  typedef typename Base::Event                         Event;
  typedef typename Base::Subcurve                      Subcurve;
  typedef typename Base::SubCurveIter                  SubCurveIter;
  typedef typename Base::SL_iterator                   SL_iterator;

   
  typedef typename Traits::X_monotone_curve_2          X_monotone_curve_2;
  typedef typename Traits::Point_2                     Point_2;


  typedef Basic_sweep_line_2<Traits,
                             Self>                     Sweep_line;

             



public:


  Sweep_line_points_visitor(OutputIerator out,
                            bool endpoints,
                            Traits* tr): m_traits(tr),
                                         m_out(out),
                                         m_includeEndPoints(endpoints)
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
    if(m_includeEndPoints ||
       event->is_intersection() ||
       event->is_weak_intersection())
      *m_out++ = event->get_point();
    return true;
  }

 void add_subcurve(X_monotone_curve_2 cv,Subcurve* sc){}

  OutputIerator get_output_iterator()
  {
    return m_out;
  }

 



protected:

    Traits*          m_traits;
    OutputIerator    m_out;
    bool             m_includeEndPoints;

  };

CGAL_END_NAMESPACE

#endif
