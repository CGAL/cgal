// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Id$ $Date$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef GPS_POLYGON_VALIDATION_VISITOR
#define GPS_POLYGON_VALIDATION_VISITOR

#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2_empty_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_traits_adaptor.h>

CGAL_BEGIN_NAMESPACE

/*! */
template <class Traits>
class Gps_polygon_validation_visitor : public Empty_visitor<Traits>
{
  typedef Gps_polygon_validation_visitor<Traits>       Self;
  typedef typename Traits::X_monotone_curve_2          X_monotone_curve_2;
  typedef typename Traits::Point_2                     Point_2;

  typedef Empty_visitor<Traits>                        Base;
  typedef typename Base::Event                         Event;
  typedef typename Base::Subcurve                      Subcurve;
  typedef typename Base::SL_iterator                   SL_iterator;

  typedef Basic_sweep_line_2<Traits, Self>             Sweep_line;

  public:

    Gps_polygon_validation_visitor(bool is_s_simple = true): 
        m_is_valid(true),
        m_is_s_simple(is_s_simple)
    {}



  template <class XCurveIterator>
  void sweep(XCurveIterator begin, XCurveIterator end)
  {
    //Perform the sweep
    reinterpret_cast<Sweep_line*>(this->m_sweep_line)->sweep(begin, end);
  }

  bool after_handle_event(Event* event,SL_iterator iter, bool flag)
  {
    if(event->is_intersection() ||
       event->is_weak_intersection() ||
       event->is_overlap())
    {
      m_is_valid = false;
      reinterpret_cast<Sweep_line*>(this->m_sweep_line) -> stop_sweep();
    }
    else
      if(m_is_s_simple && 
         (event->get_num_right_curves() + event->get_num_left_curves()) != 2)
      {
         m_is_valid = false;
         reinterpret_cast<Sweep_line*>(this->m_sweep_line) -> stop_sweep();
      }
    return true;
  }


  bool is_valid()
  {
    return m_is_valid;
  }


protected:

  bool m_is_valid;
  bool m_is_s_simple; // is strictly simple

};

template <class Traits>
class Is_valid_2
{
  typedef typename Traits::Polygon_2               Polygon_2;
  typedef typename Traits::Polygon_with_holes_2    Polygon_with_holes_2;
  typedef typename Traits::X_monotone_curve_2      X_monotone_curve_2;
  typedef typename Traits::Point_2                 Point_2;
  typedef typename Traits::Curve_const_iterator    Curve_const_iterator;
  typedef typename Traits::Construct_curves_2      Construct_curves_2;
  typedef Gps_polygon_validation_visitor<Traits>   Visitor;
  typedef Sweep_line_2<Traits, Visitor>            Sweep_line ;

public:

  bool operator()(const Polygon_2& pgn)
  {
    bool _is_closed = is_closed(pgn);
    CGAL_warning_msg(_is_closed, "The polygon's boundary is not closed.");
    if(!_is_closed)
      return false;

    bool _had_valid_orientation = had_valid_orientation(pgn);
    CGAL_warning_msg(_had_valid_orientation,
                     "The polygon has a wrong orientation.");
    if(!_had_valid_orientation)
      return false;

    bool _is_strictly_simple = is_strictly_simple(pgn);
    CGAL_warning_msg(_is_strictly_simple,
                     "The polygon is not strictly simple.");  
    if(!_is_strictly_simple)
      return false;   

    return true;
  }

  bool operator()(const Polygon_with_holes_2& pgn)
  {
    bool _is_closed = is_closed(pgn);
    CGAL_warning_msg(_is_closed, "The polygon's boundary is not closed.");
    if(!_is_closed)
      return false;

    bool _had_valid_orientation = had_valid_orientation(pgn);
    CGAL_warning_msg(_had_valid_orientation,
                     "The polygon has a wrong orientation.");
    if(!_had_valid_orientation)
      return false;

    bool _is_simple = is_simple(pgn);
    CGAL_warning_msg(_is_simple, "The polygon is not simple.");
    if(!_is_simple)
      return false;
   
    return true;
  }

protected:

  bool is_closed(const Polygon_2& pgn)
  {
    Gps_traits_adaptor<Traits>  tr;
    typename Gps_traits_adaptor<Traits>::Construct_vertex_2 ctr_v = 
      tr.construct_vertex_2_object();
    std::pair<Curve_const_iterator,
              Curve_const_iterator> itr_pair =
              tr.construct_curves_2_object()(pgn);
    Curve_const_iterator begin = itr_pair.first;
    Curve_const_iterator last = itr_pair.second;
    if(begin == last)
      return true; // empty polygon is valid
    --last;
    
    for(Curve_const_iterator itr = begin; itr != last; )
    {
      if(tr.equal_2_object()(ctr_v(*itr, 0), ctr_v(*itr ,1)))
      {
        return false;
      }
      Curve_const_iterator next = itr;
      ++next;
      if(!tr.equal_2_object()(ctr_v(*itr ,1), ctr_v(*next, 0)))
        return false;
      itr = next;
    }
    if(tr.equal_2_object()(ctr_v(*last, 0), ctr_v(*last, 1)))
    {
      return false;
    }
    if(!tr.equal_2_object()(ctr_v(*last, 1), ctr_v(*begin, 0)))
      return false;
    return true;
  }

  bool is_closed(const Polygon_with_holes_2& pgn)
  {
    Traits tr;
    typedef typename Polygon_with_holes_2::Hole_const_iterator    HCI;
    if(!is_closed(pgn.outer_boundary()))
      return false;

    for(HCI itr = pgn.holes_begin(); itr != pgn.holes_end(); ++itr)
    {
      if(!is_closed(*itr))
        return false;
    }
    return true;
  }

  bool is_strictly_simple(const Polygon_2& pgn)
  {
    Traits tr;
    std::pair<Curve_const_iterator,
            Curve_const_iterator> itr_pair =
            tr.construct_curves_2_object()(pgn);
    // Define the sweep-line types:
    // Perform the sweep and obtain the subcurves.
    Visitor     visitor;
    Sweep_line  sweep_line (&tr, &visitor);
    visitor.sweep(itr_pair.first, itr_pair.second);
    return visitor.is_valid();
  }

  bool is_simple(const Polygon_with_holes_2& pgn)
  {
    Traits tr;
    std::list<X_monotone_curve_2>  curves;
    std::pair<Curve_const_iterator,
              Curve_const_iterator> itr_pair = 
      tr.construct_curves_2_object()(pgn.outer_boundary());
    std::copy(itr_pair.first, itr_pair.second, std::back_inserter(curves));

    typename Polygon_with_holes_2::Hole_const_iterator i = pgn.holes_begin();
    for(; i != pgn.holes_end(); ++i)
    {
      itr_pair = tr.construct_curves_2_object()(*i);
      std::copy(itr_pair.first, itr_pair.second, std::back_inserter(curves));
    }
    // Define the sweep-line types:
    // Perform the sweep and obtain the subcurves.
    Visitor     visitor(false);
    Sweep_line  sweep_line (&tr, &visitor);
    visitor.sweep(curves.begin(), curves.end());
    return visitor.is_valid();
  }

  bool had_valid_orientation(const Polygon_2& pgn)
  {
    Gps_traits_adaptor<Traits>  tr;
    Construct_curves_2 ctr_curves = tr.construct_curves_2_object();
    std::pair<Curve_const_iterator, Curve_const_iterator> itr_pair = 
      ctr_curves(pgn);
    
    if(itr_pair.first == itr_pair.second)
      return true; // empty polygon

    return (tr.orientation_2_object()(itr_pair.first,
                                      itr_pair.second) == COUNTERCLOCKWISE);
  }

  bool had_valid_orientation(const Polygon_with_holes_2& pgn)
  {
    Gps_traits_adaptor<Traits>  tr;

    Construct_curves_2 ctr_curves = tr.construct_curves_2_object();
    std::pair<Curve_const_iterator, Curve_const_iterator> itr_pair = 
      ctr_curves(pgn.outer_boundary());
    if((itr_pair.first != itr_pair.second) && 
       tr.orientation_2_object()(itr_pair.first,  
                                 itr_pair.second) != COUNTERCLOCKWISE)
      return false;

    typedef typename Polygon_with_holes_2::Hole_const_iterator    HCI;
    for(HCI hit = pgn.holes_begin(); hit != pgn.holes_end(); ++hit)
    {
      itr_pair = ctr_curves(*hit);
      if((itr_pair.first !=itr_pair.second) &&
         tr.orientation_2_object()(itr_pair.first,  
                                   itr_pair.second) != CLOCKWISE)
        return false;
    }

    return true;
  }
};

CGAL_END_NAMESPACE

#endif
