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
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_POLYGON_VALIDATION_VISITOR_H
#define CGAL_GPS_POLYGON_VALIDATION_VISITOR_H

#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2_empty_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_traits_adaptor.h>

CGAL_BEGIN_NAMESPACE

/*! */
template <class ArrTraits_>
class Gps_polygon_validation_visitor : public Empty_visitor<ArrTraits_>
{
private:

  typedef ArrTraits_                                   Traits_2;

  typedef Gps_polygon_validation_visitor<Traits_2>     Self;
  typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename Traits_2::Point_2                   Point_2;

  typedef Empty_visitor<Traits_2>                      Base;
  typedef typename Base::Event                         Event;
  typedef typename Base::Subcurve                      Subcurve;
  typedef typename Base::SL_iterator                   SL_iterator;

  typedef Basic_sweep_line_2<Traits_2, Self>           Sweep_line;

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
    return (true);
  }


  bool is_valid()
  {
    return m_is_valid;
  }


protected:

  bool m_is_valid;
  bool m_is_s_simple; // is strictly simple

};

template <class ArrTraits_, class Polygon_, class PolygonWithHoles_,
          class CurveConstIterator_, class ConstructCurves_,
          class ConstructVertex_, class CheckOrientation_>
class Is_valid_2
{
private:
  
  typedef ArrTraits_                               Traits_2;
  typedef typename Traits_2::Point_2               Point_2;
  typedef typename Traits_2::X_monotone_curve_2    X_monotone_curve_2;

  typedef Polygon_                                 Polygon_2;
  typedef PolygonWithHoles_                        Polygon_with_holes_2;

  typedef CurveConstIterator_                      Curve_const_iterator;
  typedef std::pair<Curve_const_iterator,
                    Curve_const_iterator>          Cci_pair;

  typedef Gps_polygon_validation_visitor<Traits_2> Visitor;
  typedef Sweep_line_2<Traits_2, Visitor>          Sweep_line ;

public:

  typedef ConstructCurves_                         Construct_curves_2;
  typedef ConstructVertex_                         Construct_vertex_2;
  typedef CheckOrientation_                        Check_orientation_2;

private:

  // Data members:
  Construct_curves_2    construct_curves_func;
  Construct_vertex_2    construct_vertex_func;
  Check_orientation_2   check_orientation_func;

public:

  /*! Constructor. */
  Is_valid_2 (const Construct_curves_2& construct_curves,
              const Construct_vertex_2& construct_vertex,
              const Check_orientation_2& check_orientation) :
    construct_curves_func (construct_curves),
    construct_vertex_func (construct_vertex),
    check_orientation_func (check_orientation)
  {}

  /*! Check if the given polygon is valid. */
  bool operator()(const Polygon_2& pgn)
  {
    bool is_closed = _is_closed(pgn);
    CGAL_warning_msg (is_closed,
                      "The polygon's boundary is not closed.");
    if (! is_closed)
      return (false);

    bool has_valid_orientation = _has_valid_orientation(pgn);
    CGAL_warning_msg (has_valid_orientation,
                      "The polygon has a wrong orientation.");
    if (! has_valid_orientation)
      return (false);

    bool is_strictly_simple = _is_strictly_simple(pgn);
    CGAL_warning_msg (is_strictly_simple,
                      "The polygon is not strictly simple.");  
    if (! is_strictly_simple)
      return (false);   

    return (true);
  }

  /*! Check if the given polygon with holes is valid. */
  bool operator()(const Polygon_with_holes_2& pgn)
  {
    bool is_closed = _is_closed(pgn);
    CGAL_warning_msg (is_closed, 
                      "The polygon's boundary is not closed.");
    if (! is_closed)
      return (false);

    bool has_valid_orientation = _has_valid_orientation(pgn);
    CGAL_warning_msg (has_valid_orientation,
                      "The polygon has a wrong orientation.");
    if (! has_valid_orientation)
      return (false);

    bool is_simple = _is_simple(pgn);
    CGAL_warning_msg (is_simple,
                      "The polygon is not simple.");
    if (! is_simple)
      return (false);
   
    return (true);
  }

protected:

  bool _is_closed(const Polygon_2& pgn)
  {
    Cci_pair              itr_pair = construct_curves_func (pgn);
    Curve_const_iterator  begin = itr_pair.first;
    Curve_const_iterator  last = itr_pair.second;

    if (begin == last)
      return (true); // empty polygon is valid
    --last;
    
    Traits_2                    tr;
    typename Traits_2::Equal_2  equal_func = tr.equal_2_object();
    Curve_const_iterator        itr;

    for(itr = begin; itr != last; )
    {
      if (equal_func (construct_vertex_func (*itr, 0),
                      construct_vertex_func (*itr, 1)))
        return (false);

      Curve_const_iterator next = itr;
      ++next;
      if (! equal_func (construct_vertex_func (*itr, 1), 
                        construct_vertex_func (*next, 0)))
        return (false);

      itr = next;
    }

    if (equal_func (construct_vertex_func (*last, 0),
                    construct_vertex_func (*last, 1)))
      return (false);

    if (! equal_func (construct_vertex_func (*last, 1),
                      construct_vertex_func (*begin, 0)))
      return (false);

    return (true);
  }

  bool _is_closed (const Polygon_with_holes_2& pgn)
  {
    Traits_2 tr;
    if(! _is_closed (pgn.outer_boundary()))
      return (false);

    typename Polygon_with_holes_2::Hole_const_iterator    itr;

    for (itr = pgn.holes_begin(); itr != pgn.holes_end(); ++itr)
    {
      if(! _is_closed (*itr))
        return (false);
    }
    return (true);
  }

  bool _is_strictly_simple (const Polygon_2& pgn)
  {
    // Sweep the boundary curves and look for intersections.
    Cci_pair              itr_pair = construct_curves_func (pgn);
    Traits_2              traits;
    Visitor               visitor;
    Sweep_line            sweep_line (&traits, &visitor);

    visitor.sweep(itr_pair.first, itr_pair.second);
    return (visitor.is_valid());
  }

  bool _is_simple (const Polygon_with_holes_2& pgn)
  {
    // Construct a container of all boundary curves.
    Cci_pair         itr_pair = construct_curves_func (pgn.outer_boundary());
    
    std::list<X_monotone_curve_2>  curves;
    std::copy (itr_pair.first, itr_pair.second,
               std::back_inserter(curves));

    typename Polygon_with_holes_2::Hole_const_iterator  hoit;

    for (hoit = pgn.holes_begin(); hoit != pgn.holes_end(); ++hoit)
    {
      itr_pair = construct_curves_func (*hoit);
      std::copy (itr_pair.first, itr_pair.second,
                 std::back_inserter(curves));
    }

    // Perform the sweep and check fir intersections.
    Traits_2     traits;
    Visitor      visitor(false);
    Sweep_line   sweep_line (&traits, &visitor);

    visitor.sweep (curves.begin(), curves.end());
    return (visitor.is_valid());
  }

  bool _has_valid_orientation (const Polygon_2& pgn)
  {
    Cci_pair         itr_pair = construct_curves_func (pgn);
    
    if(itr_pair.first == itr_pair.second)
      return (true); // empty polygon

    return (check_orientation_func (itr_pair.first,
                                    itr_pair.second) == COUNTERCLOCKWISE);
  }

  bool _has_valid_orientation (const Polygon_with_holes_2& pgn)
  {
    // Check the orientation of the outer boundary.
    Cci_pair         itr_pair = construct_curves_func (pgn.outer_boundary());

    if ((itr_pair.first != itr_pair.second) && 
        check_orientation_func (itr_pair.first,  
                                itr_pair.second) != COUNTERCLOCKWISE)
    {
      return (false);
    }

    // Check the orientation of each of the holes.
    typename Polygon_with_holes_2::Hole_const_iterator    hoit;
    
    for (hoit = pgn.holes_begin(); hoit != pgn.holes_end(); ++hoit)
    {
      itr_pair = construct_curves_func (*hoit);

      if ((itr_pair.first !=itr_pair.second) &&
          check_orientation_func (itr_pair.first,  
                                  itr_pair.second) != CLOCKWISE)
      {
        return (false);
      }
    }

    return (true);
  }
};

CGAL_END_NAMESPACE

#endif
