
#ifndef GPS_POLYGON_VALIDATION_VISITOR
#define GPS_POLYGON_VALIDATION_VISITOR

#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2_empty_visitor.h>

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
  typedef Gps_polygon_validation_visitor<Traits>   Visitor;
  typedef Sweep_line_2<Traits, Visitor>            Sweep_line ;

public:

  bool operator()(const Polygon_2& pgn)
  {
    if(!is_closed(pgn))
      return false;
    if(!is_strictly_simple(pgn))
      return false;
    return true;
  }

  bool operator()(const Polygon_with_holes_2& pgn)
  {
    if(!is_closed(pgn))
      return false;
    if(!is_simple(pgn))
      return false;
    return true;
  }

protected:

  Point_2 source(const X_monotone_curve_2& cv)
  {
    Traits tr;
    if(tr.compare_endpoints_xy_2_object()(cv) == SMALLER)
      return (tr.construct_min_vertex_2_object()(cv));
    return (tr.construct_max_vertex_2_object()(cv));
  }

  Point_2 target(const X_monotone_curve_2& cv)
  {
    Traits tr;
    if(tr.compare_endpoints_xy_2_object()(cv) == SMALLER)
      return (tr.construct_max_vertex_2_object()(cv));
    return (tr.construct_min_vertex_2_object()(cv));
  }

  bool is_closed(const Polygon_2& pgn)
  {
    Traits tr;
    std::pair<Curve_const_iterator,
            Curve_const_iterator> itr_pair =
            tr.construct_curves_2_object()(pgn);
    Curve_const_iterator begin = itr_pair.first;
    Curve_const_iterator last = itr_pair.second;
    --last;
    
    for(Curve_const_iterator itr = begin; itr != last; )
    {
      if(tr.equal_2_object()(source(*itr), target(*itr)))
      {
        return false;
      }
      Curve_const_iterator next = itr;
      ++next;
      if(!tr.equal_2_object()(target(*itr), source(*next)))
        return false;
      itr = next;
    }
    if(tr.equal_2_object()(source(*last), target(*last)))
    {
      return false;
    }
    if(!tr.equal_2_object()(target(*last), source(*begin)))
      return false;
    return true;
  }

  bool is_closed(const Polygon_with_holes_2& pgn)
  {
    Traits tr;
    typedef typename Polygon_with_holes_2::Holes_const_iterator    HCI;
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
    std::pair<Curve_const_iterator,
            Curve_const_iterator> itr_pair =
            tr.construct_curves_2_object()(pgn);
    // Define the sweep-line types:
    // Perform the sweep and obtain the subcurves.
    Visitor     visitor(false);
    Sweep_line  sweep_line (&tr, &visitor);
    visitor.sweep(itr_pair.first, itr_pair.second);
    return visitor.is_valid();
  }
};

CGAL_END_NAMESPACE

#endif
