// Copyright (c) 2008  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s): Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein        <wein@post.tau.ac.il>
//                 Boris Kozorovitzky <boriskoz@post.tau.ac.il>
//                 Guy Zucker <guyzucke@post.tau.ac.il>

#ifndef CGAL_GPS_POLYGON_VALIDATION_2_H
#define CGAL_GPS_POLYGON_VALIDATION_2_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Boolean_set_operations_2/Gps_traits_adaptor.h>
#include <CGAL/Boolean_set_operations_2/Gps_default_dcel.h>
#include <CGAL/Boolean_set_operations_2/Gps_on_surface_base_2.h>

#include <CGAL/Arrangement_2/Arr_default_planar_topology.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/Arr_naive_point_location.h>


#include <iostream>
#include <list>
#include <iterator>


namespace CGAL {

/*Arrangement is templated with extended face dcel*/
template <typename Arrangement_2>
class ValidationOverlayTraits :
  public CGAL::Arr_default_overlay_traits<Arrangement_2>
{
public:
  typedef CGAL::Arr_default_overlay_traits<Arrangement_2>       Base;
  typedef typename Base::Face_handle_A                          Face_handle_A;
  typedef typename Base::Face_handle_B                          Face_handle_B;
  typedef typename Base::Face_handle_R                          Face_handle_R;

  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Halfedge_const_handle
    Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle
    Face_const_handle;
  typedef typename Arrangement_2::Inner_ccb_const_iterator
    Inner_ccb_const_iterator;

  /* red faces source is the arrangement of holes. The blue faces (face) are
   * caused by the PWH's outer boundary
   */
  virtual void create_face(Face_handle_A red_face, Face_handle_B blue_face,
                           Face_handle_R /*r_face*/) const
  {
    if ((red_face->contained()==true) && (blue_face->contained()==false)) {
      hole_overlap = true;
    }
  }

public:
  ValidationOverlayTraits() : hole_overlap(false) {}
  bool getHoleOverlap() {
    return hole_overlap;
  }
  void setHoleOverlap(bool b) {
    hole_overlap = b;
    return;
  }
private:
  mutable bool hole_overlap;
};

/*! \class
 * A visitor used for checking whether the edges of a polygon are
 * non-intersecting.
 */
template <class ArrTraits_>
class Gps_polygon_validation_visitor :
  public Sweep_line_empty_visitor<ArrTraits_>
{
private:
  typedef ArrTraits_                                   Traits_2;
  typedef Gps_polygon_validation_visitor<Traits_2>     Self;
  typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename Traits_2::Point_2                   Point_2;

  typedef Sweep_line_empty_visitor<Traits_2>           Base;
  typedef typename Base::Event                         Event;
  typedef typename Base::Subcurve                      Subcurve;
  typedef typename Base::Status_line_iterator          SL_iterator;

  typedef Basic_sweep_line_2<Traits_2, Self>           Sweep_line;

public:
  enum Error_code {
    ERROR_NONE = 0,
    ERROR_EDGE_INTERSECTION,
    ERROR_EDGE_VERTEX_INTERSECTION,
    ERROR_EDGE_OVERLAP,
    ERROR_VERTEX_INTERSECTION
  };

  Gps_polygon_validation_visitor(bool is_s_simple = true) :
    m_is_valid(true),
    m_is_s_simple(is_s_simple),
    m_error_code(ERROR_NONE)
  {}

  template <class XCurveIterator>
  void sweep(XCurveIterator begin, XCurveIterator end)
  {
    //Perform the sweep
    reinterpret_cast<Sweep_line*>(this->m_sweep_line)->sweep(begin, end);
  }

  bool after_handle_event(Event* event, SL_iterator, bool)
  {
    if (event->is_intersection()) {
      m_error_code = ERROR_EDGE_INTERSECTION;
      m_is_valid = false;
      reinterpret_cast<Sweep_line*>(this->m_sweep_line)->stop_sweep();
    }
    else if (event->is_weak_intersection()) {
      m_error_code = ERROR_EDGE_VERTEX_INTERSECTION;
      m_is_valid = false;
      reinterpret_cast<Sweep_line*>(this->m_sweep_line)->stop_sweep();
    }
    else if (event->is_overlap()) {
      m_error_code = ERROR_EDGE_OVERLAP;
      m_is_valid = false;
      reinterpret_cast<Sweep_line*>(this->m_sweep_line)->stop_sweep();
    } else {
      if (m_is_s_simple &&
          (event->number_of_right_curves() + event->number_of_left_curves()) !=
          2)
      {
        m_error_code = ERROR_VERTEX_INTERSECTION;
        m_is_valid = false;
        reinterpret_cast<Sweep_line*>(this->m_sweep_line)->stop_sweep();
      }
    }
    return true;
  }

  bool is_valid() const { return m_is_valid; }
  Error_code error_code() const { return m_error_code; }

protected:
  bool m_is_valid;
  bool m_is_s_simple; // is strictly simple

private:
  Error_code m_error_code;
};


//Traits_2 templates the General_polygon_set_2 Traits.
//These include types for polygon and PWH.
template <typename Traits_2>
bool is_closed_polygon(const typename Traits_2::Polygon_2& pgn,
                       const Traits_2& traits)
{
  typedef Gps_traits_adaptor<Traits_2>                  Traits_adapter_2;
  typedef typename Traits_2::Curve_const_iterator       Curve_const_iterator;
  typedef std::pair<Curve_const_iterator,Curve_const_iterator>    Cci_pair;
  typedef typename Traits_adapter_2::Construct_vertex_2 Construct_vertex_2;

  Cci_pair itr_pair = traits.construct_curves_2_object()(pgn);
  Curve_const_iterator begin = itr_pair.first;
  Curve_const_iterator end = itr_pair.second;

  if (begin == end) return true;  // An empty polygon is valid.

  Traits_adapter_2 traits_adapter(traits);
  typename Traits_2::Equal_2 equal_func = traits.equal_2_object();
  Curve_const_iterator curr, next;
  Construct_vertex_2 construct_vertex_func =
    traits_adapter.construct_vertex_2_object();
  curr = next = begin;
  ++next;

  if (next == end)
    return false; // A polygon cannot have just a single edge.

  while (next != end) {
    // Make sure that the current target equals the next source.
    if (equal_func(construct_vertex_func(*curr, 0),
                   construct_vertex_func(*curr, 1)))
      return false;

    if (! equal_func(construct_vertex_func(*curr, 1),
                     construct_vertex_func(*next, 0)))
      return false;

    // Move to the next pair of edges.
    curr = next;
    ++next;
  }

  // Make sure that the last target equals the first source.
  if (equal_func (construct_vertex_func (*curr, 0),
                  construct_vertex_func (*curr, 1)))
    return false;

  if (! equal_func (construct_vertex_func (*curr, 1),
                    construct_vertex_func (*begin, 0)))
    return false;

  return true;
}

// Previously known as is_strictly_simple
template <typename Traits_2>
bool is_simple_polygon(const typename Traits_2::Polygon_2& pgn,
                       const Traits_2& traits)
{
  typedef typename Traits_2::Curve_const_iterator    Curve_const_iterator;
  typedef std::pair<Curve_const_iterator,Curve_const_iterator>    Cci_pair;

  // Sweep the boundary curves and look for intersections.
  typedef Gps_polygon_validation_visitor<Traits_2>  Visitor;
  typedef Sweep_line_2<Traits_2, Visitor>           Sweep_line;

  Cci_pair              itr_pair = traits.construct_curves_2_object()(pgn);
  Visitor               visitor;
  Sweep_line            sweep_line (&traits, &visitor);

  visitor.sweep(itr_pair.first, itr_pair.second);
  if (!visitor.is_valid()) {
    switch (visitor.error_code()) {
     case Visitor::ERROR_NONE: break;
     case Visitor::ERROR_EDGE_INTERSECTION:
      CGAL_warning_msg(false, "The polygon boundary self intersects at edges.");
      break;

     case Visitor::ERROR_EDGE_VERTEX_INTERSECTION:
      CGAL_warning_msg(false, "The polygon boundary self (weakly) intersects.");
      break;

     case Visitor::ERROR_EDGE_OVERLAP:
      CGAL_warning_msg(false, "The polygon boundary self overlaps.");
      break;

     case Visitor::ERROR_VERTEX_INTERSECTION:
      CGAL_warning_msg(false, "The polygon boundary intersects at vertices.");
      break;
    }
    return false;
  }
  return true;
}

template <typename Traits_2>
bool has_valid_orientation_polygon(const typename Traits_2::Polygon_2& pgn,
                                   const Traits_2& traits)
{

  typedef Gps_traits_adaptor<Traits_2>              Traits_adapter_2;
  typedef typename Traits_2::Curve_const_iterator    Curve_const_iterator;
  typedef std::pair<Curve_const_iterator,Curve_const_iterator>    Cci_pair;

  Cci_pair         itr_pair = traits.construct_curves_2_object()(pgn);
  Traits_adapter_2 traits_adapter(traits);

  if(itr_pair.first == itr_pair.second)
    return true; // empty polygon

  return
    (traits_adapter.orientation_2_object()(itr_pair.first, itr_pair.second) ==
     COUNTERCLOCKWISE);
}

/* A valid polygon is :
 * 1 - Closed or empty polygon
 * 2 - Simple (previously known as strictly simple)
 * 3 - Counterclockwise oriented
 */
template <typename Traits_2>
bool is_valid_polygon(const typename Traits_2::Polygon_2& pgn,
                      const Traits_2& traits)
{
  bool closed = is_closed_polygon(pgn,traits);
  CGAL_warning_msg (closed, "The polygon's boundary is not closed.");
  if (! closed)
    return false;

  bool simple = is_simple_polygon(pgn,traits);
  CGAL_warning_msg (simple, "The polygon is not simple.");
  if (!simple)
    return false;

  bool valid_orientation = has_valid_orientation_polygon(pgn,traits);
  CGAL_warning_msg (valid_orientation,
                    "The polygon has a wrong orientation.");
  if (! valid_orientation)
    return false;

  return true;
}


template <typename Traits_2>
bool
is_closed_polygon_with_holes(const typename Traits_2::Polygon_with_holes_2& pgn,
                             const Traits_2& traits)
{
  typedef typename Traits_2::Polygon_with_holes_2       Polygon_with_holes_2;
  if (! is_closed_polygon(pgn.outer_boundary(),traits)) return false;

  typename Polygon_with_holes_2::Hole_const_iterator    itr;
  for (itr = pgn.holes_begin(); itr != pgn.holes_end(); ++itr) {
    if (! is_closed_polygon(*itr,traits)) return false;
  }
  return true;
}

// templated point location version
template<class Traits_2, class PointLocation>
bool
is_crossover_outer_boundary(const typename Traits_2::Polygon_with_holes_2& pgn,
                            const Traits_2& traits, PointLocation& pl)
{
  typedef typename Traits_2::Curve_const_iterator    Curve_const_iterator;
  typedef std::pair<Curve_const_iterator,Curve_const_iterator>    Cci_pair;

  typedef typename Traits_2::Point_2                 Point_2;
  typedef typename Traits_2::Compare_endpoints_xy_2  Compare_endpoints_xy_2;
  typedef typename Traits_2::Construct_min_vertex_2  Construct_min_vertex_2;
  typedef typename Traits_2::Construct_max_vertex_2  Construct_max_vertex_2;
  typedef CGAL::Gps_default_dcel<Traits_2>           Dcel;

  // IMPORTATNT! TODO!
  // Currently the topology traits is the bounded planar traits. This
  // should be replaced with a templated topology traits!
  typedef typename Default_planar_topology<Traits_2, Dcel>::Traits
                                                     Topology_traits;
  typedef CGAL::Gps_on_surface_base_2<Traits_2, Topology_traits>
    Polygon_set_2;

  typedef typename Polygon_set_2::Arrangement_on_surface_2 Arrangement_2;
  typedef typename Arrangement_2::Halfedge_handle          Halfedge_handle;
  typedef typename Arrangement_2::Vertex_handle            Vertex_handle;
  typedef typename Arrangement_2::Vertex_const_handle      Vertex_const_handle;
  typedef typename Traits_2::Curve_const_iterator          Curve_const_iterator;

  typename std::list<Halfedge_handle>           he_path;
  typename std::list<Halfedge_handle>::iterator he_itr;
  //functors used throughout the function
  Construct_min_vertex_2 min_functor = traits.construct_min_vertex_2_object();
  Construct_max_vertex_2 max_functor = traits.construct_max_vertex_2_object();
  Compare_endpoints_xy_2 cmp_endpoints =  traits.compare_endpoints_xy_2_object();

  Cci_pair itr_pair = traits.construct_curves_2_object()(pgn.outer_boundary());
  Curve_const_iterator  begin = itr_pair.first;
  Curve_const_iterator  end = itr_pair.second;
  if (begin == end) return true;  // An empty polygon is valid.
  // handles to consecutive curves
  Curve_const_iterator curr, next;
  curr = next = begin;
  // handles to vertices for insert. one maintains the current curve (already
  // inserted) and next curve's joint vertex.
  // the other maintains the next curve's second vertex if it already exists in
  // the arrangement.
  Vertex_handle joint_ver, second_ver;
  // closed check guarantees polygon has more than 1 curve
  ++next;
  // halfedge handle whose target is always the joint vertex between next and
  // curr.
  Halfedge_handle last_he;

  Polygon_set_2 gps(traits);
  Arrangement_2& arr = gps.arrangement();
  pl.attach(arr);

  // insert first edge lexicographically to arrangement
  // compute the joint vertex and insert to the path list a halfedge whose
  // target is the joint vertex
  last_he = CGAL::insert_non_intersecting_curve(arr, *curr);
  if  (cmp_endpoints(*curr) == SMALLER) {
    // polygon's boundary first curve is in lexicographic direction
    joint_ver = last_he->target();
    he_path.push_back(last_he);
  }
  else { // polygon's boundary first curve not lexicographic
    joint_ver = last_he->source();
    he_path.push_back(last_he->twin());
  }

  /* insert the rest of the curves to the arrangement efficiently the previous
   * closed polygon check guarantees equal_func
   * (construct_vertex_func (*curr, 1), construct_vertex_func (*next, 0)))
   */
  while (next != end) {
    CGAL::Object obj;
    Vertex_const_handle cver;
    Point_2 second_point;
    if (cmp_endpoints(*next) == SMALLER) {
      // next curve's minimum is the joint vertex. Look if it's max exists in
      // the arrangement and insert lexicographically
      second_point = max_functor(*next);
      obj = pl.locate(second_point);
      if (CGAL::assign (cver, obj)) {
        // insert where both vertices exist
        second_ver = arr.non_const_handle(cver);
        last_he = arr.insert_at_vertices( *next, joint_ver, second_ver);
      }
      else // insert from left vertex
        last_he = arr.insert_from_left_vertex ( *next,joint_ver) ;
    } else {
      // next curve's maximum vertex is the joint vertex. try to locate the
      // min vertex, and insert from right or from both vertices
      second_point = min_functor(*next);
      obj = pl.locate(second_point);
      if (CGAL::assign (cver, obj))  {
        // insert where both vertices exist
        second_ver = arr.non_const_handle(cver);
        last_he = arr.insert_at_vertices( *next, joint_ver, second_ver);
      }
      else  // insert from right vertex
        last_he = arr.insert_from_right_vertex ( *next,joint_ver) ;
    }
    // Move to the next pair of edges.
    he_path.push_back(last_he);
    joint_ver=last_he->target();
    curr = next;
    ++next;
  } //end of while

  /* We created a path of halfedges that circulates the polygon
   * counterclockwise. The polygon should lay on the left of each of these
   * half edges. If the boundary is invalid, the unbounded face should be
   * on the left of one of more than one of the halfedges.
   * The unbounded face is always to the right of the halfedges. We check if
   * all faces that lay on the right of the halfedges are equal (to the
   *"unbounded" face).
   */
  typename Arrangement_2::Face_handle fh = (*he_path.begin())->twin()->face();
  for (he_itr = he_path.begin(); he_itr != he_path.end(); he_itr++) {
    if ((*he_itr)->twin()->face() != fh)
      return false;
  }
  return true;
}

template<typename Traits_2>
bool is_crossover_outer_boundary
(const typename Traits_2::Polygon_with_holes_2& pgn, const Traits_2& traits)
{

  typedef CGAL::Gps_default_dcel<Traits_2>                      Dcel;
  // IMPORTATNT! TODO!
  // Currently the topology traits is the bounded planar traits. This
  // should be replaced with a templated topology traits!
  typedef typename Default_planar_topology<Traits_2, Dcel>::Traits
                                                                Topology_traits;

  typedef CGAL::Gps_on_surface_base_2<Traits_2, Topology_traits>
                                                                Polygon_set_2;
  typedef typename Polygon_set_2::Arrangement_on_surface_2      Arrangement_2;
  typedef CGAL::Arr_naive_point_location<Arrangement_2>         Naive_pl;

  Naive_pl pl;
  return is_crossover_outer_boundary(pgn, traits, pl);
}

// previously known as Simple
template <typename Traits_2>
bool is_relatively_simple_polygon_with_holes
(const typename Traits_2::Polygon_with_holes_2& pgn, const Traits_2& traits)
{
  typedef typename Traits_2::Curve_const_iterator    Curve_const_iterator;
  typedef std::pair<Curve_const_iterator,Curve_const_iterator>    Cci_pair;
  typedef typename Traits_2::Construct_curves_2    Construct_curves_2;

  typedef typename Traits_2::X_monotone_curve_2     X_monotone_curve_2;
  typedef Gps_polygon_validation_visitor<Traits_2>  Visitor;
  typedef Sweep_line_2<Traits_2, Visitor>           Sweep_line;
  typedef typename Traits_2::Polygon_with_holes_2   Polygon_with_holes_2;

  Construct_curves_2 construct_curves_func = traits.construct_curves_2_object();
  // Construct a container of all outer boundary curves.
  Cci_pair         itr_pair = construct_curves_func (pgn.outer_boundary());
  std::list<X_monotone_curve_2>  outer_curves;
  std::copy (itr_pair.first, itr_pair.second,
             std::back_inserter(outer_curves));
  // Create visitor and sweep to verify outer boundary is relatively simple
  Visitor      relative_visitor(false);
  Sweep_line   sweep_line (&traits, &relative_visitor);
  relative_visitor.sweep (outer_curves.begin(), outer_curves.end());
  if (!relative_visitor.is_valid()) {
    switch (relative_visitor.error_code()) {
     case Visitor::ERROR_NONE: break;
     case Visitor::ERROR_EDGE_INTERSECTION:
      CGAL_warning_msg(false, "The outer boundary self intersects at edges.");
      break;

     case Visitor::ERROR_EDGE_VERTEX_INTERSECTION:
      CGAL_warning_msg(false, "The outer boundary self (weakly) intersects.");
      break;

     case Visitor::ERROR_EDGE_OVERLAP:
      CGAL_warning_msg(false, "The outer boundary self overlaps.");
      break;

     case Visitor::ERROR_VERTEX_INTERSECTION:
      CGAL_warning_msg(false, "The outer boundary self intersects at vertices.");
      break;
    }
    return false;
  }

  // Verify every hole is simple
  typename Polygon_with_holes_2::Hole_const_iterator  hoit;
  std::list<X_monotone_curve_2>  hole_curves;
  for (hoit = pgn.holes_begin(); hoit != pgn.holes_end(); ++hoit) {
    bool simple_hole = is_simple_polygon(*hoit, traits);
    if (!simple_hole)
      return false;
  }
  return true;
}

template <typename Traits_2>
bool has_valid_orientation_polygon_with_holes
(const typename Traits_2::Polygon_with_holes_2& pgn, const Traits_2& traits)
{
  typedef Gps_traits_adaptor<Traits_2>                  Traits_adapter_2;
  typedef typename Traits_2::Curve_const_iterator       Curve_const_iterator;
  typedef std::pair<Curve_const_iterator,Curve_const_iterator>    Cci_pair;
  typedef typename Traits_2::Construct_curves_2         Construct_curves_2;

  typedef typename Traits_adapter_2::Orientation_2      Check_orientation_2;
  typedef typename Traits_2::Polygon_with_holes_2       Polygon_with_holes_2;

  Traits_adapter_2 traits_adapter(traits);

  Construct_curves_2 construct_curves_func = traits.construct_curves_2_object();
  Check_orientation_2 check_orientation_func =
    traits_adapter.orientation_2_object();
  // Check the orientation of the outer boundary.
  Cci_pair itr_pair = construct_curves_func (pgn.outer_boundary());

  if ((itr_pair.first != itr_pair.second) &&
      (check_orientation_func (itr_pair.first, itr_pair.second) !=
       COUNTERCLOCKWISE))
  {
    return false;
  }

  // Check the orientation of each of the holes.
  typename Polygon_with_holes_2::Hole_const_iterator    hoit;

  for (hoit = pgn.holes_begin(); hoit != pgn.holes_end(); ++hoit) {
    itr_pair = construct_curves_func (*hoit);

    if ((itr_pair.first != itr_pair.second) &&
        (check_orientation_func (itr_pair.first, itr_pair.second) != CLOCKWISE))
    {
      return false;
    }
  }
  return true;
}

/* Verify holes do not intersect between themselves as well with the outer
 * boundary (except intersection on a vertex which is allowed).
 *
 * This efficient implementation utilizes the general poygon set for aggregated
 * join operations for N holes which should result in a GPS that contains N
 * independent PWH.
 * Executing a difference(gps, outer boundary) should result in an empty set if
 * no holes intersect the boundary.
 *
 * An iterative use of the difference free function while iterating over the
 * holes may have an advantage in case there are numerous holes that intersect
 * the boundary and the iterative loop will be stopped after a small number of
 * iterations.
 */
template <class Traits_2>
bool are_holes_and_boundary_pairwise_disjoint
(const typename Traits_2::Polygon_with_holes_2& pwh, Traits_2& traits)
{
  typedef typename Traits_2::Curve_const_iterator    Curve_const_iterator;
  typedef std::pair<Curve_const_iterator,Curve_const_iterator>    Cci_pair;
  typedef typename Traits_2::Construct_curves_2    Construct_curves_2;

  typedef CGAL::Gps_default_dcel<Traits_2>                 Dcel;
  // IMPORTATNT! TODO!
  // Currently the topology traits is the bounded planar traits. This
  // should be replaced with a templated topology traits!
  typedef typename Default_planar_topology<Traits_2, Dcel>::Traits
                                                           Topology_traits;

  typedef CGAL::Gps_on_surface_base_2<Traits_2, Topology_traits>
    Polygon_set_2;
  typedef typename Polygon_set_2::Size                     Size;
  typedef  typename Traits_2::Polygon_2                    Polygon_2;
  typedef typename Traits_2::Polygon_with_holes_2          Polygon_with_holes_2;
  typedef typename Polygon_with_holes_2::Hole_const_iterator
    Hole_const_iterator;
  typedef typename Traits_2::X_monotone_curve_2            X_monotone_curve_2;
  typedef std::pair<Curve_const_iterator,Curve_const_iterator>
                                                           Cci_pair;
  typedef typename Traits_2::Construct_curves_2            Construct_curves_2;
  typedef typename Traits_2::Construct_general_polygon_with_holes_2
    Construct_polygon_with_holes_2;

  typedef Gps_polygon_validation_visitor<Traits_2>         Visitor;
  typedef Sweep_line_2<Traits_2, Visitor>                  Sweep_line ;
  typedef typename Polygon_set_2::Arrangement_on_surface_2 Arrangement_2;

  /* Should be perfored more efficeintly  than using sweep and than
   * difference().
   *
   * Use sweep to find intersections on the interior of curves (not on vertices)
   * and overlapping edges which are not allowed (note that 0/1 dimension
   * intersections are not detectes by do_intersect() which only returns the
   * 2D intersection polygon if exists)
   * Note that using this sweep alone allows for a hole and an edge to share
   * a vertex and intersect (like illegal input pgn_w_overlap_hole.dat in
   * validation_example)
   */
  Hole_const_iterator hoit;
  // Construct a container of all boundary curves.
  Polygon_2 pgn2 = traits.construct_outer_boundary_object()(pwh);
  Construct_curves_2    construct_curves_func;
  Cci_pair itr_pair = construct_curves_func(pgn2);

  std::list<X_monotone_curve_2>  curves;
  std::copy(itr_pair.first, itr_pair.second, std::back_inserter(curves));

  std::pair<Hole_const_iterator, Hole_const_iterator> pair =
    traits.construct_holes_object()(pwh);
  //for (hoit = pgn.holes_begin(); hoit != pgn.holes_end(); ++hoit)
  for (hoit = pair.first; hoit!=pair.second; ++hoit) {
    itr_pair = construct_curves_func (*hoit);
    std::copy (itr_pair.first, itr_pair.second, std::back_inserter(curves));
  }

  // Perform the sweep and check for curve  intersections on the interior.
  // Traits_2     traits; moved to top, needed also for boundary.
  Visitor visitor(false);
  Sweep_line sweep_line(&traits, &visitor);
  visitor.sweep(curves.begin(), curves.end());
  if (!visitor.is_valid()) return false;

  Polygon_set_2 gps(traits);
  // check for 2D  intersections of holes (holes must be disjoint except for
  // vertices)
  Size num_of_holes = 0;
  // functors for creating a pwh needed for inserting pgns into the arrangement
  // quickly
  Construct_polygon_with_holes_2 construct_pwh_functor =
    traits.construct_polygon_with_holes_2_object() ;
  for (hoit = pwh.holes_begin(); hoit != pwh.holes_end(); ++hoit) {
    Polygon_2 hole(*hoit);
    hole.reverse_orientation();
    /* gps.join() and gps.insert()requires that the polyon insrted is valid,
     * and therfore hole orientation must be reversed
     */
    bool intersect = gps.do_intersect(hole);
    if (intersect) return false;
    else {
      /* to use gps.insert(hole) it is required that the set coponents and the
       * new holes  do not intersect.
       * because the sweep detects shared edges and the do_intersect query
       * detects 2D intersections we can safely use the insert(pwh) function
       * whose performance is better than the join(pgn)
       */
      Polygon_with_holes_2 empty_pwh = construct_pwh_functor(hole);
      // traits.Construct_general_polygon_with_holes_2 (hole);
      // Polygon_with_holes_2 empty_pwh(hole);
      gps.insert(empty_pwh);
      num_of_holes++;
    }
  }
  /* not good - doesn't work if intersection at vertices is legal.
   * Size arr_num_of_holes = gps.number_of_polygons_with_holes();
   * if (num_of_holes != arr_num_of_holes)
   *   return false;
   */

  // check for intersection of holes with the outer boundary

  /* outer boundary can be relatively simple. Execution of
   * do_intersect(hole, boundary) or difference(hole,boundary) relies on
   * implementation of General polygon set which has a precondition that
   * requires valid polygon or PWH to be inserted (not just a simple polygon).
   * This helper function is utilized after checking for the PWH closure,
   * relative simplicity and orientation. Therefore it is safe to assume the
   * outer boundary is  valid PWH with no holes. We can't assume it is a valid
   * (simple) polygon.
   */

  //Polygon_with_holes_2 boundary(pwh.outer_boundary(), fit, fit);
 Polygon_with_holes_2 boundary =  construct_pwh_functor(pwh.outer_boundary());
  // Unbounded outer boundaries contain all the holes and the holes were checked
  // and are OK.
  if (boundary.is_unbounded()) return true;

  /* do_intersect predicate will not suffice as hole can be completely outside
   * the outer boundary in an (extremely strange) case
   * The gps now contains all the holes. the difference between the boundary
   * and a union of all the holes should be the empty set. For performance
   * reasons, we use a customized overlay traits and perform an arrangement
   * overlay instead of difference
   */
  ValidationOverlayTraits<Arrangement_2> valOverlayTraits;
  valOverlayTraits.setHoleOverlap(false);
  Polygon_set_2 gps2(traits);

  Arrangement_2& boundary_arr = gps2.arrangement();
  gps2._insert(boundary,boundary_arr);
  Arrangement_2& holes_arr = gps.arrangement();
  Arrangement_2 output_arr(holes_arr.geometry_traits());
  overlay(holes_arr, boundary_arr, output_arr, valOverlayTraits);
  if (valOverlayTraits.getHoleOverlap()) return false;

  /* old code that works less efficiently than the new overly traits
   * gps.validation_difference(boundary);
   * if gps is not empty at least one hole intersected the boundary
   * if (!gps.is_empty())
   *   return false;
   */
  return true;
}

/* A valid polygon with holes is :
 * 1 - Has empty or closed boundary and all the holes are closed
 * 2 - The PWH is relatively simple polygon (holes are simple...)
 * 3 - Has it's boundary oriented counterclockwise and the holes oriented
 *     clockwise
 * 4 - All the segments (boundry and holes) do not cross or intersect in their
 *     relative interior
 * 5 - The holes are on the interior of the boundary polygon if the boundary
 *     is not empty
 */
template <typename Traits_2>
bool
is_valid_polygon_with_holes(const typename Traits_2::Polygon_with_holes_2& pgn,
                            const Traits_2& traits)
{
  bool closed = is_closed_polygon_with_holes(pgn, traits);
  CGAL_warning_msg(closed,
                   "The polygon's boundary or one of its holes is not closed.");
  if (! closed) return false;

  bool relatively_simple = is_relatively_simple_polygon_with_holes(pgn, traits);
  CGAL_warning_msg (relatively_simple, "The polygon is not relatively simple.");
  if (! relatively_simple) return false;

  bool no_cross = is_crossover_outer_boundary(pgn, traits);
  CGAL_warning_msg (no_cross, "The polygon has a crossover.");
  if (!no_cross) return false;

  bool valid_orientation = has_valid_orientation_polygon_with_holes(pgn, traits);
  CGAL_warning_msg (valid_orientation, "The polygon has a wrong orientation.");
  if (! valid_orientation) return false;

  bool holes_disjoint = are_holes_and_boundary_pairwise_disjoint(pgn, traits);
  CGAL_warning_msg
    (holes_disjoint,
     "Holes of the PWH intersect amongst themselves or with outer boundary");
  if (! holes_disjoint) return false;

  return true;
}

template <typename Traits_2>
bool
is_valid_unknown_polygon(const typename Traits_2::Polygon_with_holes_2& pgn,
                         const Traits_2& traits)
{ return is_valid_polygon_with_holes(pgn, traits); }

template <typename Traits_2>
bool is_valid_unknown_polygon(const typename Traits_2::Polygon_2& pgn,
                              const Traits_2& traits)
{ return is_valid_polygon(pgn, traits); }

} //namespace CGAL

#endif
