// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//                 Ophir Setter    <ophir.setter@cs.tau.ac.il>
//                 Guy Zucker <guyzucke@post.tau.ac.il>

#ifndef CGAL_GPS_ON_SURFACE_BASE_2_IMPL_H
#define CGAL_GPS_ON_SURFACE_BASE_2_IMPL_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/iterator.h>
#include <CGAL/function_objects.h>
#include <CGAL/circulator.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_scanner.h>
#include <CGAL/Arr_accessor.h>

#include <queue>
#include <list>

namespace CGAL {

template <class Traits_, class TopTraits_, class ValidationPolicy>
void Gps_on_surface_base_2<Traits_, TopTraits_,ValidationPolicy>::
construct_polygon(Ccb_halfedge_const_circulator ccb, Polygon_2 & pgn,
                  const Traits_* tr)
{
  typedef CGAL::Ccb_curve_iterator<Arrangement_on_surface_2>
    Ccb_curve_iterator;
  Ccb_curve_iterator begin(ccb, false);
  Ccb_curve_iterator end(ccb, true);

  tr->construct_polygon_2_object()(begin, end, pgn);
}

// The comments below was written after trying to understand what the visitors
// do. There was no comment by the author of this class.
// This class is used afterwards to extract polygons from the representing
// arrangement.
// This scanner is not the same as the Gps_bfs_scanner. In this file, the
// Gps_bfs_scanner is used with Init_faces_visitor to init the faces of the
// representing arrangement.
// It seems that Gps_bfs_scanner is used for a regular bfs scan on the faces
// of the arrangements, with comparison to Arr_bfs_scanner that cares about
// inner ccbs and outer ccbs (it treats them differently).
// If this is the case, we should unite Gps_bfs_scanner with the regular
// adaptation of arrangement to boost graph.
template <class Arrangement, class OutputIterator>
class Arr_bfs_scanner
{
public:
  typedef typename Arrangement::Geometry_traits_2       Gps_traits;
  typedef typename Arrangement::Topology_traits         Gps_top_traits;
  typedef typename Gps_traits::Polygon_2                Polygon_2;
  typedef typename Gps_traits::Polygon_with_holes_2     Polygon_with_holes_2;
  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement::Face_const_iterator     Face_const_iterator;
  typedef typename Arrangement::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Arrangement::Outer_ccb_const_iterator
    Outer_ccb_const_iterator;
  typedef typename Arrangement::Inner_ccb_const_iterator
    Inner_ccb_const_iterator;


protected:
  const Gps_traits* m_traits;
  std::queue<Face_const_iterator> m_holes_q;
  std::list<Polygon_2> m_pgn_holes;
  OutputIterator m_oi;

public:
  /*! Constructor */
  Arr_bfs_scanner(const Gps_traits* tr, OutputIterator oi) :
    m_traits(tr),
    m_oi(oi)
  {}

  void scan(Arrangement& arr)
  {
    Face_const_iterator   ubf;
    for (ubf = arr.faces_begin(); ubf != arr.faces_end(); ++ubf)
    {
      if (ubf->number_of_outer_ccbs() != 0)
        continue;
      if (ubf->visited())
        continue;

      Inner_ccb_const_iterator  holes_it;
      if (!ubf->contained())
      {
        ubf->set_visited(true);
        for (holes_it = ubf->inner_ccbs_begin();
             holes_it != ubf->inner_ccbs_end(); ++holes_it)
        {
          scan_ccb (*holes_it);
        }
      }
      else
      {
        // ubf is contained -> unbounded polygon !!
        scan_contained_ubf(ubf);

      }

      while(!m_holes_q.empty())
      {
        Face_const_iterator top_f = m_holes_q.front();
        m_holes_q.pop();
        top_f->set_visited(true);
        for (holes_it = top_f->inner_ccbs_begin();
             holes_it != top_f->inner_ccbs_end(); ++holes_it)
        {
          scan_ccb(*holes_it);
        }

        //scan_uncontained_face(top_f->outer_ccb());
      }
    }

    Face_const_iterator   fit;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
    {
      fit->set_visited(false);
    }
  }

  OutputIterator output_iterator() const
  {
    return m_oi;
  }

  void scan_ccb(Ccb_halfedge_const_circulator ccb)
  {

    Polygon_2 pgn_boundary;
    Gps_on_surface_base_2<Gps_traits, Gps_top_traits>::
      construct_polygon(ccb, pgn_boundary, m_traits);

    Ccb_halfedge_const_circulator ccb_end = ccb;
    do
    {
      Halfedge_const_iterator he = ccb;
      if (!he->twin()->face()->visited())
        all_incident_faces(he->twin()->face());
      ++ccb;
    }
    while(ccb != ccb_end);
    Polygon_with_holes_2 pgn =
      m_traits->construct_polygon_with_holes_2_object()(pgn_boundary,
                                                        m_pgn_holes.begin(),
                                                        m_pgn_holes.end());
    /*Polygon_with_holes_2 pgn(pgn_boundary,
                             m_pgn_holes.begin(),
                             m_pgn_holes.end());*/
    *m_oi = pgn;
    ++m_oi;
    m_pgn_holes.clear();
  }

  void scan_contained_ubf(Face_const_iterator ubf)
  {
    CGAL_assertion(ubf->number_of_outer_ccbs() == 0 && ubf->contained());
    // ubf is contained -> unbounded polygon !!
    all_incident_faces(ubf);
    Polygon_2 boundary;
    Polygon_with_holes_2 pgn =
      m_traits->construct_polygon_with_holes_2_object()(boundary,
                                                        m_pgn_holes.begin(),
                                                        m_pgn_holes.end());
    /*Polygon_with_holes_2 pgn(boundary,
                             m_pgn_holes.begin(),
                             m_pgn_holes.end());*/
    *m_oi = pgn;
    ++m_oi;
    m_pgn_holes.clear();
  }


  void all_incident_faces(Face_const_iterator f)
  {
    CGAL_assertion(!f->visited());
    f->set_visited(true);
    if (f->number_of_outer_ccbs() != 0)
    {
      if (!f->contained())
      {
        for (Outer_ccb_const_iterator oci = f->outer_ccbs_begin();
             oci != f->outer_ccbs_end(); ++oci)
        {
          m_pgn_holes.push_back(Polygon_2());
          Gps_on_surface_base_2<Gps_traits, Gps_top_traits>::
            construct_polygon(*oci, m_pgn_holes.back(), m_traits);
        }

        m_holes_q.push(f);
      }


      for (Outer_ccb_const_iterator oci = f->outer_ccbs_begin();
           oci != f->outer_ccbs_end(); ++oci)
      {
        Ccb_halfedge_const_circulator ccb_end = *oci;
        Ccb_halfedge_const_circulator ccb_circ = ccb_end;
        do
        {
          //get the current halfedge on the face boundary
          Halfedge_const_iterator he  = ccb_circ;
          Face_const_iterator new_f = he->twin()->face();
          if (!new_f->visited())
          {
            all_incident_faces(new_f);
          }
          ++ccb_circ;
        }
        while(ccb_circ != ccb_end);
      }
    }

    if (f->contained())
    {
      Inner_ccb_const_iterator hit;
      for(hit = f->inner_ccbs_begin(); hit != f->inner_ccbs_end(); ++hit)
      {
        Ccb_halfedge_const_circulator ccb_of_hole = *hit;
        Halfedge_const_iterator he = ccb_of_hole;
        if (is_single_face(ccb_of_hole))
        {
          CGAL_assertion(!he->twin()->face()->contained());

          m_pgn_holes.push_back(Polygon_2());
          Gps_on_surface_base_2<Gps_traits, Gps_top_traits>::
            construct_polygon(he->twin()->face()->outer_ccb(),
                              m_pgn_holes.back(), m_traits);
          m_holes_q.push(he->twin()->face());
        }
        else
        {
          Ccb_halfedge_const_circulator ccb_end = ccb_of_hole;
          do
          {
            Halfedge_const_iterator he = ccb_of_hole;
            if (!he->twin()->face()->visited())
              all_incident_faces(he->twin()->face());
            ++ccb_of_hole;
          }
          while(ccb_of_hole != ccb_end);
        }
      }
    }
  }

  bool is_single_face(Ccb_halfedge_const_circulator ccb)
  {
    Ccb_halfedge_const_circulator ccb_end = ccb;
    Ccb_halfedge_const_circulator ccb_circ = ccb_end;
    Halfedge_const_iterator he = ccb;
    Face_const_iterator curr_f = he->twin()->face();
    do
    {
      //get the current halfedge on the face boundary
      Halfedge_const_iterator he  = ccb_circ;
      if (he->twin()->face() != curr_f)
        return false;
      if (he->twin()->target()->degree() != 2)
        return false;
      ++ccb_circ;
    }
    while(ccb_circ != ccb_end);
    return true;
  }
};


template <class Arrangement>
class Init_faces_visitor
{
  typedef typename Arrangement::Face_iterator             Face_iterator;
  typedef typename Arrangement::Halfedge_iterator         Halfedge_iterator;

public:

  //! discovered_face
/*! discovered_face is called by Gps_bfs_scanner when it reveals a new face
    during a BFS scan. It is important to say that I have a strong suspition
    that this place is the reason why discovered_face was once called
    "flip_face" (WTF?)
  \param old_f The face that was already revealed
  \param new_f The face that we have just now revealed
*/
  void discovered_face(Face_iterator old_f,
                       Face_iterator new_f,
                       Halfedge_iterator /*he*/)
  {
    new_f->set_contained(!old_f->contained());
  }
};

//! _insert
/*! The function inserts a polygon into an arrangement, assuming that the
    polygon is contained in one face of the arrangement.
  \param pgn The polygon to be inserted to the arrangement. pgn must be
             completely disjoint from the arrangement
  \param arr The arrangement to insert the polygon to.
*/
template <class Traits_, class TopTraits_, class ValidationPolicy>
void Gps_on_surface_base_2<Traits_, TopTraits_,ValidationPolicy>::
_insert(const Polygon_2& pgn, Arrangement_on_surface_2 & arr)
{
  typedef Arr_accessor<Arrangement_on_surface_2>                  Arr_accessor;

  Arr_accessor  accessor(arr);
  Compare_endpoints_xy_2  cmp_ends = m_traits->compare_endpoints_xy_2_object();

  std::pair<Curve_const_iterator, Curve_const_iterator> itr_pair =
    m_traits->construct_curves_2_object()(pgn);

  if (itr_pair.first == itr_pair.second)
    return;

  Curve_const_iterator curr = itr_pair.first;
  Curve_const_iterator end  = itr_pair.second;

  const Arr_parameter_space  ps_x =
    m_traits_adaptor.parameter_space_in_x_2_object()(*curr, ARR_MIN_END);
  const Arr_parameter_space  ps_y =
    m_traits_adaptor.parameter_space_in_y_2_object()(*curr, ARR_MIN_END);

  Object obj_f;
  if ((ps_x == ARR_INTERIOR) && (ps_y == ARR_INTERIOR))
  {
    Point_location pl(arr);
    obj_f = pl.locate(m_traits->construct_min_vertex_2_object()(*curr));
  }
  else
  {
    obj_f = accessor.locate_curve_end(*curr, ARR_MIN_END, ps_x, ps_y);
  }

  Face_const_handle const_f;
  // face should not be contained as the pgn is completly disjoint of the
  // arrangement.
  CGAL_assertion(CGAL::assign(const_f, obj_f) && !const_f->contained());
  CGAL::assign(const_f, obj_f);
  Face_iterator f = arr.non_const_handle(const_f);

  Halfedge_handle first_he =
    arr.insert_in_face_interior(*curr, f);
  //first_he is directed from left to right (see insert_in_face_interior)

  Halfedge_handle curr_he;
  if (cmp_ends(*curr) == CGAL::SMALLER)
  {
    // curr curve and first_he have the same direction
    curr_he = first_he;
    first_he = first_he->twin();
  }
  else
  {
    // curr curve and first_he have opposite directions
    CGAL_assertion(cmp_ends(*curr) == CGAL::LARGER);
    curr_he = first_he->twin();
  }

  Curve_const_iterator temp = curr;
  ++temp;
  if (temp == end) // a polygon with circular arcs may have only
                  // two edges (full circle for example)
  {
    /*Halfedge_handle he =
      arr.insert_at_vertices(*temp, curr_he, first_he);*/
    bool new_face_created = false;
    bool dummy_swapped_predecessors = false;
    Halfedge_handle he = accessor.insert_at_vertices_ex (curr_he,
                                                         *temp, (cmp_ends(*temp) == CGAL::SMALLER ? ARR_LEFT_TO_RIGHT : ARR_RIGHT_TO_LEFT),
                                                         first_he->next(),
                                                         new_face_created,
                                                         dummy_swapped_predecessors);
    // TODO EBEB 2012-08-06 do we have to care if order has been swapped,
    // or do we have to disallow swapping?

    CGAL_assertion(new_face_created);
    CGAL_assertion((he->face() != he->twin()->face()));

    he->face()->set_contained(true);
    return;
  }

  //The polygon has 3 or more edges
  Curve_const_iterator last = end;
  --last;
  for(++curr ; curr != last; ++curr)
  {
    const X_monotone_curve_2& curr_cv = *curr;
    if (cmp_ends(curr_cv) == CGAL::SMALLER)
      curr_he = arr.insert_from_left_vertex(curr_cv, curr_he);
    else
    {
      CGAL_assertion(cmp_ends(curr_cv) == CGAL::LARGER);
      curr_he = arr.insert_from_right_vertex(curr_cv, curr_he);
    }
  }

  const X_monotone_curve_2& last_cv = *last;
  /*Halfedge_handle last_he =
    arr.insert_at_vertices(last_cv, curr_he, first_he); */
  bool new_face_created = false;
  bool dummy_swapped_predecessors = false;
  Halfedge_handle last_he =
    accessor.insert_at_vertices_ex (curr_he,
                                    last_cv, ( cmp_ends(last_cv) == CGAL::SMALLER ? ARR_LEFT_TO_RIGHT : ARR_RIGHT_TO_LEFT),
                                    first_he->next(),
                                    new_face_created,
                                    dummy_swapped_predecessors);
  // TODO EBEB 2012-08-06 do we have to care if order has been swapped,
  // or do we have to disallow swapping?

  CGAL_assertion(new_face_created);
  CGAL_assertion((last_he->face() != last_he->twin()->face()));

  last_he->face()->set_contained(true);
}


template <class Traits_, class TopTraits_, class ValidationPolicy>
  template<class PolygonIter >
  void Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  insert(PolygonIter p_begin, PolygonIter p_end)
{
  typename std::iterator_traits<PolygonIter>::value_type pgn;
  //check validity of all polygons
  for(PolygonIter pitr = p_begin; pitr != p_end; ++pitr)
  {
    ValidationPolicy::is_valid(*pitr, *m_traits);
  }

  _insert(p_begin, p_end, pgn);
}

template <class Traits_, class TopTraits_, class ValidationPolicy>
  template<class PolygonIter, class PolygonWithHolesIter>
  void Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  insert(PolygonIter p_begin, PolygonIter p_end,
         PolygonWithHolesIter pwh_begin, PolygonWithHolesIter pwh_end)
{
  typedef std::list<X_monotone_curve_2>                  XCurveList;
  typedef Init_faces_visitor<Arrangement_on_surface_2>              My_visitor;
  typedef Gps_bfs_scanner<Arrangement_on_surface_2, My_visitor>     Arr_bfs_scanner;

  XCurveList xcurve_list;

  for( ; p_begin != p_end; ++p_begin)
  {
    ValidationPolicy::is_valid(*p_begin, *m_traits);
    _construct_curves(*p_begin, std::back_inserter(xcurve_list));
  }

  bool is_unbounded = false;
  for( ; pwh_begin != pwh_end; ++pwh_begin)
  {
    ValidationPolicy::is_valid(*pwh_begin, *m_traits);
    is_unbounded = (is_unbounded || m_traits->construct_is_unbounded_object()(*pwh_begin));
    // is_unbounded = (is_unbounded || pwh_begin->is_unbounded());
    _construct_curves(*pwh_begin, std::back_inserter(xcurve_list));
  }
  insert_non_intersecting_curves(*m_arr, xcurve_list.begin(), xcurve_list.end());

  if (is_unbounded)
  {
    for (Face_iterator fit = m_arr->faces_begin();
         fit != m_arr->faces_end(); ++fit)
    {
      if (fit->number_of_outer_ccbs() == 0)
        fit->set_contained(true);
    }
  }

  My_visitor v;
  Arr_bfs_scanner scanner(v);
  scanner.scan(*m_arr);
  _reset_faces(m_arr);
}

//insert a range of simple polygons to the arrangement
template <class Traits_, class TopTraits_, class ValidationPolicy>
  template<class PolygonIter>
  void Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  _insert(PolygonIter p_begin, PolygonIter p_end, Polygon_2 & /*pgn*/)
{
  for(PolygonIter pitr = p_begin; pitr != p_end; ++pitr)
  {
    this->_insert(*pitr, *m_arr);
  }
}

template <class Traits_, class TopTraits_, class ValidationPolicy>
  template<class PolygonIter>
  void Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  _insert(PolygonIter p_begin, PolygonIter p_end, Polygon_with_holes_2 & /*pgn*/)
{
  typedef std::list<X_monotone_curve_2>                  XCurveList;
  typedef Init_faces_visitor<Arrangement_on_surface_2>              My_visitor;
  typedef Gps_bfs_scanner<Arrangement_on_surface_2, My_visitor>     Arr_bfs_scanner;

  XCurveList xcurve_list;
  bool is_unbounded = false;
  for( ; p_begin != p_end; ++p_begin)
  {
    // is_unbounded = (is_unbounded || p_begin->is_unbounded());
    is_unbounded = (is_unbounded || m_traits->construct_is_unbounded_object()(*p_begin));
    _construct_curves(*p_begin, std::back_inserter(xcurve_list));

  }
  insert_non_intersecting_curves(*m_arr, xcurve_list.begin(), xcurve_list.end());

  if (is_unbounded)
  {
    for (Face_iterator fit = m_arr->faces_begin();
         fit != m_arr->faces_end(); ++fit)
    {
      if (fit->number_of_outer_ccbs() == 0)
        fit->set_contained(true);
    }
  }

  My_visitor v;
  Arr_bfs_scanner scanner(v);
  scanner.scan(*m_arr);
  _reset_faces(m_arr);
}

//insert non-sipmle poloygons with holes (non incident edges may have
// common vertex,  but they dont intersect at their interior
template <class Traits_, class TopTraits_, class ValidationPolicy>
  void Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  _insert(const Polygon_with_holes_2 & pgn, Arrangement_on_surface_2 & arr)
{
 // inner function not exposed to user - no validation
 // ValidationPolicy::is_valid(pgn, *m_traits);

  typedef std::list<X_monotone_curve_2>                  XCurveList;
  typedef Init_faces_visitor<Arrangement_on_surface_2>          My_visitor;
  typedef Gps_bfs_scanner<Arrangement_on_surface_2, My_visitor> Arr_bfs_scanner;

  XCurveList xcurve_list;
  _construct_curves(pgn, std::back_inserter(xcurve_list));
  insert_non_intersecting_curves(arr, xcurve_list.begin(), xcurve_list.end());

  //if (pgn.is_unbounded())
  if (m_traits->construct_is_unbounded_object()(pgn))
  {
    for (Face_iterator fit = arr.faces_begin();
         fit != arr.faces_end(); ++fit)
    {
      if (fit->number_of_outer_ccbs() == 0)
        fit->set_contained(true);
    }
  }

  My_visitor v;
  Arr_bfs_scanner scanner(v);
  scanner.scan(arr);
  _reset_faces(&arr);
}

template <class Traits_, class TopTraits_, class ValidationPolicy>
  template <class OutputIterator>
  void
  Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  _construct_curves(const Polygon_2 & pgn, OutputIterator oi)
{
  std::pair<Curve_const_iterator, Curve_const_iterator> itr_pair =
    m_traits->construct_curves_2_object()(pgn);
  std::copy (itr_pair.first, itr_pair.second, oi);
}

template <class Traits_, class TopTraits_, class ValidationPolicy>
  template <class OutputIterator>
  void Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  _construct_curves(const Polygon_with_holes_2 & pgn, OutputIterator oi)
{
  //if (!pgn.is_unbounded())
  if (!m_traits->construct_is_unbounded_object()(pgn))
  {
    const Polygon_2& pgn_boundary =
      m_traits->construct_outer_boundary_object()(pgn);
    std::pair<Curve_const_iterator, Curve_const_iterator> itr_pair =
      m_traits->construct_curves_2_object()(pgn_boundary);
    std::copy (itr_pair.first, itr_pair.second, oi);
  }
  std::pair<GP_Holes_const_iterator, GP_Holes_const_iterator> hpair =
    m_traits->construct_holes_object()(pgn);
  GP_Holes_const_iterator hit;
  for (hit = hpair.first; hit != hpair.second; ++hit)
  {
    const Polygon_2& pgn_hole = *hit;
    std::pair<Curve_const_iterator, Curve_const_iterator> itr_pair =
      m_traits->construct_curves_2_object()(pgn_hole);
    std::copy (itr_pair.first, itr_pair.second, oi);
  }
}

template <class Traits_, class TopTraits_, class ValidationPolicy>
  template <class OutputIterator>
  OutputIterator
  Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  polygons_with_holes(OutputIterator out) const
{
  typedef Arr_bfs_scanner<Arrangement_on_surface_2, OutputIterator>     Arr_bfs_scanner;
  Arr_bfs_scanner scanner(this->m_traits, out);
  scanner.scan(*(this->m_arr));
  return (scanner.output_iterator());
}


template <class Traits_, class TopTraits_, class ValidationPolicy>
  typename Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::Size
  Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  number_of_polygons_with_holes() const
{

  typedef Arr_bfs_scanner<Arrangement_on_surface_2, Counting_output_iterator>
    Arr_bfs_scanner;
  //counting_output_operator CTOR reqires a parameter
  std::size_t cc = 0;
  Arr_bfs_scanner scanner(this->m_traits, Counting_output_iterator(&cc));
  scanner.scan(*(this->m_arr));
  return (scanner.output_iterator().current_counter());
}


template <class Traits_, class TopTraits_, class ValidationPolicy>
  bool Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  locate(const Point_2& q, Polygon_with_holes_2& pgn) const
{
  Point_location pl(*m_arr);

  Object obj = pl.locate(q);
  Face_const_iterator f;
  if (CGAL::assign(f, obj))
  {
    if (!f->contained())
      return false;
  }
  else
  {
    Halfedge_const_handle he;
    if (CGAL::assign(he, obj))
    {
      if (he->face()->contained())
        f = he->face();
      else
      {
        CGAL_assertion(he->twin()->face()->contained());
        f = he->twin()->face();
      }
    }
    else
    {
      Vertex_const_handle v;
      CGAL_assertion(CGAL::assign(v, obj));
      CGAL::assign(v, obj);
      Halfedge_around_vertex_const_circulator hav = v->incident_halfedges();
      Halfedge_const_handle he = hav;
      if (he->face()->contained())
        f = he->face();
      else
      {
        CGAL_assertion(he->twin()->face()->contained());
        f = he->twin()->face();
      }
    }
  }

  typedef Oneset_iterator<Polygon_with_holes_2>    OutputItr;
  typedef Arr_bfs_scanner<Arrangement_on_surface_2, OutputItr>     Arr_bfs_scanner;

  OutputItr oi (pgn);
  Arr_bfs_scanner scanner(this->m_traits, oi);


  Ccb_halfedge_const_circulator ccb_of_pgn = get_boundary_of_polygon(f);
  this->_reset_faces();
  if (ccb_of_pgn == Ccb_halfedge_const_circulator())
  {
    // the polygon has no boundary

    // f is unbounded
    for (Face_iterator fit = m_arr->faces_begin(); fit != m_arr->faces_end();
         ++fit)
    {
      if (fit->number_of_outer_ccbs() == 0)
        scanner.scan_contained_ubf(fit);
    }
  }
  else
  {
    Halfedge_const_handle he_of_pgn = ccb_of_pgn;
    this->_reset_faces();
    he_of_pgn->face()->set_visited(true);
    scanner.scan_ccb(ccb_of_pgn);
  }

  this->_reset_faces();
  return true;
}

template <class Traits_, class TopTraits_, class ValidationPolicy>
  typename Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::Ccb_halfedge_const_circulator
  Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  get_boundary_of_polygon(Face_const_iterator f) const
{
  CGAL_assertion(!f->visited());
  f->set_visited(true);

  if (f->number_of_outer_ccbs() == 0) // (f->is_unbounded())
  {
    return Ccb_halfedge_const_circulator();
  }

  // We assume that a polygon has only one outer_ccb. This code does not handle
  // the case where there are more than 1 outer ccbs. If this is the case, we
  // need to devise a method to convert the outer ccbs to inner ccbs so we
  // will have only one outer ccb.
  if (f->number_of_outer_ccbs() > 1)
    CGAL_error_msg("Not implemented yet.");

	// Some compilers (VC 9) do not like that we directly access the ccb_circ. So we have
	// to pass through the iterator.
  Outer_ccb_const_iterator oci_temp = f->outer_ccbs_begin();
  Ccb_halfedge_const_circulator ccb_end = *oci_temp;
  Ccb_halfedge_const_circulator ccb_circ = ccb_end;
  do
  {
    //get the current halfedge on the face boundary
    Halfedge_const_iterator he  = ccb_circ;
    Face_const_iterator new_f = he->twin()->face();
    if (!new_f->visited())
    {
      if (is_hole_of_face(new_f, he) && !new_f->contained())
        return (he->twin());
      return (get_boundary_of_polygon(new_f));
    }
    ++ccb_circ;
  }
  while(ccb_circ != ccb_end);
  CGAL_error();
  return Ccb_halfedge_const_circulator();

}

template <class Traits_, class TopTraits_, class ValidationPolicy>
  bool Gps_on_surface_base_2<Traits_, TopTraits_, ValidationPolicy>::
  is_hole_of_face(Face_const_handle f, Halfedge_const_handle he) const
{
  Inner_ccb_const_iterator   holes_it;
  for (holes_it = f->inner_ccbs_begin();
       holes_it != f->inner_ccbs_end(); ++holes_it)
  {
    Ccb_halfedge_const_circulator ccb = *holes_it;
    Ccb_halfedge_const_circulator ccb_end = ccb;
    do
    {
      Halfedge_const_handle he_inside_hole = ccb;
      he_inside_hole = he_inside_hole->twin();
      if (he == he_inside_hole)
        return true;

      ++ccb;
    }
    while(ccb != ccb_end);
  }

  return false;
}

} // namespace CGAL

#endif // CGAL_GPS_UTILS_H
