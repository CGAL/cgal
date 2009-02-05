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
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_GPS_UTILS_H
#define CGAL_GPS_UTILS_H

#include <CGAL/Unique_hash_map.h>
#include <CGAL/iterator.h>
#include <CGAL/function_objects.h>
#include <CGAL/circulator.h> 
#include <CGAL/Boolean_set_operations_2/Gps_bfs_scanner.h>
#include <CGAL/Arr_accessor.h>

#include <queue>

template <class Traits_, class Dcel_>
void General_polygon_set_2<Traits_, Dcel_>::
construct_polygon(Ccb_halfedge_const_circulator ccb, Polygon_2 & pgn,
                  Traits_ * tr)
{
  typedef CGAL::Ccb_curve_iterator<Arrangement_2>    Ccb_curve_iterator;
  Ccb_curve_iterator begin(ccb, false);
  Ccb_curve_iterator end(ccb, true);

  tr->construct_polygon_2_object()(begin, end, pgn);
}


template <class Arrangement, class OutputIterator>
class Arr_bfs_scanner
{
public:
  typedef typename Arrangement::Traits_2                Gps_traits;
  typedef typename Arrangement::Dcel                    Gps_dcel;
  typedef typename Gps_traits::Polygon_2                Polygon_2;
  typedef typename Gps_traits::Polygon_with_holes_2     Polygon_with_holes_2;
  typedef typename Arrangement::Ccb_halfedge_const_circulator 
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement::Face_const_iterator     Face_const_iterator;
  typedef typename Arrangement::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Arrangement::Inner_ccb_const_iterator     Hole_const_iterator;


protected:

  Gps_traits*                            m_traits;
  std::queue<Face_const_iterator>        m_holes_q;
  std::list<Polygon_2>                   m_pgn_holes;
  OutputIterator                         m_oi;

public:

  /*! Constructor */
  Arr_bfs_scanner(Gps_traits* tr, OutputIterator oi) : m_traits(tr), m_oi(oi)
  {}
                             

  void scan(Arrangement& arr)
  {
    Face_const_iterator   ubf = arr.unbounded_face();  
    Hole_const_iterator  holes_it;
    Face_const_iterator   fit;

    if (!ubf->contained())
    {
      ubf->set_visited(true);
      for (holes_it = ubf->holes_begin();
            holes_it != ubf->holes_end(); ++holes_it)
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
      for (holes_it = top_f->holes_begin();
           holes_it != top_f->holes_end(); ++holes_it)
      {
        scan_ccb(*holes_it);
      }

      //scan_uncontained_face(top_f->outer_ccb());
    }

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
    General_polygon_set_2<Gps_traits, Gps_dcel>::
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
    Polygon_with_holes_2 pgn = m_traits->construct_polygon_with_holes_2_object()(pgn_boundary,
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
    CGAL_assertion(ubf->is_unbounded() && ubf->contained());
    // ubf is contained -> unbounded polygon !!
    all_incident_faces(ubf);
    Polygon_2 boundary;
    Polygon_with_holes_2 pgn = m_traits->construct_polygon_with_holes_2_object()(boundary,
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
    if (!f->is_unbounded())
    {
      if (!f->contained())
      {
        m_pgn_holes.push_back(Polygon_2());
        General_polygon_set_2<Gps_traits, Gps_dcel>::
          construct_polygon(f->outer_ccb(), m_pgn_holes.back(), m_traits);
        m_holes_q.push(f);
      }

    
      Ccb_halfedge_const_circulator ccb_end = f->outer_ccb();
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

    if (f->contained())
    {
      Hole_const_iterator hit;
      for(hit = f->holes_begin(); hit != f->holes_end(); ++hit)
      {
        Ccb_halfedge_const_circulator ccb_of_hole = *hit;
        Halfedge_const_iterator he = ccb_of_hole;
        if (is_single_face(ccb_of_hole))
        {
          CGAL_assertion(!he->twin()->face()->contained());
         
          m_pgn_holes.push_back(Polygon_2());
          General_polygon_set_2<Gps_traits, Gps_dcel>::
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

  void flip_face(Face_iterator f1, Face_iterator f2, Halfedge_iterator /*he*/)
  {
    f2->set_contained(!f1->contained());
  }
};
  
template <class Traits_, class Dcel_>
void General_polygon_set_2<Traits_, Dcel_>::
_insert(const Polygon_2& pgn, Arrangement_2 & arr)
{
  typedef Arr_accessor<Arrangement_2>                  Arr_accessor;

  Arr_accessor  accessor(arr);
  Compare_endpoints_xy_2  cmp_ends = m_traits->compare_endpoints_xy_2_object();

  std::pair<Curve_const_iterator, Curve_const_iterator> itr_pair =
    m_traits->construct_curves_2_object()(pgn);

  if (itr_pair.first == itr_pair.second)
    return;

  Curve_const_iterator curr = itr_pair.first;
  Curve_const_iterator end  = itr_pair.second;

  Face_iterator f;
  if (arr.is_empty())
  {
    f = arr.unbounded_face();
  }
  else
  {
    Walk_pl pl(arr);

    Object obj = pl.locate(m_traits->construct_min_vertex_2_object()(*curr));

    Face_const_iterator const_f;
    // pgn must be completely disjoint from the arrangement
    CGAL_assertion(CGAL::assign(const_f, obj) && !const_f->contained());
    CGAL::assign(const_f, obj);
    f = arr.non_const_handle(const_f);
  }
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
    bool new_face_created;
    Halfedge_handle he = accessor.insert_at_vertices_ex (*temp,
                                                         curr_he,
                                                         first_he,
                                                         cmp_ends(*temp),
                                                         new_face_created);
    CGAL_assertion(new_face_created); 
    CGAL_assertion((he->face() != he->twin()->face()) && 
                   (he->face() != arr.unbounded_face()));
    
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
  bool new_face_created;
  Halfedge_handle last_he = 
    accessor.insert_at_vertices_ex (last_cv,
                                    curr_he,
                                    first_he,
                                    cmp_ends(last_cv),
                                    new_face_created);
  CGAL_assertion(new_face_created); 
  CGAL_assertion((last_he->face() != last_he->twin()->face()) && 
                 (last_he->face() != arr.unbounded_face()));
  
  last_he->face()->set_contained(true);
}


template <class Traits_, class Dcel_>
template<class PolygonIter >
void General_polygon_set_2<Traits_, Dcel_>::
insert(PolygonIter p_begin, PolygonIter p_end)
{
  typename std::iterator_traits<PolygonIter>::value_type pgn;
   //check validity of all polygons    
   for( ; p_begin != p_end; ++p_begin)
  {
     CGAL_precondition(is_valid_unkown_polygon(*p_begin, *m_traits));
  }
  _insert(p_begin, p_end, pgn);
}

template <class Traits_, class Dcel_>
template<class PolygonIter, class PolygonWithHolesIter>
void General_polygon_set_2<Traits_, Dcel_>::
insert(PolygonIter p_begin, PolygonIter p_end,
       PolygonWithHolesIter pwh_begin, PolygonWithHolesIter pwh_end)
{
  typedef std::list<X_monotone_curve_2>                  XCurveList;
  typedef Init_faces_visitor<Arrangement_2>              My_visitor;
  typedef Gps_bfs_scanner<Arrangement_2, My_visitor>     Arr_bfs_scanner;

  XCurveList xcurve_list;
  
  for( ; p_begin != p_end; ++p_begin)
  {
      CGAL_precondition(is_valid_polygon(*p_begin, *m_traits));
    _construct_curves(*p_begin, std::back_inserter(xcurve_list));
  }

  bool is_unbounded = false;
  for( ; pwh_begin != pwh_end; ++pwh_begin)
  {
    CGAL_precondition(is_valid_polygon_with_holes(*pwh_begin, *m_traits));
     is_unbounded = (is_unbounded || m_traits->construct_is_unbounded_object()(*pwh_begin));
   // is_unbounded = (is_unbounded || pwh_begin->is_unbounded());
    _construct_curves(*pwh_begin, std::back_inserter(xcurve_list));
  }
  insert_non_intersecting_curves(*m_arr, xcurve_list.begin(), xcurve_list.end());

  if (is_unbounded)
    m_arr->unbounded_face()->set_contained(true);

  My_visitor v;
  Arr_bfs_scanner scanner(v);
  scanner.scan(*m_arr);
  _reset_faces(m_arr);
}

//insert a range of simple polygons to the arrangement
template <class Traits_, class Dcel_>
template<class PolygonIter>
void General_polygon_set_2<Traits_, Dcel_>::
_insert(PolygonIter p_begin, PolygonIter p_end, Polygon_2 & /*pgn*/)
{  
  for(PolygonIter pitr = p_begin; pitr != p_end; ++pitr)
  {
        this->_insert(*pitr, *m_arr);
  }
}

template <class Traits_, class Dcel_>
template<class PolygonIter>
void General_polygon_set_2<Traits_, Dcel_>::
_insert(PolygonIter p_begin, PolygonIter p_end, Polygon_with_holes_2 & /*pgn*/)
{  
  typedef std::list<X_monotone_curve_2>                  XCurveList;
  typedef Init_faces_visitor<Arrangement_2>              My_visitor;
  typedef Gps_bfs_scanner<Arrangement_2, My_visitor>     Arr_bfs_scanner;

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
    m_arr->unbounded_face()->set_contained(true);

  My_visitor v;
  Arr_bfs_scanner scanner(v);
  scanner.scan(*m_arr);
  _reset_faces(m_arr);
}


 //insert non-sipmle poloygons with holes (non incident edges may have
// common vertex,  but they dont intersect at their interior
template <class Traits_, class Dcel_>
void General_polygon_set_2<Traits_, Dcel_>::
_insert(const Polygon_with_holes_2 & pgn, Arrangement_2 & arr)
{
  //not needed gps.insert(PWH) has the precondition
 // CGAL_precondition(is_valid_polygon_with_holes(pgn, *m_traits));
  typedef std::list<X_monotone_curve_2>                  XCurveList;
  typedef Init_faces_visitor<Arrangement_2>              My_visitor;
  typedef Gps_bfs_scanner<Arrangement_2, My_visitor>     Arr_bfs_scanner;

  XCurveList xcurve_list;
  _construct_curves(pgn, std::back_inserter(xcurve_list));
  insert_non_intersecting_curves(arr, xcurve_list.begin(), xcurve_list.end());

  //if (pgn.is_unbounded())
  if (m_traits->construct_is_unbounded_object()(pgn))	  
    arr.unbounded_face()->set_contained(true);

  My_visitor v;
  Arr_bfs_scanner scanner(v);
  scanner.scan(arr);
  _reset_faces(&arr);
}

template <class Traits_, class Dcel_>
template <class OutputIterator>
void 
General_polygon_set_2<Traits_, Dcel_>::
_construct_curves(const Polygon_2 & pgn, OutputIterator oi)
{
    std::pair<Curve_const_iterator,
              Curve_const_iterator> itr_pair =
              m_traits->construct_curves_2_object()(pgn);
    std::copy (itr_pair.first, itr_pair.second, oi);
}

template <class Traits_, class Dcel_>
template <class OutputIterator>
void General_polygon_set_2<Traits_, Dcel_>::
_construct_curves(const Polygon_with_holes_2 & pgn, OutputIterator oi)
{
  //if (!pgn.is_unbounded())
  if (!m_traits->construct_is_unbounded_object()(pgn))
  {
    const Polygon_2& pgn_boundary = m_traits->construct_outer_boundary_object ()(pgn);
    std::pair<Curve_const_iterator,
              Curve_const_iterator> itr_pair = 
              m_traits->construct_curves_2_object()(pgn_boundary);
    std::copy (itr_pair.first, itr_pair.second, oi);
  }
  std::pair<GP_Holes_const_iterator, GP_Holes_const_iterator> hpair = m_traits->construct_holes_object()(pgn);
  GP_Holes_const_iterator hit;
  for (hit = hpair.first; hit != hpair.second; ++hit)
  {
    const Polygon_2& pgn_hole = *hit;
    std::pair<Curve_const_iterator,
              Curve_const_iterator> itr_pair =
              m_traits->construct_curves_2_object()(pgn_hole);
    std::copy (itr_pair.first, itr_pair.second, oi);
  }
}

template <class Traits_, class Dcel_>
template <class OutputIterator>
OutputIterator
General_polygon_set_2<Traits_, Dcel_>::
polygons_with_holes(OutputIterator out) const
{
  typedef Arr_bfs_scanner<Arrangement_2, OutputIterator>     Arr_bfs_scanner;
  Arr_bfs_scanner scanner(this->m_traits, out);
  scanner.scan(*(this->m_arr));
  return (scanner.output_iterator());
}


template <class Traits_, class Dcel_>
typename General_polygon_set_2<Traits_, Dcel_>::Size 
General_polygon_set_2<Traits_, Dcel_>::
number_of_polygons_with_holes() const
{
 
  typedef Arr_bfs_scanner<Arrangement_2, Counting_output_iterator>
    Arr_bfs_scanner;
  //counting_output_operator CTOR reqires a parameter  
  std::size_t *cc = new size_t();  
  Arr_bfs_scanner scanner(this->m_traits, Counting_output_iterator(cc));
  scanner.scan(*(this->m_arr));
  return (scanner.output_iterator().current_counter());
}

template <class Traits_, class Dcel_>
bool General_polygon_set_2<Traits_, Dcel_>::
locate(const Point_2& q, Polygon_with_holes_2& pgn) const
{
  Walk_pl pl(*m_arr);

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
  typedef Arr_bfs_scanner<Arrangement_2, OutputItr>     Arr_bfs_scanner;

  OutputItr oi (pgn);
  Arr_bfs_scanner scanner(this->m_traits, oi);
  
  
  Ccb_halfedge_const_circulator ccb_of_pgn = get_boundary_of_polygon(f);
  this->_reset_faces();
  if (ccb_of_pgn == Ccb_halfedge_const_circulator()) // the polygon has no boundary
  {
    // f is unbounded 
    scanner.scan_contained_ubf(m_arr->unbounded_face());
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

template <class Traits_, class Dcel_>
typename General_polygon_set_2<Traits_, Dcel_>::Ccb_halfedge_const_circulator
General_polygon_set_2<Traits_, Dcel_>::
get_boundary_of_polygon(Face_const_iterator f) const
{
  CGAL_assertion(!f->visited());
  f->set_visited(true);
  
  if (f->is_unbounded())
  {
    return Ccb_halfedge_const_circulator();
  }
  Ccb_halfedge_const_circulator ccb_end = f->outer_ccb();
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

template <class Traits_, class Dcel_>
bool General_polygon_set_2<Traits_, Dcel_>::
is_hole_of_face(Face_const_handle f, Halfedge_const_handle he) const
{
  Hole_const_iterator   holes_it;
  for (holes_it = f->holes_begin(); holes_it != f->holes_end(); ++holes_it)
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

#endif
