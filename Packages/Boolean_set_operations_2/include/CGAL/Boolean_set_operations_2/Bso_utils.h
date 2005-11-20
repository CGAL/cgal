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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef BSO_UTILS
#define BSO_UTILS


#include <list>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/iterator.h>
#include <CGAL/function_objects.h>
#include <CGAL/circulator.h> 




template <class Arrangement>
class Curve_creator : 
  public Creator_1<typename Arrangement::Halfedge,
                   typename Arrangement::X_monotone_curve_2>
{
public:

  typedef typename Arrangement::Halfedge              Halfedge;
  typedef typename Arrangement::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Arrangement::Traits_2              Traits_2;   
  typedef Creator_1<Halfedge, X_monotone_curve_2>     Base;

protected:
  Traits_2*    m_traits;



public:
  Curve_creator(Traits_2* tr): Base(),
                               m_traits(tr)
  {}


  X_monotone_curve_2 operator()(Halfedge he) const
  {
    if(he.direction() != 
       m_traits->compare_endpoints_xy_2_object()(he.curve()))
       return (m_traits->construct_opposite_2_object()(he.curve()));
    return (he.curve());
  }
};

template < class Traits_ >
void General_polygon_set_2<Traits_>::
  construct_polygon(Ccb_halfedge_circulator ccb,
                    Polygon_2&      pgn,
                    Traits_*                tr)
{
  typedef Container_from_circulator<Ccb_halfedge_circulator> Ccb_container;
  typedef typename Ccb_container::iterator                   Ccb_iterator;
  Ccb_container container(ccb);

  typedef Join_input_iterator_1<Ccb_iterator, Curve_creator<Arrangement_2> > Curve_iterator;
  Curve_creator<Arrangement_2> cv_creator(tr);

  Curve_iterator  begin(container.begin(), cv_creator);
  Curve_iterator  end  (container.end()  , cv_creator);

  tr->construct_polygon_2_object()(begin, end, pgn);
}

template <class Arrangement, class Visitor>
class Arr_bfs_scanner
{
public:

  typedef typename Arrangement::Ccb_halfedge_circulator 
                                                     Ccb_halfedge_circulator;
  typedef typename Arrangement::Face_iterator        Face_iterator;
  typedef typename Arrangement::Halfedge_iterator    Halfedge_iterator;
  typedef typename Arrangement::Holes_iterator       Holes_iterator;


protected:

  Visitor*                              m_visitor;
  std::queue<Ccb_halfedge_circulator>   m_holes_q;


public:

  /*! Constructor */
  Arr_bfs_scanner(Visitor& v) : m_visitor(&v)
  {}
                             

  void scan(Arrangement& arr)
  {
    m_visitor->visit_ccb(Ccb_halfedge_circulator(), true);
    Face_iterator ubf = arr.unbounded_face();
    m_visitor->mark_face(ubf);
    for(Holes_iterator hit = ubf->holes_begin(); hit != ubf->holes_end(); ++hit)
    {
      m_visitor->found_hole(ubf, hit);
      m_holes_q.push(*hit);
    }
    // finished the "infinite outer_ccb" of the unbounded face
    m_visitor->finish_ccb(Ccb_halfedge_circulator(), true);

    while(!m_holes_q.empty())
    {
      Ccb_halfedge_circulator ccb =  m_holes_q.front();
      scan(ccb);
      m_holes_q.pop();
    }
   
    m_visitor->reset_faces(arr.faces_begin(), arr.faces_end());
  }


  private:

  void scan(Ccb_halfedge_circulator ccb)
  {
    m_visitor->visit_ccb(ccb);
    Halfedge_iterator  hei = ccb;
    Face_iterator s_face = hei->face();
    Ccb_halfedge_circulator ccb_end = ccb;
    Ccb_halfedge_circulator ccb_circ = ccb_end;
    do
    { 
      //get the current halfedge on the face boundary
      Halfedge_iterator he  = ccb_circ;
      Face_iterator fi = he->twin()->face();
      if(!m_visitor->is_marked(fi))
      {
        all_incident_faces(s_face, fi);
      }
      ++ccb_circ;
    }
    while(ccb_circ != ccb_end);
    m_visitor->finish_ccb(ccb_circ);
  }


  void all_incident_faces(Face_iterator s_face, Face_iterator f)
  {
    CGAL_assertion(!m_visitor->is_marked(f));
     
    m_visitor->mark_face(f);
    m_visitor->flip_face(s_face, f);

    for(Holes_iterator hit = f->holes_begin(); hit != f->holes_end(); ++hit)
    {
      m_visitor->found_hole(f, hit);
      m_holes_q.push(*hit);
    }

    Ccb_halfedge_circulator ccb_end = f->outer_ccb();
    Ccb_halfedge_circulator ccb_circ = ccb_end;
    do
    { 
      //get the current halfedge on the face boundary
      Halfedge_iterator he  = ccb_circ;
      Face_iterator new_f = he->twin()->face();
      if(!m_visitor->is_marked(new_f))
      {
        all_incident_faces(f, new_f);
      }
      ++ccb_circ;
    }
    while(ccb_circ != ccb_end);
  }
};


template <class Arrangement>
class Base_visitor
{
public:

  typedef typename Arrangement::Face_iterator       Face_iterator;
  typedef typename Arrangement::Ccb_halfedge_circulator 
                                                    Ccb_halfedge_circulator;
  typedef typename Arrangement::Holes_iterator      Holes_iterator;
  
  Base_visitor(){}

  void flip_face(Face_iterator s_face, Face_iterator f)
  {}

  void visit_ccb(Ccb_halfedge_circulator ccb, bool is_infinite_ccb = false) 
  {}

  void finish_ccb(Ccb_halfedge_circulator ccb, bool is_infinite_ccb = false)
  {}

  void mark_face(Face_iterator f)
  {
    f->set_visited(true);
  }

  bool is_marked(Face_iterator f)
  {
    return (f->visited());
  }

  void found_hole(Face_iterator f, Holes_iterator hit)
  {}

  void reset_faces(Face_iterator begin, Face_iterator end)
  {
    for(; begin != end; ++begin)
    {
      begin->set_visited(false);
    }
  }
};


template <class Arrangement>
class Init_faces_visitor : public Base_visitor<Arrangement>
{
  typedef typename Arrangement::Face_iterator    Face_iterator;
public:

  void flip_face(Face_iterator s_face, Face_iterator f)
  {
    f->set_contained(!s_face->contained());
  }
};


template <class Arrangement, class OutputIterator>
class Construct_polygons_visitor : public Base_visitor<Arrangement>
{
public:

  typedef typename Arrangement::Traits_2          Traits_2;
  typedef typename Arrangement::Ccb_halfedge_circulator
                                                  Ccb_halfedge_circulator;
  typedef typename Arrangement::Halfedge_handle   Halfedge_handle;
  typedef typename Arrangement::Face_iterator     Face_iterator;
  typedef typename Arrangement::Holes_iterator    Holes_iterator;
  typedef typename Traits_2::Polygon_with_holes_2 
                                                  Polygon_with_holes_2;
  typedef typename Traits_2::Polygon_2    Polygon_2;
  typedef General_polygon_set_2<Traits_2>         Gps;

  Construct_polygons_visitor(Traits_2 & tr, OutputIterator out) : 
    m_boundary(),
    m_traits(&tr),
    m_out(out)
  {}



protected:
  Polygon_2               m_boundary;
  std::list<Polygon_2>    m_holes;
  Traits_2*                       m_traits;
  OutputIterator                  m_out;

public:

  void visit_ccb(Ccb_halfedge_circulator ccb, bool is_inf_ccb = false) 
  {
    if(is_inf_ccb)
      return;

    Halfedge_handle he = ccb;

    if(he->face()->contained())
    {
      // its a hole inside a polygon 
      return;
    }

    General_polygon_set_2<Traits_2>::construct_polygon(ccb, m_boundary, m_traits);
  }

  void found_hole(Face_iterator f, Holes_iterator hit)
  {
    if(!f->contained())
    {
      // its a polygon inside a non-contained area 
      return;
    }
    
    m_holes.push_back(Polygon_2());
    General_polygon_set_2<Traits_2>::construct_polygon(*hit, m_holes.back(), m_traits);
  }

  void finish_ccb(Ccb_halfedge_circulator ccb, bool is_inf_ccb = false)
  {
    if(is_inf_ccb)
      return;
    Halfedge_handle he = ccb;
    if(he->face()->contained())
    {
      //its a hole inside a polygon
      return;
    }

    *m_out = Polygon_with_holes_2(m_boundary, m_holes.begin(), m_holes.end());
    ++m_out;
    m_boundary = Polygon_2();  // clear the old polygon boundary
  }

  OutputIterator output_iterator()
  {
    return m_out;
  }
};


template <class Arrangement>
class Count_polygons_visitor : public Base_visitor<Arrangement>
{
public:

  typedef typename Arrangement::Traits_2          Traits_2;
  typedef typename Arrangement::Ccb_halfedge_circulator
                                                  Ccb_halfedge_circulator;
  typedef typename Arrangement::Halfedge_handle   Halfedge_handle;
  typedef typename Arrangement::Face_iterator     Face_iterator;
  typedef typename Arrangement::Holes_iterator    Holes_iterator;
  typedef typename Traits_2::Polygon_with_holes_2 
                                                  Polygon_with_holes_2;
  typedef typename Traits_2::Polygon_2    Polygon_2;

  Count_polygons_visitor(): m_num_of_p(0)
  {}

protected:
  unsigned int m_num_of_p;

public:

  void finish_ccb(Ccb_halfedge_circulator ccb, bool is_inf_ccb = false)
  {
    if(is_inf_ccb)
      return;
    Halfedge_handle he = ccb;
    if(he->face()->contained())
    {
      //its a hole inside a polygon
      return;
    }

    ++m_num_of_p;
  }

  unsigned int number_of_pgns() const
  {
    return m_num_of_p;
  }
};
  
template < class Traits_ >
void General_polygon_set_2<Traits_>::pgn2arr(const Polygon_2& pgn,
                                              Arrangement_2& arr)
{
  Compare_endpoints_xy_2  cmp_ends = 
    arr.get_traits()->compare_endpoints_xy_2_object();

  std::pair<Curve_const_iterator,
            Curve_const_iterator> itr_pair =
    arr.get_traits()->construct_curves_2_object()(pgn);

  if(itr_pair.first == itr_pair.second)
    return;

  Curve_const_iterator curr = itr_pair.first;
  Curve_const_iterator end  = itr_pair.second;

  Halfedge_handle first_he = 
    arr.insert_in_face_interior(*curr, arr.unbounded_face());
  //first_he is directed from left to right (see insert_in_face_interior)
  
  Halfedge_handle curr_he;
  if(cmp_ends(*curr) == SMALLER)
  {
    // curr curve and first_he have the same direction
    curr_he = first_he;
    first_he = first_he->twin();
  }
  else
  {
    // curr curve and first_he have opposite directions
    CGAL_assertion(cmp_ends(*curr) == LARGER);
    curr_he = first_he->twin();
  }

  Curve_const_iterator temp = curr;
  ++temp;
  if(temp == end) // a polygon with circular arcs may have only
                  // two edges (full circle for example)
  {
    Halfedge_handle he = 
      arr.insert_at_vertices(*temp, curr_he, first_he);
    if(he->face() == arr.unbounded_face())
      he->twin()->face()->set_contained(true);
    else
      he->face()->set_contained(true);
    return;
  }

  //The polygon has 3 or more edges
  Curve_const_iterator last = end;
  --last;
  for(++curr ; curr != last; ++curr)
  {
    const X_monotone_curve_2& curr_cv = *curr;
    if(cmp_ends(curr_cv) == SMALLER)
      curr_he = arr.insert_from_left_vertex(curr_cv, curr_he);
    else
    {
      CGAL_assertion(cmp_ends(curr_cv) == LARGER);
      curr_he = arr.insert_from_right_vertex(curr_cv, curr_he);
    }
  }

  const X_monotone_curve_2& last_cv = *last;
  Halfedge_handle last_he =
    arr.insert_at_vertices(last_cv, curr_he, first_he); 

  if(last_he->face() == arr.unbounded_face())
    last_he->twin()->face()->set_contained(true);
  else
    last_he->face()->set_contained(true);
}



//insert a range of simple polygons to the arrangement
template < class Traits_ >
  template< class PolygonIter >
void General_polygon_set_2<Traits_>::pgns2arr(PolygonIter p_begin,
                                              PolygonIter p_end,
                                              Arrangement_2&   arr)
{  
  typedef std::list<X_monotone_curve_2>                XCurveList;

  typedef Init_faces_visitor<Arrangement_2>              My_visitor;
  typedef Arr_bfs_scanner<Arrangement_2, My_visitor>     Arr_bfs_scanner;


  XCurveList xcurve_list;

  for(PolygonIter pitr = p_begin; pitr != p_end; ++pitr)
  {
    std::pair<Curve_const_iterator,
              Curve_const_iterator> itr_pair = 
      arr.get_traits()->construct_curves_2_object()(*pitr);
    std::copy(itr_pair.first, itr_pair.second, std::back_inserter(xcurve_list));
  }
  
  insert_non_intersecting(arr, xcurve_list.begin(), xcurve_list.end());

  My_visitor v;
  Arr_bfs_scanner scanner(v);
  scanner.scan(arr);
}



 //insert non-sipmle poloygons with holes (non incident edges may have
// common vertex,  but they dont intersect at their interior
template < class Traits_ >
void General_polygon_set_2<Traits_>::pgn_with_holes2arr (const Polygon_2& pgn,
                                                        Arrangement_2& arr)
{
  typedef std::list<X_monotone_curve_2>                XCurveList;

  typedef Init_faces_visitor<Arrangement_2>              My_visitor;
  typedef Arr_bfs_scanner<Arrangement_2, My_visitor>     Arr_bfs_scanner;

  XCurveList xcurve_list;
  
  bool is_bounded = !pgn.is_unbounded();
  if(is_bounded)
  {
    const Polygon_2& pgn_boundary = pgn.outer_boundary();
    arr.get_traits()->construct_curves_2_object()
    (pgn, std::back_inserter(xcurve_list));
  }

  for(GP_Holes_const_iterator hit = pgn.holes_begin();
      hit != pgn.holes_end();
      ++hit)
  {
    const Polygon_2& pgn_hole = *hit;
    arr.get_traits()->construct_curves_2_object()
    (pgn_hole, std::back_inserter(xcurve_list));
  }
  
  insert_non_intersecting(arr, xcurve_list.begin(), xcurve_list.end());

  if(!is_bounded)
    arr.unbounded_face()->set_contained(true);

  My_visitor v;
  Arr_bfs_scanner scanner(v);
  scanner.scan(arr);
}


 //insert a range of  non-sipmle poloygons with holes (as decribed above)
template < class Traits_ >
  template< class InputIterator >
void General_polygon_set_2<Traits_>::pgns_with_holes2arr (InputIterator begin,
                                                         InputIterator end,
                                                         Arrangement_2& arr)
{
  typedef std::list<X_monotone_curve_2>                XCurveList;

  typedef Init_faces_visitor<Arrangement_2>              My_visitor;
  typedef Arr_bfs_scanner<Arrangement_2, My_visitor>     Arr_bfs_scanner;

  XCurveList xcurve_list;

  bool is_bounded = true;
  InputIterator itr; 
  for(itr = begin; itr != end; ++itr)
  {
    if(!itr->is_unbounded())
    {
      const Polygon_2& pgn = itr->outer_boundary();
      arr.get_traits()->construct_curves_2_object()
        (pgn, std::back_inserter(xcurve_list));
    }
    else
      is_bounded = false;


    for(GP_Holes_const_iterator hit = pgn.holes_begin;
        hit != pgn.holes_end();
        ++hit)
    {
      const Polygon_2& pgn_hole = *hit;
     arr.get_traits()->construct_curves_2_object()
      (pgn_hole, std::back_inserter(xcurve_list));
    }
  }

  insert_non_intersecting(arr, xcurve_list.begin(), xcurve_list.end());

  if(!is_bounded)
    arr.unbounded_face()->set_contained(true);

  My_visitor v;
  Arr_bfs_scanner scanner(v);
  scanner.scan(arr);
}



template < class Traits_ >
  template< class OutputIterator>
  OutputIterator  General_polygon_set_2<Traits_>::
    polygons(OutputIterator out1)
{
  typedef Construct_polygons_visitor<Arrangement_2,
                                     OutputIterator>     My_visitor;
  typedef Arr_bfs_scanner<Arrangement_2, My_visitor>     Arr_bfs_scanner;

  My_visitor v(*(this->m_traits), out1);
  Arr_bfs_scanner scanner(v);
  scanner.scan(*m_arr);

  return (v.output_iterator());
}

template < class Traits_ >
typename  General_polygon_set_2<Traits_>::Size 
General_polygon_set_2<Traits_>::number_of_polygons_with_holes() const
{
  typedef Count_polygons_visitor<Arrangement_2>          My_visitor;
  typedef Arr_bfs_scanner<Arrangement_2, My_visitor>     Arr_bfs_scanner;

  My_visitor v;
  Arr_bfs_scanner scanner(v);
  scanner.scan(*m_arr);
  
  return (v.number_of_pgns());
}


#endif
