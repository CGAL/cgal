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

#ifndef GENERAL_POLYGON_SET_2_H
#define GENERAL_POLYGON_SET_2_H

#include <CGAL/Boolean_set_operations_2/Bso_dcel.h>
#include <CGAL/Boolean_set_operations_2/Bso_do_intersect_functor.h>
#include <CGAL/Boolean_set_operations_2/Bso_intersection_functor.h>
#include <CGAL/Boolean_set_operations_2/Bso_join_functor.h>
#include <CGAL/Boolean_set_operations_2/Bso_difference_functor.h>
#include <CGAL/Boolean_set_operations_2/Bso_sym_diff_functor.h>
#include <CGAL/Boolean_set_operations_2/Gps_merge.h>


#include <CGAL/Arr_walk_along_line_point_location.h>

#include <CGAL/Arr_overlay.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Object.h>
#include <CGAL/enum.h>
#include <CGAL/iterator.h> 
#include <queue>

CGAL_BEGIN_NAMESPACE

// General_polygon_set_2
template < class Traits_ >
class General_polygon_set_2 
{
public:

  typedef Traits_                                         Traits_2;

  typedef typename Traits_2::Point_2                      Point_2;
  typedef typename Traits_2::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits_2::Polygon_2                    Polygon_2;
  typedef typename Traits_2::Polygon_with_holes_2         Polygon_with_holes_2;
  typedef typename Polygon_with_holes_2::Holes_const_iterator  
                                                 GP_Holes_const_iterator;
  typedef typename Traits_2::Equal_2                      Equal_2;
  typedef typename Traits_2::Compare_xy_2                 Compare_xy_2;
  typedef typename Traits_2::Construct_min_vertex_2       Construct_min_vertex_2;
  typedef typename Traits_2::Construct_max_vertex_2       Construct_max_vertex_2;
  typedef typename Traits_2::Construct_polygon_2          Construct_polygon_2;
  typedef typename Traits_2::Construct_curves_2           Construct_curves_2;
  typedef typename Traits_2::Curve_const_iterator         Curve_const_iterator;
  typedef typename Traits_2::Compare_endpoints_xy_2       Compare_endpoints_xy_2;
  typedef typename Traits_2::Compare_endpoints_xy_2       Compare_endpoints_xy_2;
  typedef typename Traits_2::Construct_opposite_2         Construct_opposite_2;

  typedef Bso_dcel<Traits_2>                                Bso_dcel;  
  typedef Arrangement_2<Traits_2, Bso_dcel>                 Arrangement_2;

  typedef typename Arrangement_2::Face_const_iterator        Face_const_iterator;
  typedef typename Arrangement_2::Halfedge_const_iterator    Halfedge_const_iterator;
  typedef typename Arrangement_2::Vertex_const_iterator      Vertex_const_iterator;
  typedef typename Arrangement_2::Edge_const_iterator        Edge_const_iterator;
  typedef typename Arrangement_2::Holes_const_iterator       Holes_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator 
                                                      Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Face_iterator              Face_iterator;
  typedef typename Arrangement_2::Halfedge_iterator          Halfedge_iterator;
  typedef typename Arrangement_2::Vertex_iterator            Vertex_iterator;
  typedef typename Arrangement_2::Edge_iterator              Edge_iterator;
  typedef typename Arrangement_2::Holes_iterator             Holes_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_circulator 
                                                      Ccb_halfedge_circulator;
  typedef typename Arrangement_2::Face_handle                Face_handle;
  typedef typename Arrangement_2::Halfedge_handle            Halfedge_handle;
  typedef typename Arrangement_2::Vertex_handle              Vertex_handle;

  typedef typename Arrangement_2::Face_const_handle          Face_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle      Halfedge_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle        Vertex_const_handle;

  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
                                            Halfedge_around_vertex_const_circulator;

  typedef Arr_walk_along_line_point_location<Arrangement_2>  Walk_pl;

  typedef unsigned int                                       Size;


protected:

  Traits_2*             m_traits;
  bool                  m_traits_owner;

  // the underlying arrangement
  Arrangement_2*        m_arr;


public:

  // default costructor
  General_polygon_set_2() : m_traits(new Traits_2()),
                            m_traits_owner(true),
                            m_arr(new Arrangement_2(m_traits))       
  {}


  // constructor with traits object
  General_polygon_set_2(Traits_2& tr) : m_traits(&tr),
                                        m_traits_owner(false),
                                        m_arr(new Arrangement_2(m_traits)) 
  {}


  General_polygon_set_2(const General_polygon_set_2& ps):  m_traits(new Traits_2(*(ps.m_traits))),
                                                           m_traits_owner(true),
                                                           m_arr(new Arrangement_2(*(ps.m_arr)))
  {}

  General_polygon_set_2& operator=(const General_polygon_set_2& ps)
  {
    if(this == &ps)
      return (*this);

    if(m_traits_owner)
      delete m_traits;
    delete m_arr;
    m_traits = new Traits_2(*(ps.m_traits));
    m_traits_owner = true;
    m_arr = new Arrangement_2(*(ps.m_arr));
    return (*this);
  }


  explicit General_polygon_set_2(const Polygon_2& pgn) : 
    m_traits(new Traits_2()),
    m_traits_owner(true),
    m_arr(new Arrangement_2(m_traits)) 
  {
    pgn2arr(pgn, *m_arr);
  }

  explicit General_polygon_set_2(const Polygon_with_holes_2& pgn_with_holes): 
    m_traits(new Traits_2()),
    m_traits_owner(true),
    m_arr(new Arrangement_2(m_traits))
  {
    pgn_with_holes2arr(pgn_with_holes, *m_arr);
  }


  // constructor of range of polygons that can be either simple polygons
  // or polygons with holes
  // precondition: the polygons are disjoint and simple
  template <class PolygonIterator>
    General_polygon_set_2(PolygonIterator pgn_begin,
                          PolygonIterator pgn_end): 
    m_traits(new Traits_2()),
    m_traits_owner(true),
    m_arr(new Arrangement_2(m_traits)) 
  {}


  // constructor of two ranges of : the first one for simple polygons,
  // the second one for polygons with holes
  // precondition: the first range is disjoint simple polygons 
  //               the second range is fisjoint polygons with holes
  template <class PolygonIterator, class PolygonWithHolesIterator>
  General_polygon_set_2(PolygonIterator pgn_begin,
                        PolygonIterator pgn_end,
                        PolygonWithHolesIterator  pgn_with_holes_begin,
                        PolygonWithHolesIterator  pgn_with_holes_end):
    m_traits(new Traits_2()),
    m_traits_owner(true),
    m_arr(new Arrangement_2(m_traits)) 
 {}


  //destructor
  virtual ~General_polygon_set_2()
  {
    delete m_arr;

    if(m_traits_owner)
      delete m_traits;
  }


  // insert a simple polygon
  void insert(const Polygon_2& pgn)
  {
    pgn2arr(pgn, *m_arr);
  }

  // insert a polygon with holes
  void insert(const Polygon_with_holes_2& pgn_with_holes);


  // insert a range of polygons that can be either simple polygons
  // or polygons with holes
  // precondition: the polygons are disjoint and simple
  template <class PolygonIterator>
  void insert(PolygonIterator pgn_begin, PolygonIterator pgn_end);


  // insert two ranges of : the first one for simple polygons,
  // the second one for polygons with holes
  // precondition: the first range is disjoint simple polygons 
  //               the second range is fisjoint polygons with holes
  template <class PolygonIterator, class PolygonWithHolesIterator>
  void insert(PolygonIterator pgn_begin,
              PolygonIterator pgn_end,
              PolygonWithHolesIterator  pgn_with_holes_begin,
              PolygonWithHolesIterator  pgn_with_holes_end);



   // test for intersection with a simple polygon
  bool do_intersect(const Polygon_2 &pgn) const
  {
    if(this->is_empty())
      return false;
    
    Arrangement_2 second_arr;
    pgn2arr(pgn, second_arr);
    if(second_arr.is_empty())
      return false;

    Arrangement_2 res_arr;
    Bso_do_intersect_functor<Traits_2>  func;
    overlay(*m_arr, second_arr, res_arr, func);
    return func.found_intersection();
  }

  // test for intersection with a polygon with holes
  bool do_intersect(const Polygon_with_holes_2& pgn_with_holes)
  {
    if(this->is_empty())
      return false;
    Arrangement_2 second_arr;
    pgn_with_holes2arr(pgn_with_holes, second_arr);
    if(second_arr.is_empty())
      return false;

    Arrangement_2 res_arr;
    Bso_do_intersect_functor<Traits_2>  func;
    overlay(*m_arr, second_arr, res_arr, func);
    return func.found_intersection();
  }

  //test for intersection with another General_polygon_set_2 object
  bool do_intersect(const General_polygon_set_2& bops)
  {
    if(this->is_empty() || bops.is_empty())
      return false;
    
    Arrangement_2 res_arr;

    Bso_do_intersect_functor<Traits_2>  func;
    overlay(*m_arr, *(bops.m_arr), res_arr, func);
    return func.found_intersection();
  }


  // intersection with a simple polygon
  void intersection(const Polygon_2& pgn)
  {
    if(this->is_empty())
    {
      return;
    }
    Arrangement_2 second_arr;
    pgn2arr(pgn, second_arr);
    if(second_arr.is_empty())
    {
      this->clear();
      return;
    }
    Arrangement_2 res_arr;

    Bso_intersection_functor<Traits_2> func(m_traits);
    overlay(*m_arr, second_arr, res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }

  // intersection with a polygon with holes
  void intersection(const Polygon_with_holes_2& pgn_with_holes)
  {
    if(this->is_empty())
    {
      return;
    }
    Arrangement_2 second_arr;
    Arrangement_2 res_arr;
    pgn_with_holes2arr(pgn_with_holes, second_arr);
    if(second_arr.is_empty())
    {
      if(! second_arr.unbounded_face()->contained())
        this ->clear();
     
      return;
    }


    Bso_intersection_functor<Traits_2>  func(m_traits);
    overlay(*m_arr, second_arr, res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }

  //intersection with another General_polygon_set_2 object
  void intersection(const General_polygon_set_2& bops)
  {
    Arrangement_2 res_arr;

    Bso_intersection_functor<Traits_2> func(m_traits);
    overlay(*m_arr, *(bops.m_arr), res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }


  // join with a simple polygon
  void join(const Polygon_2& pgn)
  {
    if(this->is_plane())
    {
      return;
    }

    if(this->is_empty())
    {
      pgn2arr(pgn, *m_arr);
      return;
    }

    Arrangement_2 second_arr;
    Arrangement_2 res_arr;
    pgn2arr(pgn, second_arr);

    Bso_join_functor<Traits_2>  func(m_traits);
    overlay(*m_arr, second_arr, res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }

  // join with a polygon with holes
  void join(const Polygon_with_holes_2& pgn_with_holes)
  {
     if(this->is_plane())
    {
      return;
    }

    if(this->is_empty())
    {
      pgn_with_holes2arr(pgn_with_holes, *m_arr);
      return;
    }

    Arrangement_2 second_arr;
    Arrangement_2 res_arr;
    pgn_with_holes2arr(pgn_with_holes, second_arr);

    Bso_join_functor<Traits_2>  func(m_traits);
    overlay(*m_arr, second_arr, res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }

  //join with another General_polygon_set_2 object
  void join(const General_polygon_set_2& bops)
  {
    Arrangement_2 res_arr;

    Bso_join_functor<Traits_2>  func(m_traits);
    overlay(*m_arr, *(bops.m_arr), res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }


   // difference with a simple polygon
  void difference (const Polygon_2& pgn)
  {
    if(this->is_empty())
    {
      return;
    }

    Arrangement_2 second_arr;
    Arrangement_2 res_arr;
    pgn2arr(pgn, second_arr);
    if(second_arr.is_empty())
    {
      return;
    }

    Bso_difference_functor<Traits_2>  func(m_traits);
    overlay(*m_arr, second_arr, res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }

  // difference with a polygon with holes
  void difference (const Polygon_with_holes_2& pgn_with_holes)
  {
    if(this->is_empty())
    {
      return;
    }

    Arrangement_2 second_arr;
    Arrangement_2 res_arr;
    pgn_with_holes2arr(pgn_with_holes, second_arr);
    if(second_arr.is_empty())
    {
      return;
    }

    Bso_difference_functor<Traits_2>  func(m_traits);
    overlay(*m_arr, second_arr, res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }

  //difference with another General_polygon_set_2 object
  void difference (const General_polygon_set_2& bops)
  {
   Arrangement_2 res_arr;

    Bso_difference_functor<Traits_2>  func(m_traits);
    overlay(*m_arr,  *(bops.m_arr), res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }


  // symmetric_difference with a simple polygon
  void symmetric_difference(const Polygon_2& pgn)
  {
    Arrangement_2 second_arr;
    Arrangement_2* res_arr = new Arrangement_2(m_traits);
    pgn2arr(pgn, second_arr);
    Bso_sym_diff_functor<Traits_2>  func;
    overlay(*m_arr, second_arr, *res_arr, func);
    std::vector<Halfedge_handle>& he_vec = func.get_spare_edges();
    for(unsigned int i=0; i<he_vec.size(); ++i)
    {
      res_arr->remove_edge(he_vec[i]);
    }

    delete m_arr; // delete the previous arrangement
    m_arr = res_arr;
  }

  // symmetric_difference with a polygon with holes
  void symmetric_difference(const Polygon_with_holes_2& pgn_with_holes)
  {
    Arrangement_2 second_arr;
    Arrangement_2* res_arr = new Arrangement_2(m_traits);
    pgn_with_holes2arr(pgn_with_holes, second_arr);

    Bso_sym_diff_functor<Traits_2>  func;
   
    overlay(*m_arr, second_arr, *res_arr, func);
    std::vector<Halfedge_handle>& he_vec = func.get_spare_edges();
    for(unsigned int i=0; i<he_vec.size(); ++i)
    {
      res_arr->remove_edge(he_vec[i]);
    }
    
    delete m_arr; // delete the previous arrangement
    m_arr = res_arr;
  }

  //symmetric_difference with another General_polygon_set_2 object
  void symmetric_difference(const General_polygon_set_2& bops)
  {
    Arrangement_2* res_arr = new Arrangement_2(m_traits);
    Bso_sym_diff_functor<Traits_2>  func;
    overlay(*m_arr, *(bops.m_arr), *res_arr, func);
    std::vector<Halfedge_handle>& he_vec = func.get_spare_edges();
    for(unsigned int i=0; i<he_vec.size(); ++i)
    {
      res_arr->remove_edge(he_vec[i]);
    }

    delete m_arr; // delete the previous arrangement
    m_arr = res_arr;
  }

  void complement()
  {
    for(Face_iterator fit = m_arr->faces_begin();
        fit != m_arr->faces_end();
        ++fit)
    {
      fit->set_contained(!fit->contained());
    }

    Construct_opposite_2 ctr_opp = m_traits->construct_opposite_2_object();
    for(Edge_iterator eit = m_arr->edges_begin();
        eit != m_arr->edges_end();
        ++eit)
    {
      Halfedge_handle he = eit;
      const X_monotone_curve_2& cv = he->curve();
      m_arr->modify_edge(he, ctr_opp(cv));
    }
  }

  Size number_of_polygons_with_holes() const;

  Traits_2& traits()
  {
    return *m_traits;
  }

  bool is_empty() const
  {
    return (m_arr->is_empty() && ! m_arr->unbounded_face()->contained());
  }

  bool is_plane() const
  {
    return (m_arr->is_empty() &&  m_arr->unbounded_face()->contained());
  }

  void clear()
  {
    m_arr->clear();
  }

  
  bool is_inside(const Point_2& q)
  {
    Walk_pl pl(*m_arr);

    Object obj = pl.locate(q);
    Face_const_iterator f;
    if(CGAL::assign(f, obj))
      return (f->contained());
    return true;
  }


  // returns the location of the query point
  bool locate(const Point_2& q, Polygon_with_holes_2& pgn) const;

  //advanced function: get const reference to the arrangement
  const Arrangement_2& arrangement() const
  {
    return *m_arr;
  }

  bool is_valid() const
  {
    if(! is_valid(m_arr))
      return false;

    for(Edge_const_iterator eci = m_arr->edges_begin();
        eci != m_arr->edges_end();
        ++eci)
    {
      Halfedge_const_handle he = eci;
      if(he->face() == he->twin()->face())
      {
        return false;
      }
      if(he->face()->contained() == he->twin()->face()->contained())
      {
        return false;
      }
    }
    return true;
  }


  //static void pgn2arr (const Polygon_2& pgn, Arrangement_2& arr);

  template< class PolygonIter >
  void pgns2arr(PolygonIter p_begin, PolygonIter p_end, Arrangement_2& arr);
  
  //void pgn_with_holes2arr (const Polygon_with_holes_2& pgn, Arrangement_2& arr);

  
  template< class InputIterator >
  void pgns_with_holes2arr (InputIterator begin, 
                            InputIterator end,
                            Arrangement_2& arr);



  // get the simple polygons, takes O(n)
  template <class OutputIterator>
  OutputIterator polygons_with_holes(OutputIterator out) const;

 // join a range of polygons
 template <class InputIterator>
 void join(InputIterator begin, InputIterator end, unsigned int k = 5)
{
  typename std::iterator_traits<InputIterator>::value_type pgn;
  this->join(begin, end, pgn, k);
  this->remove_redundant_edges();
}

// join range of simple polygons
template <class InputIterator>
inline void join(InputIterator begin,
                 InputIterator end,
                 Polygon_2& pgn,
                 unsigned int k=5)
{
  std::vector<Arrangement_2*> arr_vec (std::distance(begin, end) + 1);
 
  arr_vec[0] = this->m_arr;
  unsigned int i = 1;
  for(InputIterator itr = begin; itr!=end; ++itr, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn2arr(*itr, *arr_vec[i]);
  }

  Join_merge<Arrangement_2> join_merge;
  divide_and_conquer(0, arr_vec.size()-1, arr_vec, k, join_merge);
  
  //the result arrangement is at index 0
  this->m_arr = arr_vec[0];
}


//join range of polygons with holes
template <class InputIterator>
inline void join(InputIterator begin,
                 InputIterator end,
                 Polygon_with_holes_2& pgn,
                 unsigned int k=5)
{
  std::vector<Arrangement_2*> arr_vec (std::distance(begin, end) + 1);
  arr_vec[0] = this->m_arr;
 
  unsigned int i = 1;
  for(InputIterator itr = begin; itr!=end; ++itr, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn_with_holes2arr(*itr, *arr_vec[i]);
  }

  Join_merge<Arrangement_2> join_merge;
  divide_and_conquer(0, arr_vec.size()-1, arr_vec, k, join_merge);

  //the result arrangement is at index 0
  this->m_arr = arr_vec[0];
}

template <class InputIterator1, class InputIterator2>
inline void join(InputIterator1 begin1,
                 InputIterator1 end1,
                 InputIterator2 begin2,
                 InputIterator2 end2,
                 unsigned int k=5)
{
  std::vector<Arrangement_2*> arr_vec (std::distance(begin1, end1)+
                                       std::distance(begin2, end2)+1);
 
  arr_vec[0] = this->m_arr;
  unsigned int i = 1;

  for(InputIterator1 itr1 = begin1; itr1!=end1; ++itr1, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn2arr(*itr1, *arr_vec[i]);
  }

  for(InputIterator2 itr = begin2; itr2!=end2; ++itr2, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn_with_holes2arr(*itr2, *arr_vec[i]);
  }

  Join_merge<Arrangement_2> join_merge;
  divide_and_conquer(0, arr_vec.size()-1, arr_vec, k, join_merge);

  //the result arrangement is at index 0
  this->m_arr = arr_vec[0];
  this->remove_redundant_edges();
}


// intersect range of polygins
template <class InputIterator>
inline void intersection(InputIterator begin, InputIterator end)
{
  typename std::iterator_traits<InputIterator>::value_type pgn;
  this->intersection(begin, end, pgn);
  this->remove_redundant_edges();
}


// intersect range of simple polygons
template <class InputIterator>
inline void intersection(InputIterator begin,
                         InputIterator end,
                         Polygon_2& pgn)
{
  std::vector<Arrangement_2*> arr_vec (std::distance(begin, end) + 1);
  arr_vec[0] = this->m_arr;
  unsigned int i = 1;
 
  for(InputIterator itr = begin; itr!=end; ++itr, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn2arr(*itr, *arr_vec[i]);
  }

  Intersection_merge<Arrangement_2> intersection_merge;
  divide_and_conquer(0, arr_vec.size()-1, arr_vec, k, intersection_merge);

  //the result arrangement is at index 0
  this->m_arr = arr_vec[0];
}

//intersect range of polygons with holes
template <class InputIterator>
inline void intersection(InputIterator begin,
                         InputIterator end,
                         Polygon_with_holes_2& pgn)
{
  std::vector<Arrangement_2*> arr_vec (std::distance(begin, end) + 1);
  arr_vec[0] = this->m_arr;
  unsigned int i = 1;
 
  for(InputIterator itr = begin; itr!=end; ++itr, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn_with_holes2arr(*itr, *arr_vec[i]);
  }

  Intersection_merge<Arrangement_2> intersection_merge;
  divide_and_conquer(0, arr_vec.size()-1, arr_vec, k, intersection_merge);

  //the result arrangement is at index 0
  this->m_arr = arr_vec[0];
}


template <class InputIterator1, class InputIterator2>
inline void intersection(InputIterator1 begin1,
                         InputIterator1 end1,
                         InputIterator2 begin2,
                         InputIterator2 end2)
{
  std::vector<Arrangement_2*> arr_vec (std::distance(begin1, end1)+
                                       std::distance(begin2, end2)+1);
  arr_vec[0] = this->m_arr;
  unsigned int i = 1;
 
  for(InputIterator1 itr1 = begin1; itr1!=end1; ++itr1, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn2arr(*itr1, *arr_vec[i]);
  }

  for(InputIterator2 itr = begin2; itr2!=end2; ++itr2, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn_with_holes2arr(*itr2, *arr_vec[i]);
  }

  Intersection_merge<Arrangement_2> intersection_merge;
  divide_and_conquer(0, arr_vec.size()-1, arr_vec, k, intersection_merge);

  //the result arrangement is at index 0
  this->m_arr = arr_vec[0];
  this->remove_redundant_edges();
}



// symmetric_difference of a range of polygons (similar to xor)
template <class InputIterator>
inline void symmetric_difference(InputIterator begin, InputIterator end)
{
  typename std::iterator_traits<InputIterator>::value_type pgn;
  this->symmetric_difference(begin, end, pgn);
  this->remove_redundant_edges();
}


// intersect range of simple polygons
template <class InputIterator>
inline void symmetric_difference(InputIterator begin,
                                 InputIterator end,
                                 Polygon_2& pgn)
{
  std::vector<Arrangement_2*> arr_vec (std::distance(begin, end) + 1);
  arr_vec[0] = this->m_arr;
  unsigned int i = 1;
 
  for(InputIterator itr = begin; itr!=end; ++itr, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn2arr(*itr, *arr_vec[i]);
  }

  Xor_merge<Arrangement_2> xor_merge;
  divide_and_conquer(0, arr_vec.size()-1, arr_vec, k, xor_merge);

  //the result arrangement is at index 0
  this->m_arr = arr_vec[0];
}

//intersect range of polygons with holes
template <class InputIterator>
inline void symmetric_difference(InputIterator begin,
                                 InputIterator end,
                                 Polygon_with_holes_2& pgn)
{
  std::vector<Arrangement_2*> arr_vec (std::distance(begin, end) + 1);
  arr_vec[0] = this->m_arr;
  unsigned int i = 1;
 
  for(InputIterator itr = begin; itr!=end; ++itr, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn_with_holes2arr(*itr, *arr_vec[i]);
  }

  Xor_merge<Arrangement_2> xor_merge;
  divide_and_conquer(0, arr_vec.size()-1, arr_vec, k, xor_merge);

  //the result arrangement is at index 0
  this->m_arr = arr_vec[0];
}


template <class InputIterator1, class InputIterator2>
inline void symmetric_difference(InputIterator1 begin1,
                                 InputIterator1 end1,
                                 InputIterator2 begin2,
                                 InputIterator2 end2)
{
  std::vector<Arrangement_2*> arr_vec (std::distance(begin1, end1)+
                                       std::distance(begin2, end2)+1);
  arr_vec[0] = this->m_arr;
  unsigned int i = 1;
 
  for(InputIterator1 itr1 = begin1; itr1!=end1; ++itr1, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn2arr(*itr1, *arr_vec[i]);
  }

  for(InputIterator2 itr = begin2; itr2!=end2; ++itr2, ++i)
  {
    arr_vec[i] = new Arrangement_2(m_traits);
    pgn_with_holes2arr(*itr2, *arr_vec[i]);
  }

  Xor_merge<Arrangement_2> xor_merge;
  divide_and_conquer(0, arr_vec.size()-1, arr_vec, k, xor_merge);

  //the result arrangement is at index 0
  this->m_arr = arr_vec[0];
  this->remove_redundant_edges();
}



static void construct_polygon(Ccb_halfedge_const_circulator ccb,
                              Polygon_2&              pgn,
                              Traits_2*               tr);


bool is_hole_of_face(Face_const_handle f,
                     Halfedge_const_handle he) const;

Ccb_halfedge_const_circulator get_boundary_of_polygon(Face_const_iterator f) const;


void remove_redundant_edges()
{
    for(Edge_iterator itr = m_arr->edges_begin(); 
      itr != m_arr->edges_end(); )
  {
    Halfedge_handle he = itr;
    if(he->face()->contained() && he->twin()->face()->contained())
    {
      Edge_iterator next = itr;
      ++next;
      m_arr->remove_edge(he);
      itr = next;
    }
    else
      ++itr;
  }
}

template <class Merge>
void divide_and_conquer(unsigned int lower,
                        unsigned int upper,
                        std::vector<Arrangement_2*>& arr_vec,
                        unsigned int k,
                        Merge merge_func)
{
  if((upper - lower) < k)
  {
    merge_func(lower, upper, 1, arr_vec);
    return;
  }

  unsigned int sub_size = ((upper - lower + 1) / k);
  
  unsigned int i=0;
  unsigned int curr_lower = lower;

  for(; i<k-1; ++i, curr_lower += sub_size )
  {
    divide_and_conquer(curr_lower,
                       curr_lower + sub_size-1,
                       arr_vec,
                       k,
                       merge_func);
  }
  divide_and_conquer(curr_lower, upper,arr_vec, k, merge_func);
  merge_func(lower, curr_lower, sub_size ,arr_vec);

  return;
}

};

#include <CGAL/Boolean_set_operations_2/Bso_utils.h>


CGAL_END_NAMESPACE

#endif
