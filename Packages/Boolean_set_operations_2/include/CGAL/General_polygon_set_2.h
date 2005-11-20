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



  General_polygon_set_2(const Polygon_2& pgn) : 
    m_traits(new Traits_2()),
    m_traits_owner(true),
    m_arr(new Arrangement_2(m_traits)) 
  {
    pgn2arr(pgn, *m_arr);
  }

  General_polygon_set_2(const Polygon_with_holes_2& pgn_with_holes): 
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
  General_polygon_set_2(PolygonIterator pgn_begin, PolygonIterator pgn_end);


  // constructor of two ranges of : the first one for simple polygons,
  // the second one for polygons with holes
  // precondition: the first range is disjoint simple polygons 
  //               the second range is fisjoint polygons with holes
  template <class PolygonIterator, class PolygonWithHolesIterator>
  General_polygon_set_2(PolygonIterator pgn_begin,
                       PolygonIterator pgn_end,
                       PolygonWithHolesIterator  pgn_with_holes_begin,
                       PolygonWithHolesIterator  pgn_with_holes_end);


  //destructor
  ~General_polygon_set_2()
  {
    if(m_traits_owner)
      delete m_arr;
  }


  // insert a simple polygon
  void insert(const Polygon_2& pgn);

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
    Arrangement_2 second_arr;
    Arrangement_2 res_arr;
    pgn2arr(pgn.vertices_begin(), pgn.vertices_end(), second_arr);

    Bso_do_intersect_functor<Traits_2>  func;
    overlay(*m_arr, second_arr, res_arr, func);
    return func.found_intersection();
  }

  // test for intersection with a polygon with holes
  bool do_intersect(const Polygon_with_holes_2& pgn_with_holes)
  {
    Arrangement_2 second_arr;
    Arrangement_2 res_arr;
    pgn_with_holes2arr(pgn_with_holes, second_arr);

    Bso_do_intersect_functor<Traits_2>  func;
    overlay(*m_arr, second_arr, res_arr, func);
    return func.found_intersection();
  }

  //test for intersection with another General_polygon_set_2 object
  bool do_intersect(const General_polygon_set_2& bops)
  {
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
    pgn_with_holes2arr(pgn.vertices_begin(), pgn.vertices_end(), second_arr);
    if(second_arr->is_empty())
    {
      if(! second_arr->unbounded_face()->contained())
        this ->clear();
     
      return;
    }


    Bso_intersection_functor<Traits_2>  func;
    overlay(*m_arr, second_arr, res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }

  //intersection with another General_polygon_set_2 object
  void intersection(const General_polygon_set_2& bops)
  {
    Arrangement_2 res_arr;

    Bso_intersection_functor<Traits_2> func(m_traits);
    overlay(*m_arr, second_arr, *(bops.m_arr), func);
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
  void join(const Polygon_with_holes_2& pgn_with_holes);

  //join with another General_polygon_set_2 object
  void join(const General_polygon_set_2& bops)
  {
    Arrangement_2 res_arr;

    Bso_join_functor<Traits_2>  func(m_traits);
    overlay(*m_arr, second_arr, *(bops.m_arr), func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }


   // difference with a simple polygon
  void difference (const Polygon_2& pgn)
  {
    boolean_operation<Bso_difference_functor<Traits_2> >(pgn);
  }

  // difference with a polygon with holes
  void difference (const Polygon_with_holes_2& pgn_with_holes);

  //difference with another General_polygon_set_2 object
  void difference (const General_polygon_set_2& bops)
  {
   Arrangement_2 res_arr;

    Bso_difference_functor<Traits_2>  func(m_traits);
    overlay(*m_arr, second_arr, *(bops.m_arr), func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }


  // symmetric_difference with a simple polygon
  void symmetric_difference(const Polygon_2& pgn)
  {
    boolean_operation<Bso_sym_diff_functor<Traits_2> >(pgn);
  }

  // symmetric_difference with a polygon with holes
  void symmetric_difference(const Polygon_with_holes_2& pgn_with_holes);

  //symmetric_difference with another General_polygon_set_2 object
  void symmetric_difference(const General_polygon_set_2& bops)
  {
    boolean_operation<Bso_sym_diff_functor<Traits_2> >(bops);
  }

  void complement()
  {
    for(Face_iterator fit = m_arr->faces_begin();
        fit != m_arr->faces_end();
        ++fit)
    {
      fit->set_contained(!fit->contained());
    }
  }




  // get the simple polygons, takes O(n)
  template <class OutputIterator>
  OutputIterator polygons(OutputIterator out);

  
  Size number_of_curved_polygons_with_holes() const
  {
    Size i = 0;
    for(Face_const_iterator fit = m_arr->faces_begin();
        fit != m_arr->faces_end();
        ++fit)
    {
      if(fit->contained())
        ++i;
    }

    return i;
  }

  std::pair<Size, Size> number_of_polygons() const
  {
    Size i1 = 0;
    Size i2 = 0;

    for(Face_const_iterator fit = m_arr->faces_begin();
        fit != m_arr->faces_end();
        ++fit)
    {
      if(fit->contained())
      {
        if(is_simple(fit))
          ++i1;
        else
          ++i2;
      }
    }

    return std::make_pair(i1, i2);
  }


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

  
  // Returns the symbolic constant ON_BOUNDED_SIDE, ON_BOUNDARY or
  // ON_UNBOUNDED_SIDE, depending on where point q is.
  Bounded_side  bounded_side(const Point_2& q)
  {
    Walk_pl pl(*m_arr);

    Object obj = pl.locate(q);
    Face_const_iterator f;
    if(CGAL::assign(f, obj))
    {
      if(f->contained())
        return ON_BOUNDED_SIDE;
      else
        return ON_UNBOUNDED_SIDE;
    }
    else
    {
      return ON_BOUNDARY;
    }
  }


  // returns the location of the query point
  Object locate(const Point_2& q);

  //advanced function: get const reference to the arrangement
  const Arrangement_2& arrangement() const
  {
    return *m_arr;
  }

  bool is_valid() const
  {
    if(! m_arr->is_valid())
      return false;

    for(Edge_const_iterator eci = m_arr->edges_begin();
        eci != m_arr->edges_end();
        ++eci)
    {
      Haldedge_handle he = eci;
      if(he->face() == he->twin()->face()
        return false;
      if(he->face()->contained() == he->twin()->face()->contained)
        return false;
    }
    return true;
  }

  private:
  
  template <class Bops_functor>
  void boolean_operation(const Polygon_2& pgn)
  {
    Arrangement_2 second_arr;
    Arrangement_2 res_arr;
    pgn2arr(pgn, second_arr);

    Bops_functor  func(m_traits);
    overlay(*m_arr, second_arr, res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }


  template <class Bops_functor>
  void boolean_operation(const General_polygon_set_2& bops)
  {
    Arrangement_2 res_arr;

    Bops_functor  func(m_traits);
    overlay(*m_arr, *(bops.m_arr), res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = func.result_arr();
  }

  public:

  static void pgn2arr (const Polygon_2& pgn, Arrangement_2& arr);

  template< class PolygonIter >
  void pgns2arr(PolygonIter p_begin, PolygonIter p_end, Arrangement_2& arr);
  
  void pgn_with_holes2arr (const Polygon_2& pgn, Arrangement_2& arr);

  
  template< class InputIterator >
  void pgns_with_holes2arr (InputIterator begin, 
                            InputIterator end,
                            Arrangement_2& arr);



  public:
  static void construct_polygon(Ccb_halfedge_circulator ccb,
                         Polygon_2& pgn,
                         Traits_2* tr);

};

#include <CGAL/Boolean_set_operations_2/Bso_utils.h>

template <class PolygonIter, class OutputIterator, class Traits>
OutputIterator aggregate_union(PolygonIter begin,
                               PolygonIter end,
                               OutputIterator out,
                               Traits tr)
{
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Polygon_2            Polygon_2;
  typedef typename Traits::Polygon_with_holes_2  
                                                 Polygon_with_holes_2;
  typedef typename Polygon_with_holes_2::Holes_const_iterator  
                                                 GP_Holes_const_iterator;

  typedef Bso_dcel<Traits>                                Bso_dcel;  
  typedef Arrangement_2<Traits, Bso_dcel>                 Arrangement_2;
  typedef typename Arrangement_2::Edge_iterator            Edge_iterator;
  typedef typename Arrangement_2::Face_iterator            Face_iterator;

  std::vector<Arrangement_2*> arr_vec (std::distance(begin, end));
 
  unsigned int i = 0;
  for(PolygonIter itr = begin; itr!=end; ++itr, ++i)
  {
    arr_vec[i] = new Arrangement_2(&tr);
    General_polygon_set_2<Traits>::pgn2arr(*itr, *arr_vec[i]);
  }

  aggregate_union_rec(0, arr_vec.size()-1, arr_vec);

  std::cout<<"|V| = " <<arr_vec[0]->number_of_vertices()<<"\n";
  std::cout<<"|E| = " <<arr_vec[0]->number_of_edges()<<"\n";
  std::cout<<"|F| = " <<arr_vec[0]->number_of_faces()<<"\n";
  
  std::ofstream outFile("union.arr");
  // Creates an ofstream object named outFile
  if (! outFile.is_open()) // Always test file open
  {
    std::cout << "Error opening input file" << std::endl;
    return out;
  }
  outFile << (*arr_vec[0]);
  return out;
  /*typedef Construct_polygons_visitor<Arrangement_2,
                                     OutputIterator>     My_visitor;
  typedef Arr_bfs_scanner<Arrangement_2, My_visitor>     Arr_bfs_scanner;

  My_visitor v(tr, out);
  Arr_bfs_scanner scanner(v);
  scanner.scan(*arr_vec[0]);

  return (v.output_iterator());*/


}

template <class Arrangement>
void aggregate_union_rec(unsigned int lower,
                         unsigned int upper,
                         std::vector<Arrangement*>& arr_vec)
{
  if(lower == upper)
  {
    return;
  }

  if(upper - lower == 1)
  {
    union_merge(lower, upper, arr_vec);
    return;
  }
  unsigned int mid = (upper + lower)/2;
  aggregate_union_rec(lower, mid, arr_vec);
  aggregate_union_rec(mid + 1, upper, arr_vec);

  union_merge(lower, mid+1, arr_vec);

  return;
}


template <class Arrangement>
void union_merge(unsigned int i,
                 unsigned int j,
                 std::vector<Arrangement*>& arr_vec)
{
  typedef typename Arrangement::Traits_2    Traits_2;
  Bso_join_functor<Traits_2>  func(arr_vec[i]->get_traits());
  Arrangement res;
  overlay(*(arr_vec[i]), *(arr_vec[j]), res, func);
  delete arr_vec[i]; // delete the previous arrangement
  delete arr_vec[j]; // delete the previous arrangement

  arr_vec[i] = func.result_arr();
}
 

//////////////////////////////////////////////////////////////////////



template <class PolygonIter, class OutputIterator, class Traits>
OutputIterator aggregate_intersection(PolygonIter begin,
                                       PolygonIter end,
                                       OutputIterator out,
                                       Traits tr)
{
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Polygon_2            Polygon_2;
  typedef typename Traits::Polygon_with_holes_2  
                                                 Polygon_with_holes_2;
  typedef typename Polygon_with_holes_2::Holes_const_iterator  
                                                 GP_Holes_const_iterator;

  typedef Bso_dcel<Traits>                                Bso_dcel;  
  typedef Arrangement_2<Traits, Bso_dcel>                 Arrangement_2;
  typedef typename Arrangement_2::Edge_iterator            Edge_iterator;
  typedef typename Arrangement_2::Face_iterator            Face_iterator;

  std::vector<Arrangement_2*> arr_vec (std::distance(begin, end));
 
  unsigned int i = 0;
  for(PolygonIter itr = begin; itr!=end; ++itr, ++i)
  {
    arr_vec[i] = new Arrangement_2(&tr);
    General_polygon_set_2<Traits>::pgn2arr(*itr, *arr_vec[i]);
  }

  aggregate_intersection_rec(0, arr_vec.size()-1, arr_vec);

  std::cout<<"|V| = " <<arr_vec[0]->number_of_vertices()<<"\n";
  std::cout<<"|E| = " <<arr_vec[0]->number_of_edges()<<"\n";
  std::cout<<"|F| = " <<arr_vec[0]->number_of_faces()<<"\n";
  
  return out;
  /*typedef Construct_polygons_visitor<Arrangement_2,
                                     OutputIterator>     My_visitor;
  typedef Arr_bfs_scanner<Arrangement_2, My_visitor>     Arr_bfs_scanner;

  My_visitor v(tr, out);
  Arr_bfs_scanner scanner(v);
  scanner.scan(*arr_vec[0]);

  return (v.output_iterator());*/


}

template <class Arrangement>
void aggregate_intersection_rec(unsigned int lower,
                                 unsigned int upper,
                                 std::vector<Arrangement*>& arr_vec)
{
  if(lower == upper)
  {
    return;
  }

  if(upper - lower == 1)
  {
    intersection_merge(lower, upper, arr_vec);
    return;
  }
  unsigned int mid = (upper + lower)/2;
  aggregate_intersection_rec(lower, mid, arr_vec);
  aggregate_intersection_rec(mid + 1, upper, arr_vec);

  intersection_merge(lower, mid+1, arr_vec);

  return;
}


template <class Arrangement>
void intersection_merge(unsigned int i,
                        unsigned int j,
                        std::vector<Arrangement*>& arr_vec)
{
  typedef typename Arrangement::Traits_2    Traits_2;
  
  if(arr_vec[i]->is_empty() || arr_vec[j]->is_empty())
    delete arr_vec[j];
  else
  {
    Bso_intersection_functor<Traits_2>  func(arr_vec[i]->get_traits());
    Arrangement res;
    overlay(*(arr_vec[i]), *(arr_vec[j]), res, func);
    delete arr_vec[i]; // delete the previous arrangement
    delete arr_vec[j]; // delete the previous arrangement
    arr_vec[i] = func.result_arr();
  }
}


CGAL_END_NAMESPACE

#endif
