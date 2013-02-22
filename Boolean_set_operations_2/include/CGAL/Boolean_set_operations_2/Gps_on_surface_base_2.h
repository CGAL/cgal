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
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ophir Setter    <ophir.setter@cs.tau.ac.il>
//                 Guy Zucker <guyzucke@post.tau.ac.il> 


#ifndef CGAL_GPS_ON_SURFACE_BASE_2_H
#define CGAL_GPS_ON_SURFACE_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/enum.h>
#include <CGAL/iterator.h> 
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_do_intersect_functor.h>
#include <CGAL/Boolean_set_operations_2/Gps_intersection_functor.h>
#include <CGAL/Boolean_set_operations_2/Gps_join_functor.h>
#include <CGAL/Boolean_set_operations_2/Gps_difference_functor.h>
#include <CGAL/Boolean_set_operations_2/Gps_sym_diff_functor.h>
#include <CGAL/Boolean_set_operations_2/Gps_merge.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_simplifier.h>
#include <CGAL/Boolean_set_operations_2/Ccb_curve_iterator.h>

/*!
  \file   Gps_on_surface_base_2.h
  \brief  A class that allows Boolean set operations.
  This class is the base class for General_polygon_set_on_surface_2 and 
  recieves extra template parameter which allows different validation 
  policies. If you do not want validation then use the default validation
  policy. A different validation policy example can be found in 
  General_polygon_set_on_surface_2.
*/


namespace CGAL {

namespace Boolean_set_operation_2_internal
{
  struct NoValidationPolicy
  {
   /*! is_valid - Checks if a Traits::Polygon_2 OR 
    * Traits::Polygon_with_holes_2 are valid.
    * In this validation policy we do NOT do anything.
    */ 
    template <class Polygon, class Traits>
    inline static void is_valid(const Polygon&, Traits&) {}
  };
}

//! General_polygon_set_on_surface_2 
/*! This class is the base class for General_polygon_set_on_surface_2 and 
    recieves extra template parameter which allows different validation 
    policies. If you do not want validation then use the default validation
    policy. A different validation policy example can be found in 
    General_polygon_set_on_surface_2.
 */
template <class Traits_, class TopTraits_, 
          class ValidationPolicy = 
          Boolean_set_operation_2_internal::NoValidationPolicy>
class Gps_on_surface_base_2
{
public:
  typedef Traits_                                      Traits_2;
  typedef TopTraits_                                   Topology_traits;
  typedef typename Traits_2::Polygon_2                 Polygon_2;
  typedef typename Traits_2::Polygon_with_holes_2      Polygon_with_holes_2;
  typedef CGAL::Arrangement_on_surface_2<Traits_2, Topology_traits>
                                                       Arrangement_on_surface_2;
  typedef typename Arrangement_on_surface_2::Size      Size;

private:
  typedef Arrangement_on_surface_2                     Aos_2;

  typedef Gps_on_surface_base_2 <
    Traits_2, Topology_traits, ValidationPolicy>       Self;
  typedef typename Traits_2::Point_2                   Point_2;
  typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
  
  typedef typename Polygon_with_holes_2::Hole_const_iterator
    GP_Holes_const_iterator;
  typedef typename Traits_2::Curve_const_iterator      Curve_const_iterator;
  typedef typename Traits_2::Compare_endpoints_xy_2
                                                       Compare_endpoints_xy_2;
  typedef typename Traits_2::Construct_opposite_2      Construct_opposite_2;

  typedef typename Aos_2::Face_const_iterator          Face_const_iterator;
  typedef typename Aos_2::Halfedge_const_iterator      Halfedge_const_iterator;
  typedef typename Aos_2::Vertex_const_iterator        Vertex_const_iterator;
  typedef typename Aos_2::Edge_const_iterator          Edge_const_iterator;
  typedef typename Aos_2::Outer_ccb_const_iterator     Outer_ccb_const_iterator;
  typedef typename Aos_2::Inner_ccb_const_iterator     Inner_ccb_const_iterator;
  typedef typename Aos_2::Ccb_halfedge_const_circulator 
    Ccb_halfedge_const_circulator;
  typedef typename Aos_2::Face_iterator                Face_iterator;
  typedef typename Aos_2::Halfedge_iterator            Halfedge_iterator;
  typedef typename Aos_2::Vertex_iterator              Vertex_iterator;
  typedef typename Aos_2::Edge_iterator                Edge_iterator;
  typedef typename Aos_2::Outer_ccb_iterator           Outer_ccb_iterator;
  typedef typename Aos_2::Inner_ccb_iterator           Inner_ccb_iterator;
  typedef typename Aos_2::Ccb_halfedge_circulator      Ccb_halfedge_circulator;
  typedef typename Aos_2::Face_handle                  Face_handle;
  typedef typename Aos_2::Halfedge_handle              Halfedge_handle;
  typedef typename Aos_2::Vertex_handle                Vertex_handle;

  typedef typename Aos_2::Face_const_handle            Face_const_handle;
  typedef typename Aos_2::Halfedge_const_handle        Halfedge_const_handle;
  typedef typename Aos_2::Vertex_const_handle          Vertex_const_handle;

  typedef typename Aos_2::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;

  typedef std::pair<Aos_2 *, 
                    std::vector<Vertex_handle> *>      Arr_entry;

  typedef typename Arrangement_on_surface_2::
    Topology_traits::Default_point_location_strategy   Point_location;

protected:

  // Traits* should be removed and only m_traits should be used.
  // If you, who reads this text, have time, replace m_traits 
  // with m_traits_adaptor and try to do something about m_traits_owner.
  Traits_2*                                  m_traits;
  CGAL::Arr_traits_adaptor_2<Traits_2>       m_traits_adaptor;
  bool                                       m_traits_owner;

  // the underlying arrangement
  Aos_2*        m_arr;


public:

  // default costructor
  Gps_on_surface_base_2() : m_traits(new Traits_2()),
                            m_traits_adaptor(*m_traits),
                            m_traits_owner(true),
                            m_arr(new Aos_2(m_traits))       
  {}


  // constructor with traits object
  Gps_on_surface_base_2(Traits_2& tr) : m_traits(&tr),
                                        m_traits_adaptor(*m_traits),
                                        m_traits_owner(false),
                                        m_arr(new Aos_2(m_traits)) 
  {}


  Gps_on_surface_base_2(const Self& ps) :
    m_traits(new Traits_2(*(ps.m_traits))),
    m_traits_adaptor(*m_traits),
    m_traits_owner(true),
    m_arr(new Aos_2(*(ps.m_arr)))
  {}

  
  Gps_on_surface_base_2& operator=(const Self& ps)
  {
    if (this == &ps)
      return (*this);

    if (m_traits_owner)
      delete m_traits;
    delete m_arr;
    m_traits = new Traits_2(*(ps.m_traits));
    m_traits_adaptor = CGAL::Arr_traits_adaptor_2<Traits_2>(*m_traits);
    m_traits_owner = true;
    m_arr = new Aos_2(*(ps.m_arr));
    return (*this);
  }


  explicit Gps_on_surface_base_2(const Polygon_2& pgn) : 
    m_traits(new Traits_2()),
    m_traits_adaptor(*m_traits),
    m_traits_owner(true),
    m_arr(new Aos_2(m_traits)) 
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _insert(pgn, *m_arr);
  }

  explicit Gps_on_surface_base_2(const Polygon_with_holes_2& pgn_with_holes): 
    m_traits(new Traits_2()),
    m_traits_adaptor(*m_traits),
    m_traits_owner(true),
    m_arr(new Aos_2(m_traits))
  {
    ValidationPolicy::is_valid(pgn_with_holes,*m_traits);
    _insert(pgn_with_holes, *m_arr);
  }

protected:
  Gps_on_surface_base_2(Aos_2* arr) : m_traits(new Traits_2()),
                                              m_traits_adaptor(*m_traits),
                                              m_traits_owner(true),
                                              m_arr(arr)
   {}

public:
  //destructor
  virtual ~Gps_on_surface_base_2()
  {
    delete m_arr;

    if (m_traits_owner)
      delete m_traits;
  }

  void simplify(const Polygon_2& pgn, Polygon_with_holes_2& res)
  {
    typedef Gps_polygon_simplifier<Aos_2>  Simplifier;

    Aos_2*  arr = new Aos_2();

    Simplifier simp(*arr, *m_traits);
    simp.simplify(pgn);
    _remove_redundant_edges(arr);
    Self gps(arr);
    gps._reset_faces();
  
    typedef Oneset_iterator<Polygon_with_holes_2>    OutputItr;
    OutputItr oi (res);
    gps.polygons_with_holes(oi);
  }

  // insert a simple polygon
  void insert(const Polygon_2& pgn)
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _insert(pgn, *m_arr);
  }

  // insert a polygon with holes
  void insert(const Polygon_with_holes_2& pgn_with_holes)
  {
    ValidationPolicy::is_valid(pgn_with_holes, *m_traits);
    _insert(pgn_with_holes, *m_arr);
  }
  
  // insert a range of polygons that can be either simple polygons
  // or polygons with holes
  // precondition: the polygons are disjoint and simple
  template <typename PolygonIterator>
  void insert(PolygonIterator pgn_begin, PolygonIterator pgn_end);


  // insert two ranges of : the first one for simple polygons,
  // the second one for polygons with holes
  // precondition: the first range is disjoint simple polygons 
  //               the second range is disjoint polygons with holes
  template <typename PolygonIterator, typename PolygonWithHolesIterator>
  void insert(PolygonIterator pgn_begin, PolygonIterator pgn_end,
              PolygonWithHolesIterator pgn_with_holes_begin,
              PolygonWithHolesIterator pgn_with_holes_end);

  // test for intersection with a simple polygon
  bool do_intersect(const Polygon_2 &pgn) const
  {
    ValidationPolicy::is_valid(pgn,*m_traits);
    Self other(pgn);
    return (do_intersect(other));
  }

  // test for intersection with a polygon with holes
  bool do_intersect(const Polygon_with_holes_2& pgn_with_holes) const
  {
    ValidationPolicy::is_valid(pgn_with_holes, *m_traits);
    Self other(pgn_with_holes);
    return (do_intersect(other));
  }

  //test for intersection with another Gps_on_surface_base_2 object
  bool do_intersect(const Self& other) const
  {
    if (this->is_empty() || other.is_empty())
      return false;

    if (this->is_plane() || other.is_plane())
      return true;
    
    Aos_2 res_arr;

    Gps_do_intersect_functor<Aos_2>  func;
    overlay(*m_arr, *(other.m_arr), res_arr, func);
    return func.found_reg_intersection();
  }

  // intersection with a simple polygon
  void intersection(const Polygon_2& pgn)
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _intersection(pgn);
  }

  // intersection with a polygon with holes
  void intersection(const Polygon_with_holes_2& pgn)
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _intersection(pgn);
  }

  //intersection with another Gps_on_surface_base_2 object
  void intersection(const Self& other)
  {
    _intersection(other);
  }

  void intersection(const Self& gps1, const Self& gps2)
  {
    this->clear();
    _intersection(*(gps1.m_arr), *(gps2.m_arr), *(this->m_arr));
  }


  // join with a simple polygon
  void join(const Polygon_2& pgn)
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _join(pgn);
  }

  // join with a polygon with holes
  void join(const Polygon_with_holes_2& pgn)
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _join(pgn);
  }

  //join with another Gps_on_surface_base_2 object
  void join(const Self& other)
  {
    _join(other);
  }

  void join(const Self& gps1, const Self& gps2)
  {
    this->clear();
    _join(*(gps1.m_arr), *(gps2.m_arr), *(this->m_arr));
  }

  // difference with a simple polygon
  void difference (const Polygon_2& pgn)
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _difference(pgn);
  }

  // difference with a polygon with holes
  void difference (const Polygon_with_holes_2& pgn)
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _difference(pgn);
  }
  
  //difference with another Gps_on_surface_base_2 object
  void difference (const Self& other)
  {
    _difference(other);
  }

  void difference(const Self& gps1, const Self& gps2)
  {
    this->clear();
    _difference(*(gps1.m_arr), *(gps2.m_arr), *(this->m_arr));
  }


  // symmetric_difference with a simple polygon
  void symmetric_difference(const Polygon_2& pgn)
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _symmetric_difference(pgn);
  }

  // symmetric_difference with a polygon with holes
  void symmetric_difference(const Polygon_with_holes_2& pgn)
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _symmetric_difference(pgn);
  }

  //symmetric_difference with another Gps_on_surface_base_2 object
  void symmetric_difference(const Self& other)
  {
    _symmetric_difference(other);
  }

  void symmetric_difference(const Self& gps1, const Self& gps2)
  {
    this->clear();
    _symmetric_difference(*(gps1.m_arr), *(gps2.m_arr), *(this->m_arr));
  }


  void complement()
  {
    this->_complement(m_arr);
  }

  void complement(const Self& other)
  {
    *(this->m_arr) = *(other.m_arr);
    this->complement();
  }

  void fix_curves_direction()
  {
    _fix_curves_direction(*m_arr);
  }
         
  Size number_of_polygons_with_holes() const;

  Traits_2& traits()
  {
    return *m_traits;
  }

  const Traits_2& traits() const
  {
    return *m_traits;
  }

  bool is_empty() const
  {
    // We have to check that all the faces of an empty arrangement are not
    // conained in the polygon set (there can be several faces in an empty
    // arrangement, dependant on the topology traits.
    // The point is that if the arrangement is "empty" (meaning that no curve
    // or point were inserted and that it is in its original state) then 
    // all the faces (created by the topology traits) should have the same
    // result for contained() --- from Boolean operations point of view there
    // can not be an empty arrangement which has serveral faces with different
    // attributes.
    return (m_arr->is_empty() && !m_arr->faces_begin()->contained());
  }

  bool is_plane() const
  {
    // Same comment as in "is_empty" above, just with adjustments.
    return (m_arr->is_empty() &&  m_arr->faces_begin()->contained());
  }

  void clear()
  {
    m_arr->clear();
  }

  
  Oriented_side oriented_side(const Point_2& q) const
  {
    Point_location pl(*m_arr);

    Object obj = pl.locate(q);
    Face_const_iterator f;
    if (CGAL::assign(f, obj))
    {
      if (f->contained())
        return ON_POSITIVE_SIDE;

      return ON_NEGATIVE_SIDE ;
    }  
    return ON_ORIENTED_BOUNDARY ;
  }

  Oriented_side oriented_side(const Polygon_2& pgn) const
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    Self other(pgn);
    return (oriented_side(other));
  }

  Oriented_side oriented_side(const Polygon_with_holes_2& pgn) const
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    Self other(pgn);
    return (oriented_side(other));
  }

  Oriented_side oriented_side(const Self& other) const
  {
    if (this->is_empty() || other.is_empty())
      return ON_NEGATIVE_SIDE;

    if (this->is_plane() || other.is_plane())
      return ON_POSITIVE_SIDE;
    
    Aos_2 res_arr;

    Gps_do_intersect_functor<Aos_2>  func;
    overlay(*m_arr, *(other.m_arr), res_arr, func);
    if (func.found_reg_intersection())
      return ON_POSITIVE_SIDE;
    
    if (func.found_boundary_intersection())
      return ON_ORIENTED_BOUNDARY;

    return ON_NEGATIVE_SIDE;
  }


  // returns the location of the query point
  bool locate(const Point_2& q, Polygon_with_holes_2& pgn) const;

  /*! Obtain a const reference to the underlying arrangement
   * \return the underlying arrangement.
   */
  const Aos_2& arrangement() const
  {
    return *m_arr;
  }

  /*! Obtain a reference to the underlying arrangement
   * \return the underlying arrangement.
   */
  Aos_2& arrangement()
  {
    return *m_arr;
  }
  
protected:
  
  bool _is_valid(Aos_2& arr) {
    if (!CGAL::is_valid(arr))
      return false;

    Compare_endpoints_xy_2 cmp_endpoints =
      m_traits->compare_endpoints_xy_2_object();

    for (Edge_const_iterator eci = arr.edges_begin();
         eci != arr.edges_end();
         ++eci)
    {
      Halfedge_const_handle he = eci;
      if (he->face() == he->twin()->face())
      {
        return false;
      }
      if (he->face()->contained() == he->twin()->face()->contained())
      {
        return false;
      }

      const X_monotone_curve_2&  cv = he->curve();
      const bool                 is_cont = he->face()->contained();
      const Comparison_result    he_res = 
        ((Arr_halfedge_direction)he->direction() == ARR_LEFT_TO_RIGHT) ?
        SMALLER : LARGER;
      const bool                 has_same_dir = (cmp_endpoints(cv) == he_res);

      if ((is_cont && !has_same_dir) || (!is_cont && has_same_dir)) {
        return false;
      }
    }
    return true;
  }

public:

  /*! */
  bool is_valid()
  {
    return _is_valid(*this->m_arr);
  }

  // get the simple polygons, takes O(n)
  template <typename OutputIterator>
  OutputIterator polygons_with_holes(OutputIterator out) const;

  // test for intersection of a range of polygons
  template <typename InputIterator>
  bool do_intersect(InputIterator begin, InputIterator end, unsigned int k = 5)
  {
    Self other(*this);
    other.intersection(begin, end, k);
    return (other.is_empty());
  }

  template <typename InputIterator1, typename InputIterator2>
  bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                    InputIterator2 begin2, InputIterator2 end2,
                    unsigned int k = 5)
  {
    Self other(*this);
    other.intersection(begin1, end1, begin2, end2, k);
    return (other.is_empty());
  }

  // join a range of polygons
  template <typename InputIterator>
  void join(InputIterator begin, InputIterator end, unsigned int k = 5)
  {
    typename std::iterator_traits<InputIterator>::value_type pgn;
    this->join(begin, end, pgn, k);
    this->remove_redundant_edges();
    this->_reset_faces();
  }

  // join range of simple polygons
  // 5 is the magic number in which we switch to a sweep-based algorithm
  // instead of a D&C algorithm. This point should be further studies, as
  // it is hard to believe that this is the best value for all applications.
  template <typename InputIterator>
  inline void join(InputIterator begin, InputIterator end, Polygon_2&,
                   unsigned int k = 5)
  {
    std::vector<Arr_entry> arr_vec (std::distance(begin, end) + 1);
 
    arr_vec[0].first = this->m_arr;
    unsigned int i = 1;
    for (InputIterator itr = begin; itr != end; ++itr, ++i)
    {
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr, *(arr_vec[i].first));
    }

    Join_merge<Aos_2> join_merge;
    _build_sorted_vertices_vectors (arr_vec);
    _divide_and_conquer(0, static_cast<unsigned int>(arr_vec.size()-1), arr_vec, k, join_merge);
  
    //the result arrangement is at index 0
    this->m_arr = arr_vec[0].first;
    delete arr_vec[0].second;
  }

  //join range of polygons with holes (see previous comment about k=5).
  template <typename InputIterator>
  inline void join(InputIterator begin, InputIterator end,
                   Polygon_with_holes_2&, unsigned int k = 5)
  {
    std::vector<Arr_entry> arr_vec (std::distance(begin, end) + 1);
    arr_vec[0].first = this->m_arr;
 
    unsigned int i = 1;
    for (InputIterator itr = begin; itr!=end; ++itr, ++i)
    {
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr, *(arr_vec[i].first));
    }

    Join_merge<Aos_2> join_merge;
    _build_sorted_vertices_vectors (arr_vec);
    _divide_and_conquer(0, static_cast<unsigned int>(arr_vec.size()-1), arr_vec, k, join_merge);

    //the result arrangement is at index 0
    this->m_arr = arr_vec[0].first;
    delete arr_vec[0].second;
  }

  // (see previous comment about k=5).
  template <typename InputIterator1, typename InputIterator2>
  inline void join(InputIterator1 begin1, InputIterator1 end1,
                   InputIterator2 begin2, InputIterator2 end2,
                   unsigned int k = 5)
  {
    std::vector<Arr_entry> arr_vec (std::distance(begin1, end1)+
                                    std::distance(begin2, end2)+1);
 
    arr_vec[0].first = this->m_arr;
    unsigned int i = 1;
    
    for (InputIterator1 itr1 = begin1; itr1!=end1; ++itr1, ++i)
    {
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr1, *(arr_vec[i].first));
    }
    
    for (InputIterator2 itr2 = begin2; itr2!=end2; ++itr2, ++i)
    {
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr2, *(arr_vec[i].first));
    }

    Join_merge<Aos_2> join_merge;
    _build_sorted_vertices_vectors (arr_vec);
    _divide_and_conquer(0, static_cast<unsigned int>(arr_vec.size()-1), arr_vec, k, join_merge);

    //the result arrangement is at index 0
    this->m_arr = arr_vec[0].first;
    delete arr_vec[0].second;
    this->remove_redundant_edges();
    this->_reset_faces();
  }


  // intersect range of polygins (see previous comment about k=5).
  template <typename InputIterator>
  inline void intersection(InputIterator begin, InputIterator end,
                           unsigned int k = 5)
  {
    typename std::iterator_traits<InputIterator>::value_type pgn;
    this->intersection(begin, end, pgn, k);
    this->remove_redundant_edges();
    this->_reset_faces();
  }
  
  
  // intersect range of simple polygons
  template <typename InputIterator>
  inline void intersection(InputIterator begin, InputIterator end,
                           Polygon_2&, unsigned int k)
  {
    std::vector<Arr_entry> arr_vec (std::distance(begin, end) + 1);
    arr_vec[0].first = this->m_arr;
    unsigned int i = 1;
    
    for (InputIterator itr = begin; itr!=end; ++itr, ++i)
    {
      ValidationPolicy::is_valid((*itr), *m_traits);
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr, *(arr_vec[i].first));
    }
    
    Intersection_merge<Aos_2> intersection_merge;
    _build_sorted_vertices_vectors (arr_vec);
    _divide_and_conquer(0, static_cast<unsigned int>(arr_vec.size()-1), arr_vec, k, intersection_merge);
    
    //the result arrangement is at index 0
    this->m_arr = arr_vec[0].first;
    delete arr_vec[0].second;
  }
  
  //intersect range of polygons with holes
  template <typename InputIterator>
  inline void intersection(InputIterator begin, InputIterator end,
                           Polygon_with_holes_2&, unsigned int k)
  {
    std::vector<Arr_entry> arr_vec (std::distance(begin, end) + 1);
    arr_vec[0].first = this->m_arr;
    unsigned int i = 1;
    
    for (InputIterator itr = begin; itr!=end; ++itr, ++i)
    {
      ValidationPolicy::is_valid((*itr), *m_traits);
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr, *(arr_vec[i].first));
    }
    
    Intersection_merge<Aos_2> intersection_merge;
    _build_sorted_vertices_vectors (arr_vec);
    _divide_and_conquer(0, static_cast<unsigned int>(arr_vec.size()-1), arr_vec, k, intersection_merge);
    
    //the result arrangement is at index 0
    this->m_arr = arr_vec[0].first;
    delete arr_vec[0].second;
  }
  
  
  template <typename InputIterator1, typename InputIterator2>
  inline void intersection(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2,
                           unsigned int k = 5)
  {
    std::vector<Arr_entry> arr_vec (std::distance(begin1, end1)+
                                    std::distance(begin2, end2)+1);
    arr_vec[0].first = this->m_arr;
    unsigned int i = 1;
    
    for (InputIterator1 itr1 = begin1; itr1!=end1; ++itr1, ++i)
    {
      ValidationPolicy::is_valid(*itr1, *m_traits);
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr1, *(arr_vec[i].first));
    }
    
    for (InputIterator2 itr2 = begin2; itr2!=end2; ++itr2, ++i)
    {
      ValidationPolicy::is_valid(*itr2,*m_traits);
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr2, *(arr_vec[i].first));
    }
    
    Intersection_merge<Aos_2> intersection_merge;
    _build_sorted_vertices_vectors (arr_vec);
    _divide_and_conquer(0, static_cast<unsigned int>(arr_vec.size()-1), arr_vec, k, intersection_merge);
    
    //the result arrangement is at index 0
    this->m_arr = arr_vec[0].first;
    delete arr_vec[0].second;
    this->remove_redundant_edges();
    this->_reset_faces();
  }
  
  
  
  // symmetric_difference of a range of polygons (similar to xor)
  // (see previous comment about k=5).
  template <typename InputIterator>
    inline void symmetric_difference(InputIterator begin, InputIterator end,
                                     unsigned int k = 5)
  {
    typename std::iterator_traits<InputIterator>::value_type pgn;
    this->symmetric_difference(begin, end, pgn, k);
    this->remove_redundant_edges();
    this->_reset_faces();
  }
  
  
  // intersect range of simple polygons (see previous comment about k=5).
  template <typename InputIterator>
  inline void symmetric_difference(InputIterator begin, InputIterator end,
                                   Polygon_2&, unsigned int k = 5)
  {
    std::vector<Arr_entry> arr_vec (std::distance(begin, end) + 1);
    arr_vec[0].first = this->m_arr;
    unsigned int i = 1;
    
    for (InputIterator itr = begin; itr!=end; ++itr, ++i)
    {
      ValidationPolicy::is_valid(*itr,*m_traits);
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr, *(arr_vec[i].first));
    }
    
    Xor_merge<Aos_2> xor_merge;
    _build_sorted_vertices_vectors (arr_vec);
    _divide_and_conquer(0, static_cast<unsigned int>(arr_vec.size()-1), arr_vec, k, xor_merge);
    
    //the result arrangement is at index 0
    this->m_arr = arr_vec[0].first;
    delete arr_vec[0].second;
  }
  
  //intersect range of polygons with holes (see previous comment about k=5).
  template <typename InputIterator>
    inline void symmetric_difference(InputIterator begin, InputIterator end,
                                     Polygon_with_holes_2&, unsigned int k = 5)
  {
    std::vector<Arr_entry> arr_vec (std::distance(begin, end) + 1);
    arr_vec[0].first = this->m_arr;
    unsigned int i = 1;
    
    for (InputIterator itr = begin; itr!=end; ++itr, ++i)
    {
      ValidationPolicy::is_valid(*itr,*m_traits);
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr, *(arr_vec[i].first));
    }
    
    Xor_merge<Aos_2> xor_merge;
    _build_sorted_vertices_vectors (arr_vec);
    _divide_and_conquer(0, static_cast<unsigned int>(arr_vec.size()-1), arr_vec, k, xor_merge);
    
    //the result arrangement is at index 0
    this->m_arr = arr_vec[0].first;
    delete arr_vec[0].second;
  }
  
  // (see previous comment about k=5).
  template <typename InputIterator1, typename InputIterator2>
  inline void symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                                   InputIterator2 begin2, InputIterator2 end2,
                                   unsigned int k = 5)
  {
    std::vector<Arr_entry> arr_vec (std::distance(begin1, end1)+
                                    std::distance(begin2, end2)+1);
    arr_vec[0].first = this->m_arr;
    unsigned int i = 1;
    
    for (InputIterator1 itr1 = begin1; itr1!=end1; ++itr1, ++i)
    {
      ValidationPolicy::is_valid(*itr1, *m_traits);
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr1, *(arr_vec[i].first));
    }
    
    for (InputIterator2 itr2 = begin2; itr2!=end2; ++itr2, ++i)
    {
      ValidationPolicy::is_valid(*itr2, *m_traits);
      arr_vec[i].first = new Aos_2(m_traits);
      _insert(*itr2, *(arr_vec[i].first));
    }
    
    Xor_merge<Aos_2> xor_merge;
    _build_sorted_vertices_vectors (arr_vec);
    _divide_and_conquer(0, static_cast<unsigned int>(arr_vec.size()-1), arr_vec, k, xor_merge);
    
    //the result arrangement is at index 0
    this->m_arr = arr_vec[0].first;
    delete arr_vec[0].second;
    this->remove_redundant_edges();
    this->_reset_faces();
  }
  
  static void construct_polygon(Ccb_halfedge_const_circulator ccb,
                                Polygon_2 & pgn, Traits_2 * tr);
  
  bool is_hole_of_face(Face_const_handle f, Halfedge_const_handle he) const;
  
  Ccb_halfedge_const_circulator
  get_boundary_of_polygon(Face_const_iterator f) const;
  
  void remove_redundant_edges()
  {
    this->_remove_redundant_edges(m_arr);
  }
  
protected:
  
  void _remove_redundant_edges(Aos_2* arr)
  {
    for (Edge_iterator itr = arr->edges_begin(); itr != arr->edges_end(); )
    {
      Halfedge_handle he = itr;
      if (he->face()->contained() == he->twin()->face()->contained())
      {
        Edge_iterator next = itr;
        ++next;
        arr->remove_edge(he);
        itr = next;
      }
      else
        ++itr;
    }
  }
  
  class Less_vertex_handle
  {
    typename Traits_2::Compare_xy_2     comp_xy;
    
  public:
    
    Less_vertex_handle (const typename Traits_2::Compare_xy_2& cmp) :
    comp_xy (cmp)
    {}
    
    bool operator() (Vertex_handle v1, Vertex_handle v2) const
    {
      return (comp_xy (v1->point(), v2->point()) == SMALLER);
    }
  };
  

  void _complement(Aos_2* arr)
  {
    for (Face_iterator fit = arr->faces_begin();
         fit != arr->faces_end();
         ++fit)
    {
      fit->set_contained(!fit->contained());
    }

    Construct_opposite_2 ctr_opp = m_traits->construct_opposite_2_object();
    for (Edge_iterator eit = arr->edges_begin();
         eit != arr->edges_end();
         ++eit)
    {
      Halfedge_handle he = eit;
      const X_monotone_curve_2& cv = he->curve();
      arr->modify_edge(he, ctr_opp(cv));
    }
  }

  //fix the directions of the curves (given correct marked face)
  // it should be called mostly after  symmetric_difference.
  void _fix_curves_direction(Aos_2& arr)
  {
    Compare_endpoints_xy_2 cmp_endpoints =
      arr.geometry_traits()->compare_endpoints_xy_2_object();
    Construct_opposite_2 ctr_opp = 
      arr.geometry_traits()->construct_opposite_2_object();

    for (Edge_iterator eit = arr.edges_begin();
         eit != arr.edges_end();
         ++eit)
    {
      Halfedge_handle            he = eit;
      const X_monotone_curve_2&  cv = he->curve();
      const bool                 is_cont = he->face()->contained();
      const Comparison_result    he_res = 
        ((Arr_halfedge_direction)he->direction() == ARR_LEFT_TO_RIGHT) ?
        SMALLER : LARGER;
      const bool                 has_same_dir = (cmp_endpoints(cv) == he_res);
      
      if ((is_cont && !has_same_dir) || (!is_cont && has_same_dir)) {
        arr.modify_edge(he, ctr_opp(cv));
      }
    }
  }

  void _build_sorted_vertices_vectors (std::vector<Arr_entry>& arr_vec)
  {
    Less_vertex_handle    comp (m_traits->compare_xy_2_object());
    Aos_2                 *p_arr;
    Vertex_iterator       vit;
    const std::size_t     n = arr_vec.size();
    std::size_t           i, j;
    
    for (i = 0; i < n; i++)
    {
      // Allocate a vector of handles to all vertices in the current
      // arrangement.
      p_arr = arr_vec[i].first;
      arr_vec[i].second = new std::vector<Vertex_handle>;
      arr_vec[i].second->resize (p_arr->number_of_vertices());
      
      for (j = 0, vit = p_arr->vertices_begin();
           vit != p_arr->vertices_end();
           j++, ++vit)
      {
        (*(arr_vec[i].second))[j] = vit;
      }
      
      // Sort the vector.
      std::sort (arr_vec[i].second->begin(), arr_vec[i].second->end(), comp);
    }
  }
  
  template <class Merge>
  void _divide_and_conquer (unsigned int lower, unsigned int upper,
                            std::vector<Arr_entry>& arr_vec,
                            unsigned int k, Merge merge_func)
  {
    if ((upper - lower) < k)
    {
      merge_func(lower, upper, 1, arr_vec);
      return;
    }
    
    unsigned int sub_size = ((upper - lower + 1) / k);
    unsigned int i = 0;
    unsigned int curr_lower = lower;
    
    for (; i<k-1; ++i, curr_lower += sub_size )
    {
      _divide_and_conquer(curr_lower, curr_lower + sub_size-1, arr_vec, k,
                          merge_func);
    }
    _divide_and_conquer (curr_lower, upper,arr_vec, k, merge_func);
    merge_func (lower, curr_lower, sub_size ,arr_vec);
    
    return;
  }
  
  // mark all faces as non-visited
  void _reset_faces() const
  {
    _reset_faces(m_arr);
  }
  
  void _reset_faces(Aos_2* arr) const
  {
    Face_const_iterator fit = arr->faces_begin();
    for ( ; fit != arr->faces_end(); ++fit)
    {
      fit->set_visited(false);
    }
  }


  void _insert(const Polygon_2& pgn, Aos_2& arr);
  
  // The function below is public because
  // are_holes_and_boundary_pairwise_disjoint of Gps_polygon_validation is
  // using it.
  // I have tried to define it as friend function, but with no success
  // (probably did something wrong with templates and friend.) Besides,
  // it was like this before I touched it, so I did not have the energy.
public:  
  void _insert(const Polygon_with_holes_2& pgn, Aos_2& arr);
  
protected:
  template<typename PolygonIter>
  void _insert(PolygonIter p_begin, PolygonIter p_end, Polygon_2& pgn);
  
  template<typename PolygonIter>
  void _insert(PolygonIter p_begin, PolygonIter p_end,
               Polygon_with_holes_2& pgn);
  
  template <typename OutputIterator>
  void _construct_curves(const Polygon_2& pgn, OutputIterator oi);
  
  template <typename OutputIterator>
  void _construct_curves(const Polygon_with_holes_2& pgn, OutputIterator oi);
  
  
  bool _is_empty(const Polygon_2& pgn) const
  {
    const std::pair<Curve_const_iterator, Curve_const_iterator>& itr_pair = 
      m_traits->construct_curves_2_object()(pgn);
    return (itr_pair.first == itr_pair.second);
  }
  
  bool _is_empty(const Polygon_with_holes_2& ) const
  {
    return (false);
  }
  
  bool _is_plane(const Polygon_2& ) const
  {
    return (false);
  }
  
  bool _is_plane(const Polygon_with_holes_2& pgn) const
  {
    //typedef typename  Traits_2::Is_unbounded  Is_unbounded;
    bool unbounded = m_traits->construct_is_unbounded_object()(pgn);
    std::pair<GP_Holes_const_iterator, 
      GP_Holes_const_iterator> pair = 
      m_traits->construct_holes_object()(pgn);
    return (unbounded && (pair.first == pair.second));
    //used to return
    //  (pgn.is_unbounded() && (pgn.holes_begin() == pgn.holes_end()))
  }
  
  void _intersection(const Aos_2& arr)
  {
    Aos_2* res_arr = new Aos_2(m_traits);
    Gps_intersection_functor<Aos_2> func;
    overlay(*m_arr, arr, *res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = res_arr;
    remove_redundant_edges();
    //fix_curves_direction(); // not needed for intersection
    CGAL_assertion(is_valid());
  }
  
  void _intersection(const Aos_2& arr1, const Aos_2& arr2, Aos_2& res) 
  {
    Gps_intersection_functor<Aos_2> func;
    overlay(arr1, arr2, res, func);
    _remove_redundant_edges(&res);
    //_fix_curves_direction(res); // not needed for intersection
    CGAL_assertion(_is_valid(res));
  }
  
  template <class Polygon_>
  void _intersection(const Polygon_& pgn)
  {
    if (_is_empty(pgn))
      this->clear();
    if (_is_plane(pgn)) return;
    if (this->is_empty()) return;
    if (this->is_plane())
    {
      Aos_2* arr = new Aos_2(m_traits);
      _insert(pgn, *arr);
      delete (this->m_arr);
      this->m_arr = arr;
      return;
    }
    
    Aos_2 second_arr;
    _insert(pgn, second_arr);
    _intersection(second_arr);
  }
  
  void _intersection(const Self& other)
  {
    if (other.is_empty())
    {
      m_arr->clear();
      return;
    }
    if (other.is_plane()) return;
    if (this->is_empty()) return;
    if (this->is_plane())
    {
      *(this->m_arr) = *(other.m_arr);
      return;
    }
    
    _intersection(*(other.m_arr));
  }
  
  void _join(const Aos_2& arr)
  {
    Aos_2* res_arr = new Aos_2(m_traits);
    Gps_join_functor<Aos_2> func;
    overlay(*m_arr, arr, *res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = res_arr;
    remove_redundant_edges();
    //fix_curves_direction(); // not needed for join
    CGAL_assertion(is_valid());
  }
  
  void _join(const Aos_2& arr1, const Aos_2& arr2, Aos_2& res) 
  {
    Gps_join_functor<Aos_2> func;
    overlay(arr1, arr2, res, func);
    _remove_redundant_edges(&res);
    //_fix_curves_direction(res); // not needed for join
    CGAL_assertion(_is_valid(res));
  }
  
  template <class Polygon_>
  void _join(const Polygon_& pgn)
  {
    if (_is_empty(pgn)) return;
    if (_is_plane(pgn))
    {
      this->clear();
      
      // Even in an empty arrangement there can be several faces
      // (because of the topology traits).
      for (Face_iterator fit = this->m_arr->faces_begin();
           fit != this->m_arr->faces_end(); ++fit)
        fit->set_contained(true);
      return;
    }
    if (this->is_empty())
    {
      Aos_2* arr = new Aos_2(m_traits);
      _insert(pgn, *arr);
      delete (this->m_arr);
      this->m_arr = arr;
      return;
    }
    if (this->is_plane()) return;
    
    Aos_2 second_arr;
    _insert(pgn, second_arr);
    _join(second_arr);
  }
  
  
  void _join(const Self& other)
  {
    if (other.is_empty()) return;
    if (other.is_plane())
    {
      this->clear();

      // Even in an empty arrangement there can be several faces
      // (because of the topology traits).
      for (Face_iterator fit = this->m_arr->faces_begin();
           fit != this->m_arr->faces_end(); ++fit)
        fit->set_contained(true);
      return;
    }
    if (this->is_empty())
    {
      *(this->m_arr) = *(other.m_arr);
      return;
    }
    if (this->is_plane()) return;
    _join(*(other.m_arr));
  }
  
  void _difference(const Aos_2& arr)
  {
    Aos_2* res_arr = new Aos_2(m_traits);
    Gps_difference_functor<Aos_2> func;
    overlay(*m_arr, arr, *res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = res_arr;
    remove_redundant_edges();
    fix_curves_direction();
    CGAL_assertion(is_valid());
  }
  
  void _difference(const Aos_2& arr1, const Aos_2& arr2, Aos_2& res) 
  {
    Gps_difference_functor<Aos_2> func;
    overlay(arr1, arr2, res, func);
    _remove_redundant_edges(&res);
    _fix_curves_direction(res);
    CGAL_assertion(_is_valid(res));
  }
  
  template <class Polygon_>
  void _difference(const Polygon_& pgn)
  {
    if (_is_empty(pgn)) return;
    if (_is_plane(pgn))
    {
      this->clear();
      return;
    }
    if (this->is_empty()) return;    
    if (this->is_plane())
    {
      Aos_2* arr = new Aos_2(m_traits);
      _insert(pgn, *arr);
      delete (this->m_arr);
      this->m_arr = arr;
      this->complement();
      return;
    }
    
    Aos_2 second_arr;
    _insert(pgn, second_arr);
    _difference(second_arr);
  }
  
  
  void _difference(const Self& other)
  {
    if (other.is_empty()) return;
    if (other.is_plane())
    {
      this->clear();
      return;
    }
    if (this->is_empty()) return;
    if (this->is_plane())
    {
      *(this->m_arr) = *(other.m_arr);
      this->complement();
      return;
    }
    
    _difference(*(other.m_arr));
  }
  
  void _symmetric_difference(const Aos_2& arr)
  {
    Aos_2* res_arr = new Aos_2(m_traits);
    Gps_sym_diff_functor<Aos_2> func;
    overlay(*m_arr, arr, *res_arr, func);
    delete m_arr; // delete the previous arrangement
    
    m_arr = res_arr;
    remove_redundant_edges();
    fix_curves_direction();
    CGAL_assertion(is_valid());
  }
  
  void _symmetric_difference(const Aos_2& arr1, const Aos_2& arr2, Aos_2& res) 
  {
    Gps_sym_diff_functor<Aos_2> func;
    overlay(arr1, arr2, res, func);
    _remove_redundant_edges(&res);
    _fix_curves_direction(res);
    CGAL_assertion(_is_valid(res));
  }
  
  template <class Polygon_>
  void _symmetric_difference(const Polygon_& pgn)
  {
    if (_is_empty(pgn)) return;
    
    if (_is_plane(pgn))
    {
      this->complement();
      return;
    }
    if (this->is_empty())
    {
      Aos_2* arr = new Aos_2(m_traits);
      _insert(pgn, *arr);
      delete (this->m_arr);
      this->m_arr = arr;
      return;
    }
    
    if (this->is_plane())
    {
      Aos_2* arr = new Aos_2(m_traits);
      _insert(pgn, *arr);
      delete (this->m_arr);
      this->m_arr = arr;
      this->complement();
      return;
    }
    
    Aos_2 second_arr;
    _insert(pgn, second_arr);
    _symmetric_difference(second_arr);
  }
  
  
  void _symmetric_difference(const Self& other)
  {
    if (other.is_empty()) return;
    
    if (other.is_plane())
    {
      this->complement();
      return;
    }
    if (this->is_empty())
    {
      *(this->m_arr) = *(other.m_arr);
      return;
    }
    
    if (this->is_plane())
    {
      *(this->m_arr) = *(other.m_arr);
      this->complement();
      return;
    }
    
    _symmetric_difference(*(other.m_arr));
  }
};

#include <CGAL/Boolean_set_operations_2/Gps_on_surface_base_2_impl.h>

} //namespace CGAL

#endif // CGAL_GPS_ON_SURFACE_BASE_2_H
