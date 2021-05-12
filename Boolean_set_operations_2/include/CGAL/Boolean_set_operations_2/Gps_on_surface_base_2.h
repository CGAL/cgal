// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ophir Setter    <ophir.setter@cs.tau.ac.il>
//                 Guy Zucker <guyzucke@post.tau.ac.il>


#ifndef CGAL_GPS_ON_SURFACE_BASE_2_H
#define CGAL_GPS_ON_SURFACE_BASE_2_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

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
#include <CGAL/Union_find.h>


/*!
  \file   Gps_on_surface_base_2.h
  \brief  A class that allows Boolean set operations.
  This class is the base class for General_polygon_set_on_surface_2 and
  receives extra template parameter which allows different validation
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
    inline static void is_valid(const Polygon&, const Traits&) {}
  };
}

//! General_polygon_set_on_surface_2
/*! This class is the base class for General_polygon_set_on_surface_2 and
    receives extra template parameter which allows different validation
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
  const Traits_2* m_traits;
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
  Gps_on_surface_base_2(const Traits_2& tr) : m_traits(&tr),
                                        m_traits_adaptor(*m_traits),
                                        m_traits_owner(false),
                                        m_arr(new Aos_2(m_traits))
  {}

  // Copy constructor
  Gps_on_surface_base_2(const Self& ps) :
    m_traits(new Traits_2(*(ps.m_traits))),
    m_traits_adaptor(*m_traits),
    m_traits_owner(true),
    m_arr(new Aos_2(*(ps.m_arr)))
  {}

  // Asignment operator
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

  // Constructor
  explicit Gps_on_surface_base_2(const Polygon_2& pgn) :
    m_traits(new Traits_2()),
    m_traits_adaptor(*m_traits),
    m_traits_owner(true),
    m_arr(new Aos_2(m_traits))
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _insert(pgn, *m_arr);
  }

  // Constructor
  explicit Gps_on_surface_base_2(const Polygon_2& pgn, const Traits_2& tr) :
    m_traits(&tr),
    m_traits_adaptor(*m_traits),
    m_traits_owner(false),
    m_arr(new Aos_2(m_traits))
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    _insert(pgn, *m_arr);
  }

  // Constructor
  explicit Gps_on_surface_base_2(const Polygon_with_holes_2& pgn_with_holes) :
    m_traits(new Traits_2()),
    m_traits_adaptor(*m_traits),
    m_traits_owner(true),
    m_arr(new Aos_2(m_traits))
  {
    ValidationPolicy::is_valid(pgn_with_holes,*m_traits);
    _insert(pgn_with_holes, *m_arr);
  }

  // Constructor
  explicit Gps_on_surface_base_2(const Polygon_with_holes_2& pgn_with_holes,
                                 const Traits_2& tr) :
    m_traits(&tr),
    m_traits_adaptor(*m_traits),
    m_traits_owner(false),
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
  bool do_intersect(const Polygon_2& pgn) const
  {
    ValidationPolicy::is_valid(pgn, *m_traits);
    Self other(pgn, *m_traits);
    return (do_intersect(other));
  }

  // test for intersection with a polygon with holes
  bool do_intersect(const Polygon_with_holes_2& pgn_with_holes) const
  {
    ValidationPolicy::is_valid(pgn_with_holes, *m_traits);
    Self other(pgn_with_holes, *m_traits);
    return (do_intersect(other));
  }

  //test for intersection with another Gps_on_surface_base_2 object
  bool do_intersect(const Self& other) const
  {
    if (this->is_empty() || other.is_empty()) return false;
    if (this->is_plane() || other.is_plane()) return true;
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

  // Traits_2& traits()
  // {
  //   return *m_traits;
  // }

  const Traits_2& traits() const { return *m_traits; }

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
                                Polygon_2 & pgn, const Traits_2* tr);

  bool is_hole_of_face(Face_const_handle f, Halfedge_const_handle he) const;

  Ccb_halfedge_const_circulator
  get_boundary_of_polygon(Face_const_iterator f) const;

  void remove_redundant_edges()
  {
    this->_remove_redundant_edges(m_arr);
  }

protected:

  bool is_redundant(Halfedge_handle he)
  {
    return he->face()->contained() == he->twin()->face()->contained();
  }

  typename Aos_2::Dcel::Halfedge*
  _halfedge(Halfedge_handle h)
  {
    return &(*h);
  }

  void set_flag_of_halfedges_of_final_argt(Halfedge_handle h, int flag)
  {
    Halfedge_handle start=h;
    do{
      h->set_flag(flag);
      h=h->next();
      while (is_redundant(h))
        h=h->twin()->next();
    } while(start!=h);
  }

  void _remove_redundant_edges(Aos_2* arr)
  {
    // const integer for handling the status of halfedges
    // during the flooding algorithm to tag halfedges as
    // on an inner or outer ccb in the final arrangement
    static const int ON_INNER_CCB=0;
    static const int ON_OUTER_CCB=1;
    static const int NOT_VISITED=-1;
    static const int NEW_CCB_ASSIGNED=2;

    // Consider the faces incident to a redundant edge and use a union-find
    // algorithm to group faces in set that will be merged by the removal
    // of redundant edges. Then only the master of the set will be kept.
    // Here we also collect edges that needs to be removed.
    typedef Union_find<typename Aos_2::Dcel::Face*> UF_faces;
    std::vector<typename UF_faces::handle> face_handles;
    UF_faces uf_faces;
    std::vector< typename Aos_2::Dcel::Halfedge* > edges_to_remove;
    bool all_edges_are_redundant=true;

    for (Edge_iterator itr = arr->edges_begin(); itr != arr->edges_end(); ++itr)
    {
      Halfedge_handle he = itr;
      he->set_flag(NOT_VISITED);
      he->twin()->set_flag(NOT_VISITED);

      // put in the same set faces that will be merged when removing redundant edges
      if ( is_redundant(he) )
      {
        typename Aos_2::Dcel::Face* f1=&(*he->face()),
                                  * f2=&(*he->twin()->face());
        if (f1->id_not_set()){
          f1->set_id(face_handles.size());
          face_handles.push_back( uf_faces.make_set( f1 ) );
        }
        if (f2->id_not_set()){
          f2->set_id(face_handles.size());
          face_handles.push_back( uf_faces.make_set( f2 ) );
        }

        uf_faces.unify_sets(face_handles[f1->id()], face_handles[f2->id()]);
        edges_to_remove.push_back( _halfedge(he) );
      }
      else
        all_edges_are_redundant=false;
    }
    // the code in this function assumes there is only one unbounded face
    // (in the if below and in the part to keep the unbounded face even if
    //  not the master of its set)
    CGAL_assertion(std::distance(arr->unbounded_faces_begin(),
                                 arr->unbounded_faces_end()) == 1);

    if (all_edges_are_redundant){
      bool is_contained=arr->unbounded_faces_begin()->contained();
      arr->clear();
      arr->unbounded_faces_begin()->set_contained(is_contained);
      return;
    }

    // nothing needs to be done
    if (edges_to_remove.empty() ) return;

  // Start tagging ccbs
    // For all halfedge that is part of a face that will be subject to a merge
    // due to the removal of redundant edges, we now flag whether the halfedge
    // will be part of an outer ccb or an inner ccb in the final arrangement
    // (i.e. the arrangement after the removal of redundant edges).
    // We use the following invariant:
    //  - a halfedge on an inner ccb will be on an inner ccb in the final argt
    //  - the twin of a halfedge on an inner ccb will be on an outer ccb in the final argt
    //  - there is only one outer ccb per bounded face
    //  - the unbounded face has no ccb

    // first collect all non-redundant halfedges
    std::vector<Halfedge_handle> halfedges_that_was_on_an_outer_ccb;
    // bitset indicating if the outer_ccb of a face was already set
    std::vector<bool> face_outer_ccb_set(face_handles.size(),false);
    for (Halfedge_iterator itr = arr->halfedges_begin(); itr != arr->halfedges_end(); ++itr)
    {
      Halfedge_handle h = itr;
      if (is_redundant(h))
      {
        // mark redundant edges as we will reuse ccb, thus breaking the function is_redundant()
        // needed for "update halfedge ccb pointers"
        h->set_flag(NEW_CCB_ASSIGNED);
        h->twin()->set_flag(NEW_CCB_ASSIGNED);
      }
      else{
        // tag halfedges of modified faces that are on an inner ccb
        // or twin of a halfedge on an inner ccb.
        if (h->flag()!=NOT_VISITED) continue;
        if(h->is_on_inner_ccb())
        {
          //visit inner ccb of h in the final arrangement
          if (!h->face()->id_not_set())
            set_flag_of_halfedges_of_final_argt(h,ON_INNER_CCB);
          CGAL_assertion(h->twin()->is_on_outer_ccb());
          if ( h->twin()->flag()!=NOT_VISITED ||
               h->twin()->face()->id_not_set()) continue;
          //visit outer ccb of h in the final arrangement
          set_flag_of_halfedges_of_final_argt(h->twin(),ON_OUTER_CCB);
          std::size_t master_id=
            (*uf_faces.find(face_handles[h->twin()->face()->id()]))->id();
          face_outer_ccb_set[master_id]=true;
        }
        else{
          if (!h->face()->id_not_set())
            halfedges_that_was_on_an_outer_ccb.push_back(h);
        }
      }
    }

    bool something_was_updated;
    // iterative step to propagate changes layer by layer
    do{
      something_was_updated=false;
      // update the bitset using the bit value of the master set
      // and also set the bit of the unbounded cc to 1 (as it has no unbounded ccb)
      for(typename UF_faces::iterator it=uf_faces.begin(),
                                      it_end=uf_faces.end(); it!=it_end; ++it)
      {
        if (face_outer_ccb_set[(*it)->id()]) continue;
        typename UF_faces::handle master=uf_faces.find(it);
        //remove faces that are not the master of their set (but the unbounded face)
        if ((*it)->is_unbounded())
          face_outer_ccb_set[(*master)->id()]=true;
        if ( master!=it)
          face_outer_ccb_set[(*it)->id()]=face_outer_ccb_set[(*master)->id()];
      }

      // update halfedge flag according to the flag of the twin halfedge
      // or if the outer ccb of the cc was set
      for(Halfedge_handle h : halfedges_that_was_on_an_outer_ccb)
      {
        if (h->flag()!=NOT_VISITED) continue;
        std::size_t face_master_id=(*uf_faces.find(face_handles[h->face()->id()]))->id();
        if (h->twin()->flag()==ON_INNER_CCB){
          set_flag_of_halfedges_of_final_argt(h,ON_OUTER_CCB);
          face_outer_ccb_set[face_master_id]=true;
          something_was_updated=true;
        }
        else
        {
          if (face_outer_ccb_set[face_master_id]){
            set_flag_of_halfedges_of_final_argt(h,ON_INNER_CCB);
            something_was_updated=true;
          }
        }
      }
    }
    while(something_was_updated);
    // last loop, if some tags are not set it means that they are the only ccb
    // of the face and that they have to be the outer ccb
    for(Halfedge_handle h : halfedges_that_was_on_an_outer_ccb)
    {
      if (h->flag()!=NOT_VISITED) continue;
      std::size_t face_master_id=(*uf_faces.find(face_handles[h->face()->id()]))->id();
      set_flag_of_halfedges_of_final_argt(h,ON_OUTER_CCB);
      face_outer_ccb_set[face_master_id]=true;
    }
    // at this position there might be some bits in face_outer_ccb_set not set
    // but they are corresponding to the unbounded face
  // End tagging ccbs

    // update the next/prev relationship around vertices kept incident
    // to at least one edge to remove. We link non redundant halfedges together.
    //We also collect vertices to remove at the same time.
    std::vector< typename Aos_2::Dcel::Vertex* > vertices_to_remove;
    for(Vertex_iterator vi=arr->vertices_begin(), vi_end=arr->vertices_end(); vi!=vi_end; ++vi)
    {
      typename Aos_2::Dcel::Vertex* v_ptr=&(*vi);
      Halfedge_handle h_start=vi->incident_halfedges(), h=h_start;

      std::vector<Halfedge_handle> non_redundant_edges;
      bool found_no_redundant=true;
      do{
        if( !is_redundant(h) )
          non_redundant_edges.push_back(h);
        else{
          found_no_redundant=false;
        }
        h=h->next()->twin();
      }while(h!=h_start);

      // if only redundant edges are incident to the vertex, then the
      // vertex will be removed and nothing needs to be done.
      if (non_redundant_edges.empty()){
        vertices_to_remove.push_back(v_ptr);
        continue;
      }
      //if the vertex neighbor is already correct, then continue
      if (found_no_redundant) continue;

      std::size_t nb_edges=non_redundant_edges.size();
      CGAL_assertion( nb_edges >= 2);

      non_redundant_edges.push_back(non_redundant_edges.front());

      //update vertex halfedge
      v_ptr->set_halfedge(_halfedge(non_redundant_edges.back()));
      for (std::size_t i=0; i<nb_edges; ++i)
      {
        Halfedge_handle h1 = non_redundant_edges[i], h2=non_redundant_edges[i+1];
        if ( h1->next()->twin()!=h2)
          _halfedge(h1)->set_next(_halfedge(h2->twin()));
      }
    }

    //collect faces to remove and update unbounded face flag
    std::vector< typename Aos_2::Dcel::Face*> faces_to_remove;
    std::vector< typename Aos_2::Dcel::Outer_ccb* > outer_ccbs_to_remove;
    std::vector< typename Aos_2::Dcel::Inner_ccb* > inner_ccbs_to_remove;
    for(typename UF_faces::iterator it=uf_faces.begin(),
                                    it_end=uf_faces.end(); it!=it_end; ++it)
    {
      typename UF_faces::handle master=uf_faces.find(it);
      //remove faces that are not the master of their set (but the unbounded face)
      if ( master!=it)
      {
        // force to keep the unbounded face
        if ((*it)->is_unbounded())
        {
          face_handles[(*master)->id()]=it;
          faces_to_remove.push_back(*master);
        }
        else
          faces_to_remove.push_back(*it);
      }

      //collect for reuse/removal all inner and outer ccbs
      for(void* ptr : (*it)->_outer_ccbs())
        outer_ccbs_to_remove.push_back( static_cast<typename Aos_2::Dcel::Halfedge*>(ptr)->outer_ccb() );
      for(void* ptr : (*it)->_inner_ccbs())
        inner_ccbs_to_remove.push_back( static_cast<typename Aos_2::Dcel::Halfedge*>(ptr)->inner_ccb() );
      (*it)->_outer_ccbs().clear();
      (*it)->_inner_ccbs().clear();
    }

    // accessor for  low-level arrangement fonctionalities
    CGAL::Arr_accessor<Aos_2> accessor(*arr);
    // the face field of outer and inner ccb are used in the loop to access the old face an halfedge
    // used to contribute to. These two vectors are used to delay the association to the new face to
    // avoid overwriting a field that is still needed
    typedef std::pair<typename Aos_2::Dcel::Outer_ccb*, typename Aos_2::Dcel::Face*> Outer_ccb_and_face;
    typedef std::pair<typename Aos_2::Dcel::Inner_ccb*, typename Aos_2::Dcel::Face*> Inner_ccb_and_face;
    std::vector<Outer_ccb_and_face> outer_ccb_and_new_face_pairs;
    std::vector<Inner_ccb_and_face> inner_ccb_and_new_face_pairs;
    // update halfedge ccb pointers
    for (Halfedge_iterator itr = arr->halfedges_begin(); itr != arr->halfedges_end(); ++itr)
    {
      Halfedge_handle h = itr;
      CGAL_assertion(h->face() != Face_handle());
      if (h->face()->id_not_set()) continue;
      CGAL_assertion(h->flag()!=NOT_VISITED);

      // either a redundant edge or an edge of an already handled ccb
      if ( h->flag()==NEW_CCB_ASSIGNED ) continue;

      CGAL_assertion( h->flag()==ON_INNER_CCB || h->flag()==ON_OUTER_CCB );

      typename Aos_2::Dcel::Face* f=&(*h->face());

      if (!f->id_not_set())
      {
        // we use the master of the set as face, but we force to keep the unbounded face,
        // thus this hack
        f = *(face_handles[
                (*uf_faces.find(face_handles[f->id()]))->id()
              ]);
        if (h->flag()==ON_INNER_CCB)
        {
          bool reuse_inner_ccb = !inner_ccbs_to_remove.empty();
          typename Aos_2::Dcel::Inner_ccb* inner_ccb = !reuse_inner_ccb?
            accessor.new_inner_ccb():inner_ccbs_to_remove.back();
          if ( reuse_inner_ccb ) inner_ccbs_to_remove.pop_back();

          Halfedge_handle hstart=h;
          do{
            _halfedge(h)->set_inner_ccb(inner_ccb);
            h->set_flag(NEW_CCB_ASSIGNED);
            h=h->next();
          }while(hstart!=h);
          f->add_inner_ccb(inner_ccb,_halfedge(h));
          inner_ccb->set_halfedge(_halfedge(h));
          if (!reuse_inner_ccb)
            inner_ccb->set_face(f);
          else
            inner_ccb_and_new_face_pairs.push_back( std::make_pair(inner_ccb, f) );
        }
        else{
          // create a new outer ccb if none is available
          typename Aos_2::Dcel::Outer_ccb* outer_ccb;
          if (!outer_ccbs_to_remove.empty())
          {
            outer_ccb = outer_ccbs_to_remove.back();
            outer_ccbs_to_remove.pop_back();
          }
          else{
            outer_ccb = accessor.new_outer_ccb();
            outer_ccb->set_face(f);
          }
          Halfedge_handle hstart=h;
          do{
            _halfedge(h)->set_outer_ccb(outer_ccb);
            h->set_flag(NEW_CCB_ASSIGNED);
            h=h->next();
          }while(hstart!=h);
          f->add_outer_ccb(outer_ccb,_halfedge(h));
          outer_ccb->set_halfedge(_halfedge(h));
          outer_ccb_and_new_face_pairs.push_back( std::make_pair(outer_ccb, f) );
        }
      }
    }

    // now set the new face for all ccbs
    for(Outer_ccb_and_face& ccb_and_face : outer_ccb_and_new_face_pairs)
      ccb_and_face.first->set_face(ccb_and_face.second);
    for(Inner_ccb_and_face& ccb_and_face : inner_ccb_and_new_face_pairs)
      ccb_and_face.first->set_face(ccb_and_face.second);

    //remove no longer used edges, vertices and faces
    accessor.delete_vertices( vertices_to_remove );
    accessor.delete_edges( edges_to_remove );
    accessor.delete_faces( faces_to_remove );
    accessor.delete_outer_ccbs( outer_ccbs_to_remove );
    accessor.delete_inner_ccbs( inner_ccbs_to_remove );

    for (typename Aos_2::Face_iterator fit=arr->faces_begin(),
                                       end=arr->faces_end();
                                       fit!=end; ++fit)
    {
      fit->reset_id(); // reset the id that will be no longer used for this face
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

} //namespace CGAL

#include <CGAL/Boolean_set_operations_2/Gps_on_surface_base_2_impl.h>

#include <CGAL/enable_warnings.h>

#endif // CGAL_GPS_ON_SURFACE_BASE_2_H
