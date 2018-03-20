// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
//                 (based on old version by Oren Nechushtan and Iddo Hanniel)

#ifndef CGAL_ARR_TRAPEZOID_RIC_POINT_LOCATION_H
#define CGAL_ARR_TRAPEZOID_RIC_POINT_LOCATION_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * Definition of the Arr_trapezoid_ric_point_location<Arrangement> template.
 */

#include <CGAL/Arr_point_location_result.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_point_location/Trapezoidal_decomposition_2.h>
#include <CGAL/Arr_point_location/Td_traits.h>
#include <CGAL/Arr_observer.h>

namespace CGAL {

/*! \class
 * A class that answers point-location and queries
 * on a planar arrangement using the trapezoid_ric algorithm.
 * The Arrangement parameter corresponds to an arrangement instantiation.
 */
template <typename Arrangement_>
class Arr_trapezoid_ric_point_location : public Arr_observer <Arrangement_> {
public:
  //type of arrangement on surface
  typedef Arrangement_                          Arrangement_on_surface_2;

  //type of geometry traits
  typedef typename Arrangement_on_surface_2::Geometry_traits_2		
                                                Geometry_traits_2;
  //type of traits adaptor
  typedef typename Arrangement_on_surface_2::Traits_adaptor_2		
                                                Traits_adaptor_2;
  //type of vertex handle
  typedef typename Arrangement_on_surface_2::Vertex_handle
                                                Vertex_handle;
  //type of vertex const handle
  typedef typename Arrangement_on_surface_2::Vertex_const_handle	
                                                Vertex_const_handle;
  //type of halfedge handle
  typedef typename Arrangement_on_surface_2::Halfedge_handle	
                                                Halfedge_handle;
  //type of halfedge const handle
  typedef typename Arrangement_on_surface_2::Halfedge_const_handle	
                                                Halfedge_const_handle;
  //type of face const handle
  typedef typename Arrangement_on_surface_2::Face_const_handle		
                                                Face_const_handle;
  //type of edge const iterator
  typedef typename Arrangement_on_surface_2::Edge_const_iterator	
                                                Edge_const_iterator;
  //type of isolated vertex const iterator
  typedef typename Arrangement_on_surface_2::Isolated_vertex_const_iterator
                                       Isolated_vertex_const_iterator;
  //type of point
  typedef typename Geometry_traits_2::Point_2             Point_2;

  //type of x-monotone curve
  typedef typename Geometry_traits_2::X_monotone_curve_2  
                                                X_monotone_curve_2;
  
  //type of trapezoidal decomposition traits class
  typedef CGAL::Td_traits<Traits_adaptor_2, Arrangement_on_surface_2> 
                                                Td_traits;
  //type of trapezoidal decomposition class
  typedef Trapezoidal_decomposition_2<Td_traits>    
                                                Trapezoidal_decomposition;
  
  //!types of Td_map_item-s
  typedef typename Trapezoidal_decomposition::Td_map_item       
                                                Td_map_item;
  typedef typename Trapezoidal_decomposition::Td_active_vertex  
                                                Td_active_vertex;
  typedef typename Trapezoidal_decomposition::Td_active_fictitious_vertex  
                                                Td_active_fictitious_vertex;
  typedef typename Trapezoidal_decomposition::Td_active_edge  
                                                Td_active_edge;
  typedef typename Trapezoidal_decomposition::Td_active_trapezoid  
                                                Td_active_trapezoid;
  //!type of side tags
  typedef typename Traits_adaptor_2::Left_side_category   
                                                Left_side_category;
  typedef typename Traits_adaptor_2::Bottom_side_category 
                                                Bottom_side_category;
  typedef typename Traits_adaptor_2::Top_side_category    
                                                Top_side_category;
  typedef typename Traits_adaptor_2::Right_side_category  
                                                Right_side_category;
 
protected:
  typedef Arr_point_location_result<Arrangement_on_surface_2>        Result;
  typedef typename Result::Type                           Result_type;

public:
  // Support cpp11::result_of
  typedef Result_type                                     result_type;

protected:
  //type of trapezoidal decomposition class
  typedef Trapezoidal_decomposition             TD;

  typedef typename Arr_are_all_sides_oblivious_tag< 
                     Left_side_category, Bottom_side_category, 
                     Top_side_category, Right_side_category >::result
    Are_all_sides_oblivious_category;

  // Data members:
  const Traits_adaptor_2 *m_traits;  // Its associated traits object.
  TD  td;       // instance of trapezoidal decomposition
  bool m_with_guarantees;
  //for the notification functions
  X_monotone_curve_2  m_cv_before_split; 
  Halfedge_handle     m_he_after_merge;
  //X_monotone_curve_2        m_cv_before_merge1;
  //X_monotone_curve_2        m_cv_before_merge2;

  template<typename T>
  Result_type make_result(T t) const { return Result::make_result(t); }
  inline Result_type default_result() const { return Result::default_result(); }

public:
  /*! Default constructor. */
  Arr_trapezoid_ric_point_location (bool with_guarantees = true, 
                           double depth_thrs = CGAL_TD_DEFAULT_DEPTH_THRESHOLD, 
                           double size_thrs = CGAL_TD_DEFAULT_SIZE_THRESHOLD) 
    : m_traits (NULL), m_with_guarantees(with_guarantees)
  {
    td.set_with_guarantees(with_guarantees);
    td.depth_threshold(depth_thrs);
    td.size_threshold(size_thrs);
  }

  /*! Constructor given an arrangement. */
  Arr_trapezoid_ric_point_location (const Arrangement_on_surface_2& arr, 
                           bool with_guarantees = true, 
                           double depth_thrs = CGAL_TD_DEFAULT_DEPTH_THRESHOLD, 
                           double size_thrs = CGAL_TD_DEFAULT_SIZE_THRESHOLD) :
    Arr_observer<Arrangement_on_surface_2> 
              (const_cast<Arrangement_on_surface_2 &>(arr)),
    m_with_guarantees(with_guarantees)
  {
    m_traits = static_cast<const Traits_adaptor_2*> (arr.geometry_traits());
    td.set_with_guarantees(with_guarantees);
    td.init_arrangement_and_traits(&arr);
    td.depth_threshold(depth_thrs);
    td.size_threshold(size_thrs);
    _construct_td();
  }

  /*! Destructor. */
  ~Arr_trapezoid_ric_point_location () { }

  /*! defines whether the underlying search structure guarantees logarithmic 
   *   query time and linear size */
  void with_guarantees (bool with_guarantees)
  {
    //if with_guarantees was changed from false to true - reconstruct 
    //  the search structure with guarantees
    td.set_with_guarantees(with_guarantees);
    if (with_guarantees && !m_with_guarantees)
    {
      td.clear();
      _construct_td();
    }
    m_with_guarantees = with_guarantees;
  }
  
  /*! returns the depth of the underlying search structure 
   *    (the longest path in the DAG)
   */
  unsigned long depth() //longest_dag_path()
  {
    return td.largest_leaf_depth() + 1;
  }

  /*! returns the longest query path in the underlying search structure */
  unsigned long longest_query_path_length()
  {
    return td.longest_query_path_length();
  }

#ifdef CGAL_TD_DEBUG
  //void  locate_and_print (std::ostream& out, const Point_2& p) const
  //{
  //  td.locate_and_print(out, p);
  //}

  void print_dag(std::ostream& out) const
  {
    td.print_dag(out);
  }
#endif

  /*!
   * Locate the arrangement feature containing the given point.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type locate(const Point_2& p) const;

  /*!
   * Locate the arrangement feature which a upward vertical ray emanating from
   * the given point hits.
   * \param p The query point.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either an empty object or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type ray_shoot_up(const Point_2& p) const
  { return (_vertical_ray_shoot(p, true)); }

  /*!
   * Locate the arrangement feature which a downward vertical ray emanating
   * from the given point hits.
   * \param p The query point.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either an empty object or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type ray_shoot_down(const Point_2& p) const
  { return (_vertical_ray_shoot(p, false)); }

  /// \name Notification functions, inherited and overloaded from the
  //        base observer.
  //@{

  virtual void before_assign (const Arrangement_on_surface_2& arr)
  {
    td.clear();
    m_traits = static_cast<const Traits_adaptor_2*> (arr.geometry_traits());
    td.init_arrangement_and_traits(&arr, false);
  }

  virtual void after_assign ()
  { 
    _construct_td();
  }

  virtual void before_clear ()
  {
    td.clear();
  }

  virtual void after_clear ()
  {
    _construct_td();
  }

  virtual void before_attach (const Arrangement_on_surface_2& arr)
  {
    td.clear();
    m_traits = static_cast<const Traits_adaptor_2*> (arr.geometry_traits());
    td.init_arrangement_and_traits(&arr);
  }

  virtual void after_attach ()
  {
    _construct_td();
  }

  virtual void before_detach ()
  {
    td.clear();
  }

  virtual void after_create_edge (Halfedge_handle e)
  {
    td.insert(e);
  }

  //TODO IDIT OREN: what can be done in order to avoid the need 
  //to save the original curve is to find the common endpoint of the 
  //two new halfedges, locate it in the trapezoid in order to find the 
  //curve it lies on, which is the curve that was split, and then remove 
  //this curve.
  virtual void before_split_edge (Halfedge_handle e,
				                          Vertex_handle /* v */,
                                  const X_monotone_curve_2&  /* cv1 */ ,
                                  const X_monotone_curve_2&  /* cv2 */ )
  {
  
    ////MICHAL: commented due to inefficient depth update, remove and insert instead
    ////save the curve for the "after" function.
    //m_cv_before_split = e->curve();
    //td.before_split_edge(m_cv_before_split, cv1, cv2); 
    
    td.remove(e);
  }

  virtual void after_split_edge (Halfedge_handle e1,
                                 Halfedge_handle e2)
  {
    //MICHAL: commented due to inefficient depth update, remove and insert instead
    //td.split_edge(m_cv_before_split,e1,e2);
    
    td.insert(e1);
    td.insert(e2);
  }

  virtual void before_merge_edge (Halfedge_handle e1,
                                  Halfedge_handle e2,
                                  const X_monotone_curve_2& cv)
  {
    //save the halfedge handle for the "after" function.
    m_he_after_merge = e1;
    td.merge_edge (e1, e2, cv);
  }

  virtual void after_merge_edge (Halfedge_handle e)
  {
    td.after_merge_edge(e, m_he_after_merge);
  }

  virtual void before_remove_edge (Halfedge_handle e)
  {
    //called before combinatoric deletion
    td.remove(e);
  }
  //@}

public:

//#ifdef CGAL_TD_DEBUG
//  void debug()
//  {
//    td.debug();
//  }
//#endif

protected:

  /*! Construct the trapezoidal decomposition. */
  void _construct_td ()
  {
    td.clear();
 
    std::vector<Halfedge_const_handle> he_container; 
    Edge_const_iterator eit;
    Halfedge_const_handle he_cst;
    Arrangement_on_surface_2 *arr = this->arrangement();
    //collect the arrangement halfedges
    for (eit = arr->edges_begin(); eit != arr->edges_end(); ++eit)
    {
      he_cst = eit;
      he_container.push_back(he_cst);
    }
    //container insertion
    td.insert(he_container.begin(), he_container.end()); 
  }

  /*! gets the unbounded face that contains the point when the trapezoid is
   * unbounded
   * \param tr The unbounded trapezoid whose face we should get
   * \param p  The query point.
   * \param Arr_all_sides_oblivious_tag
   * \return A Face_const_handle representing the arrangement unbounded face in
   *         which the point p lies
   */ 
  Face_const_handle _get_unbounded_face (const Td_map_item& tr,
                                         const Point_2& p, 
                                         Arr_all_sides_oblivious_tag) const;

  /*! gets the unbounded face that contains the point when the trapezoid is
   * unbounded
   * \param tr The unbounded trapezoid whose face we should get
   * \param p  The query point.
   * \param Arr_not_all_sides_oblivious_tag
   * \return A Face_const_handle representing the arrangement unbounded face in which 
   *         the point p lies
   */ 
  Face_const_handle _get_unbounded_face (const Td_map_item& tr,
                                         const Point_2& p, 
                                         Arr_not_all_sides_oblivious_tag) const;

  /*!
   * Locate the arrangement feature which a vertical ray emanating from the
   * given point hits, considering isolated vertices.
   * \param p The query point.
   * \param shoot_up Indicates whether the ray is directed upward or downward.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either a Halfedge_const_handle,
   *         a Vertex_const_handle or an empty object.
   */
  result_type _vertical_ray_shoot(const Point_2& p, bool shoot_up) const;

  /*! In vertical ray shoot, when the closest halfedge is found
   * (or unbounded face)
   * we check the isolated vertices inside the face to check whether there
   * is an isolated vertex right above/below the query point.
   */ 
  result_type
  _check_isolated_for_vertical_ray_shoot
                             (Halfedge_const_handle halfedge_found, 
                              const Point_2& p, bool shoot_up,
                              const Td_map_item& tr) const;
};

} //namespace CGAL

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Arr_trapezoid_ric_pl_impl.h>

#include <CGAL/enable_warnings.h>

#endif
