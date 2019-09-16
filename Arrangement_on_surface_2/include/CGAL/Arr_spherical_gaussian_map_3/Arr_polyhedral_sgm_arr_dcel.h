// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s): Efi Fogel         <efif@post.tau.ac.il>
//            Naama mayer       <naamamay@post.tau.ac.il>

#ifndef CGAL_ARR_POLYHEDRAL_SGM_ARR_DCEL_H
#define CGAL_ARR_POLYHEDRAL_SGM_ARR_DCEL_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/basic.h>

namespace CGAL {

/*! Extend the arrangement vertex */
template <class Point_2>
class Arr_polyhedral_sgm_arr_vertex : public CGAL::Arr_vertex_base<Point_2> {
public:
  /*! Constructor */
  Arr_polyhedral_sgm_arr_vertex() {}
};

/*! Extend the arrangement halfedge */
template <class X_monotone_curve_2>
class Arr_polyhedral_sgm_arr_halfedge :
  public CGAL::Arr_halfedge_base<X_monotone_curve_2>
{
private:
  /*! A mask of the ids of the original arrangements that contributed the
   * halfedge while performing the minkowski sum.
   * \todo This should be made optional. It is relevant only if the polytope
   * is the result of a Minkowski sum operation, and it is needed only by the 
   * drawing routines.
   */
  unsigned int m_arr_mask;
  
public:
  /*! Constructor */
  Arr_polyhedral_sgm_arr_halfedge() : m_arr_mask(0x0) {}

  /*! Add a arrangement to the mask of the original arrangements in the
   * minkowski sum.
   * \param arr_id the id of the added arrangement
   */
  void add_arr(unsigned int id) { m_arr_mask |= 0x1 << id; }

  /*! Return true iff a given arrangement contributed this halfedge
   * while performing the minkowski sum
   */
  bool is_arr(unsigned int id) const { return m_arr_mask & (0x1 << id); }

  /*! Obtain the mask of the ids of the original arrangements that contributed
   * the halfedge while performing the minkowski sum
   */
  unsigned int arr_mask() const { return m_arr_mask; }

  /*! Set the arr of an edge with a value.
  * \param arr_id the id to set to.
  */
  void set_arr(unsigned int id) { m_arr_mask  = id; }
};

/*! Extend the arrangement face */
template <class Point_3>
class Arr_polyhedral_sgm_arr_face : public CGAL::Arr_face_base {
private:
  /*! The original point of the polyhedron */
  Point_3 m_point;

  /*! Indicates that the point has been set already */
  bool m_is_set;

public:
  /*! Constructor */
  Arr_polyhedral_sgm_arr_face() : m_is_set(false) { }
    
  /*! Set the 3D point of the original polyhedron */
  void set_point(const Point_3 & point)
  {
    m_point = point;
    m_is_set = true;
  }

  /*! Obtain the 3D point of the original polyhedron */
  const Point_3 & point() const { return m_point; }

  /*! \brief returns true iff the point has been set already */
  bool is_set() const { return m_is_set; }

  /*! \brief resets the flag  */
  void set_is_set(bool flag) { m_is_set = flag; }
};

/*! A new dcel builder with SGM features */
template <class Traits>
class Arr_polyhedral_sgm_arr_dcel :
  public CGAL::Arr_dcel_base<Arr_polyhedral_sgm_arr_vertex<typename Traits::Point_2>,
                             Arr_polyhedral_sgm_arr_halfedge<typename Traits::X_monotone_curve_2>,
                             Arr_polyhedral_sgm_arr_face<typename Traits::Point_3> >
{
public:
  /*! Constructor */
  Arr_polyhedral_sgm_arr_dcel() {}
};

} //namespace CGAL

#endif
