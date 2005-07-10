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
// $Revision$
// $Name$
//
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_POLYHEDRAL_CGM_ARR_DCEL_H
#define CGAL_POLYHEDRAL_CGM_ARR_DCEL_H

#include <CGAL/basic.h>

#include "CGAL/Cgm_arr_dcel.h"

CGAL_BEGIN_NAMESPACE

/*! Extend the arrangement vertex */
template <class Point_2>
class Polyhedral_cgm_arr_vertex : public Cgm_arr_vertex<Point_2> {
public:
  /*! Constructor */
  Polyhedral_cgm_arr_vertex() {}
};

/*! Extend the arrangement halfedge */
template <class X_monotone_curve_2>
class Polyhedral_cgm_arr_halfedge : public Cgm_arr_halfedge<X_monotone_curve_2> {
public:
  /*! Constructor */
  Polyhedral_cgm_arr_halfedge() {}
};

/*! Extend the arrangement face */
template <class Point_3>
class Polyhedral_cgm_arr_face : public Cgm_arr_face {
private:
  /*! The original point of the polyhedron */
  Point_3 m_point;

  /*! Indicates that the point has been set already */
  bool m_is_set;

public:
  /*! Constructor */
  Polyhedral_cgm_arr_face() : m_is_set(false) { }
    
  /*! Set the 3D point of the original polyhedron */
  void set_point(const Point_3 & point)
  {
    m_point = point;
    m_is_set = true;
  }

  /*! Obtain the 3D point of the original polyhedron */
  const Point_3 & get_point() const { return m_point; }

  /*! \brief returns true iff the point has been set already */
  bool get_is_set() const { return m_is_set; }

  /*! \brief resets the flag  */
  void set_is_set(bool flag) { m_is_set = flag; }
};

/*! A new dcel builder with CGM features */
template <class Traits>
class Polyhedral_cgm_arr_dcel :
  public CGAL::Arr_dcel<Polyhedral_cgm_arr_vertex<typename Traits::Point_2>,
                        Polyhedral_cgm_arr_halfedge<typename Traits::X_monotone_curve_2>,
                        Polyhedral_cgm_arr_face<typename Traits::Point_3> >
{
public:
  /*! Constructor */
  Polyhedral_cgm_arr_dcel() {}
};

CGAL_END_NAMESPACE

#endif
