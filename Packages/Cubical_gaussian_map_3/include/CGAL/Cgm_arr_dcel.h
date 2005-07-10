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

#ifndef CGAL_CGM_ARR_DCEL_H
#define CGAL_CGM_ARR_DCEL_H

#include <CGAL/basic.h>
// #include <CGAL/Arr_default_dcel.h>

CGAL_BEGIN_NAMESPACE

/*! Extend the arrangement vertex */
template <class Point_2>
class Cgm_arr_vertex : public CGAL::Arr_vertex_base<Point_2> {
public:
  /*! Vertex location: */
  typedef enum { None, Interior, Edge, Corner } Vertex_location;

private:
  /*! Indicates whether the vertex is interior or on the boundary */
  Vertex_location m_location;

  /*! Indicates whether the vertex is real (as opposed to synthetic) */
  bool m_is_real;

  /*! The id of the containg unit-cube face */
  unsigned int m_face_id;
    
  /*! A pointer to a vertex handle in an adjacent arrangement */
  void * m_adjacent_vertex;

public:
  /*! Constructor */
  Cgm_arr_vertex() :
    m_location(None),
    m_is_real(false),
    m_face_id((unsigned int) -1),
    m_adjacent_vertex(0)
  {}

  /*! Obtain the vertex location */
  Vertex_location get_location() const { return m_location; }

  /*! Set the vertex location */
  void set_location(Vertex_location loc) { m_location = loc; }

  /*! Obtain the flag that indicates whether the vertex is real */
  bool get_is_real() const { return m_is_real; }
    
  /*! Set the flag that indicates whether the vertex is real */
  void set_is_real(bool value) { m_is_real = value; }

  /*! Obtain the index of the containing unit-cube face */
  unsigned int get_face_id() const { return m_face_id; }

  /*! Set the index of the containing unit-cube face */
  void set_face_id(unsigned int face_id) { m_face_id = face_id; }
                   
  /*! Obtain the incident vertex in an adjacent arrangement */    
  void * get_adjacent_vertex() const { return m_adjacent_vertex; }

  /*! Set the incident vertex in an adjacent arrangement */
  void set_adjacent_vertex(void * vertex) { m_adjacent_vertex = vertex; }
};

/*! Extend the arrangement halfedge */
template <class X_monotone_curve_2>
class Cgm_arr_halfedge : public CGAL::Arr_halfedge_base<X_monotone_curve_2> {
private:
  /*! Indicates whether the halfedge is real (as opposed to synthetic) */
  bool m_is_real;

  /*! A mask of the ids of the original arrangements that contributed the
   * halfedge while performing the minkowski sum
   */
  unsigned int m_org_arr_mask;

public:
  Cgm_arr_halfedge() : m_is_real(false), m_org_arr_mask(0x0) {}

  /*! Return true if the halfedge is real */
  bool get_is_real() const { return m_is_real; }

  /*! Set the flag that indicates whether the halfedge is real */
  void set_is_real(bool flag) { m_is_real = flag; }
  
  /*! Add a arrangement to the mask of the original arrangements in the
   * minkowski sum.
   * \param arr_id the id of the added arrangement
   */
  void add_org_arr(unsigned int arr_id) { m_org_arr_mask |= 0x1 << arr_id; }

  /*! Return true iff a given arrangement contributed this halfedge
   * while performing the minkowski sum
   */
  bool is_org_arr(unsigned int arr_id) const
  { return m_org_arr_mask & (0x1 << arr_id); }

  /*! Obtain the mask of the ids of the original arrangements that contributed
   * the halfedge while performing the minkowski sum
   */
  unsigned int get_org_arr_mask() const { return m_org_arr_mask; }
};
  
/*! Extend the arrangement face */

class Cgm_arr_face : public CGAL::Arr_face_base {
public:
  /*! Constructor */
  Cgm_arr_face() {}
};

/*! A new dcel builder with CGM features */
template <class Traits>
class Cgm_arr_dcel :
  public CGAL::Arr_dcel<Cgm_arr_vertex<typename Traits::Point_2>,
                        Cgm_arr_halfedge<typename Traits::X_monotone_curve_2>,
                        Cgm_arr_face>
{
public:
  /*! Constructor */
  Cgm_arr_dcel() {}
};

CGAL_END_NAMESPACE

#endif
