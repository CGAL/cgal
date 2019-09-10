// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein    <wein@post.tau.ac.il>
//                 Efi Fogel   <efif@post.tau.ac.il>

#ifndef CGAL_ARR_EXTENDED_DCEL_H
#define CGAL_ARR_EXTENDED_DCEL_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The definition of the extended DCEL classes.
 */

#include <CGAL/Arr_dcel_base.h>

namespace CGAL {

/*!
 * \class
 * An extended DCEL vertex with auxiliary data field.
 */
template <class VertexBase_, typename VertexData_>
class Arr_extended_vertex : public VertexBase_
{
  typedef Arr_extended_vertex<VertexBase_, VertexData_>  Self;
  typedef VertexBase_                                    Base;

public:
  typedef VertexData_                                     Data;

private:

  Data       m_data;       // The auxiliary data field.

public:

  /*! Get the auxiliary data (const version). */
  const Data& data () const
  {
    return (m_data);
  }

  /*! Get the auxiliary data (non-const version). */
  Data& data ()
  {
    return (m_data);
  }

  /*! Set the auxiliary data. */
  void set_data (const Data& data)
  {
    m_data = data;
  }

  /*! Assign from another vertex. */
  virtual void assign (const Base& v)
  {
    Base::assign (v);

    const Self&  ex_v = static_cast<const Self&>(v);
    m_data = ex_v.m_data;
  }
};

/*!
 * \class
 * An extended DCEL halfedge with auxiliary data field.
 */
template <class HalfedgeBase_, typename HalfedgeData_>
class Arr_extended_halfedge : public HalfedgeBase_
{
  typedef Arr_extended_halfedge<HalfedgeBase_, HalfedgeData_>  Self;
  typedef HalfedgeBase_                                        Base;

public:
  typedef HalfedgeData_                                        Data;

private:

  Data       m_data;       // The auxiliary data field.

public:

  /*! Get the auxiliary data (const version). */
  const Data& data () const
  {
    return (m_data);
  }

  /*! Get the auxiliary data (non-const version). */
  Data& data ()
  {
    return (m_data);
  }

  /*! Set the auxiliary data. */
  void set_data (const Data& data)
  {
    m_data = data;
  }

  /*! Assign from another vertex. */
  virtual void assign (const Base& he)
  {
    Base::assign (he);

    const Self&  ex_he = static_cast<const Self&>(he);
    m_data = ex_he.m_data;
  }
};

/*!
 * \class
 * An extended DCEL face with auxiliary data field.
 */
template <class FaceBase_, typename FaceData_>
class Arr_extended_face : public FaceBase_
{
  typedef Arr_extended_face<FaceBase_, FaceData_>  Self;
  typedef FaceBase_                                Base;

public:
  typedef FaceData_                                Data;

private:

  Data       m_data;       // The auxiliary data field.

public:

  /*! Get the auxiliary data (const version). */
  const Data& data () const
  {
    return (m_data);
  }

  /*! Get the auxiliary data (non-const version). */
  Data& data ()
  {
    return (m_data);
  }

  /*! Set the auxiliary data. */
  void set_data (const Data& data)
  {
    m_data = data;
  }

  /*! Assign from another vertex. */
  virtual void assign (const Base& f)
  {
    Base::assign (f);

    const Self&  ex_f = static_cast<const Self&>(f);
    m_data = ex_f.m_data;
  }
};

/*! \class
 * A DCEL class whose faces are extended with an auxiliary data field.
 * The Traits parameter corresponds to a geometric traits class, which 
 * defines the Point_2 and X_monotone_curve_2 types.
 * The FaceData parameter specifies the object type stored with each face.
 */
template <class Traits_, typename FaceData_,
	  class VertexBase_ = Arr_vertex_base<typename Traits_::Point_2>,
	  class HalfedgeBase_ =
   	             Arr_halfedge_base<typename Traits_::X_monotone_curve_2>,
	  class FaceBase_ = Arr_face_base>
class Arr_face_extended_dcel :
  public Arr_dcel_base<VertexBase_,
		       HalfedgeBase_,
		       Arr_extended_face<FaceBase_, FaceData_> >
{
public:

  typedef FaceData_                    Face_data;

  /*! \struct
   * An auxiliary structure for rebinding the DCEL with a new traits class.
   */
  template<typename T>
  class rebind
  {
    typedef typename VertexBase_::template rebind
                        <typename T::Point_2>              Rebind_vertex;
    typedef typename Rebind_vertex::other                  Vertex_base;
    typedef typename HalfedgeBase_::template rebind
                        <typename T::X_monotone_curve_2>   Rebind_halfedge;
    typedef typename Rebind_halfedge::other                Halfedge_base;
   
  public:

    typedef Arr_face_extended_dcel<T,
				   Face_data,
				   Vertex_base,
				   Halfedge_base,
				   FaceBase_>              other;
  };

  /*! Default constructor. */
  Arr_face_extended_dcel ()
  {}

  /*! Destructor. */
  virtual ~Arr_face_extended_dcel ()
  {}
};

/*! \class
 * A DCEL class whose features are extended with auxiliary data fields.
 * The Traits parameter corresponds to a geometric traits class, which 
 * defines the Point_2 and X_monotone_curve_2 types.
 * The VertexData, HalfedgeData and FaceData parameter specify the object types
 * stored with each vertex, halfegde and face, respectively.
 */
template <class Traits_,
	  typename VertexData_, typename HalfedgeData_, typename FaceData_,
	  class VertexBase_ = Arr_vertex_base<typename Traits_::Point_2>,
	  class HalfedgeBase_ =
   	             Arr_halfedge_base<typename Traits_::X_monotone_curve_2>,
	  class FaceBase_ = Arr_face_base>
class Arr_extended_dcel :
  public Arr_dcel_base<Arr_extended_vertex<VertexBase_, VertexData_>,
		       Arr_extended_halfedge<HalfedgeBase_, HalfedgeData_>,
		       Arr_extended_face<FaceBase_, FaceData_> >
{
public:

  typedef VertexData_                  Vertex_data;
  typedef HalfedgeData_                Halfedge_data;
  typedef FaceData_                    Face_data;

  /*! \struct
   * An auxiliary structure for rebinding the DCEL with a new traits class.
   */
  template<typename T>
  class rebind
  {
    typedef typename VertexBase_::template rebind
                        <typename T::Point_2>              Rebind_vertex;
    typedef typename Rebind_vertex::other                  Vertex_base;
    typedef typename HalfedgeBase_::template rebind
                        <typename T::X_monotone_curve_2>   Rebind_halfedge;
    typedef typename Rebind_halfedge::other                Halfedge_base;
   
  public:

    typedef Arr_extended_dcel<T,
			      Vertex_data,
			      Halfedge_data,
			      Face_data,
			      Vertex_base,
			      Halfedge_base,
			      FaceBase_>                   other;
  };

  /*! Default constructor. */
  Arr_extended_dcel ()
  {}

  /*! Destructor. */
  virtual ~Arr_extended_dcel ()
  {}
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif 
