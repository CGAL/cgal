// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein    <wein@post.tau.ac.il>
//            Efi Fogel   <efif@post.tau.ac.il>

#ifndef CGAL_ARR_EXTENDED_DCEL_H
#define CGAL_ARR_EXTENDED_DCEL_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The definition of the extended DCEL classes.
 */

#include <CGAL/Arr_dcel_base.h>

namespace CGAL {

/*! \class
 * An extended DCEL vertex with auxiliary data field.
 */
template <typename VertexBase, typename VertexData>
class Arr_extended_vertex : public VertexBase {
  using Vertex_base = VertexBase;
  using Vertex_data = VertexData;

  using Self = Arr_extended_vertex<Vertex_base, Vertex_data>;

public:
  typedef Vertex_data                                    Data;

private:
  Data m_data;       // The auxiliary data field.

public:
  /*! Get the auxiliary data (const version). */
  const Data& data() const { return m_data; }

  /*! Get the auxiliary data (non-const version). */
  Data& data() { return m_data; }

  /*! Set the auxiliary data. */
  void set_data(const Data& data) { m_data = data; }

  /*! Assign from another vertex. */
  virtual void assign(const Vertex_base& v) {
    Vertex_base::assign(v);
    const Self& ex_v = static_cast<const Self&>(v);
    m_data = ex_v.m_data;
  }

  template <typename Point_>
  struct rebind {
    using Point_2 = Point_;
    using other = typename Vertex_base::template rebind<Point_2>;
  };
};

/*! \class
 * An extended DCEL halfedge with auxiliary data field.
 */
template <typename HalfedgeBase, typename HalfedgeData>
class Arr_extended_halfedge : public HalfedgeBase {
  using Halfedge_base = HalfedgeBase;
  using Halfedge_data = HalfedgeData;

  using Self = Arr_extended_halfedge<Halfedge_base, Halfedge_data>;

public:
  typedef Halfedge_data                                       Data;

private:
  Data m_data;       // The auxiliary data field.

public:
  /*! Get the auxiliary data (const version). */
  const Data& data() const { return m_data; }

  /*! Get the auxiliary data (non-const version). */
  Data& data() { return m_data; }

  /*! Set the auxiliary data. */
  void set_data(const Data& data) { m_data = data; }

  /*! Assign from another halfedge. */
  virtual void assign(const Halfedge_base& he) {
    Halfedge_base::assign(he);
    const Self& ex_he = static_cast<const Self&>(he);
    m_data = ex_he.m_data;
  }

  template <typename XMonotoneCurve>
  struct rebind {
    using X_monotonote_curve_2 = XMonotoneCurve;
    using other = typename Halfedge_base::template rebind<X_monotonote_curve_2>;
  };
};

/*! \class
 * An extended DCEL face with auxiliary data field.
 */
template <typename FaceBase, typename FaceData>
class Arr_extended_face : public FaceBase {
  using Face_base = FaceBase;
  using Face_data = FaceData;

  using Self = Arr_extended_face<Face_base, Face_data>;

public:
  typedef Face_data                               Data;

private:
  Data m_data;       // The auxiliary data field.

public:
  /*! Get the auxiliary data (const version). */
  const Data& data() const { return m_data; }

  /*! Get the auxiliary data (non-const version). */
  Data& data() { return m_data; }

  /*! Set the auxiliary data. */
  void set_data(const Data& data) { m_data = data; }

  /*! Assign from another face. */
  virtual void assign(const Face_base& f) {
    Face_base::assign(f);
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
template <typename Traits_, typename FaceData,
          typename VertexBase = Arr_vertex_base<typename Traits_::Point_2>,
          typename HalfedgeBase =
            Arr_halfedge_base<typename Traits_::X_monotone_curve_2>,
          typename FaceBase = Arr_face_base>
class Arr_face_extended_dcel :
  public Arr_dcel_base<VertexBase, HalfedgeBase,
                       Arr_extended_face<FaceBase, FaceData>> {
public:
  using Face_base = FaceBase;
  using Face_data = FaceData;

  /*! \struct
   * An auxiliary structure for rebinding the DCEL with a new traits class.
   */
  template <typename T>
  class rebind {
  private:
    using Pnt = typename T::Point_2;
    using Xcv = typename T::X_monotone_curve_2;
    using Rebind_vertex = typename VertexBase::template rebind<Pnt>;
    using Vertex_other = typename Rebind_vertex::other;
    using Rebind_halfedge = typename HalfedgeBase::template rebind<Xcv>;
    using Halfedge_other = typename Rebind_halfedge::other;

  public:
    using other = Arr_face_extended_dcel<T, Face_data, Vertex_other,
                                         Halfedge_other, Face_base>;
  };

  /*! Default constructor. */
  Arr_face_extended_dcel() {}

  /*! Destructor. */
  virtual ~Arr_face_extended_dcel() {}
};

/*! \class
 * A DCEL class whose features are extended with auxiliary data fields.
 * The Traits parameter corresponds to a geometric traits class, which
 * defines the Point_2 and X_monotone_curve_2 types.
 * The VertexData, HalfedgeData and FaceData parameter specify the object types
 * stored with each vertex, halfedge and face, respectively.
 */
template <typename Traits_,
          typename VertexData, typename HalfedgeData, typename FaceData,
          typename VertexBase = Arr_vertex_base<typename Traits_::Point_2>,
          typename HalfedgeBase =
            Arr_halfedge_base<typename Traits_::X_monotone_curve_2>,
          typename FaceBase = Arr_face_base>
class Arr_extended_dcel :
  public Arr_dcel_base<Arr_extended_vertex<VertexBase, VertexData>,
                       Arr_extended_halfedge<HalfedgeBase, HalfedgeData>,
                       Arr_extended_face<FaceBase, FaceData>> {
public:
  using Vertex_data = VertexData;
  using Halfedge_data = HalfedgeData;
  using Face_data = FaceData;
  using Vertex_base = VertexBase;
  using Halfedge_base = HalfedgeBase;
  using Face_base = FaceBase;

  /*! \struct
   * An auxiliary structure for rebinding the DCEL with a new traits class.
   */
  template <typename T>
  struct rebind {
  private:
    using Pnt = typename T::Point_2;
    using Xcv = typename T::X_monotone_curve_2;
    using Rebind_vertex = typename VertexBase::template rebind<Pnt>;
    using Vertex_other = typename Rebind_vertex::other;
    using Rebind_halfedge = typename HalfedgeBase::template rebind<Xcv>;
    using Halfedge_other = typename Rebind_halfedge::other;

  public:
    using other = Arr_extended_dcel<T,
                                    Vertex_data, Halfedge_data, Face_data,
                                    Vertex_other, Halfedge_other, Face_base>;
  };

  /*! Default constructor. */
  Arr_extended_dcel() {}

  /*! Destructor. */
  virtual ~Arr_extended_dcel() {}
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
