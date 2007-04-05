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
// $URL$
// $Id$
// 
//
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_POLYHEDRAL_CGM_POLYHEDRON_3_H
#define CGAL_POLYHEDRAL_CGM_POLYHEDRON_3_H

/*! \file
 * Related definition of a Polyhedron_3 data structure, an instance of which
 * can be used to initialize a Polyhedral_cgm data structure.
 * A user can construct a Polyhedral_cgm data structure in two ways as follows:
 * 1. Providing a sequence of (geometric) points in three dimensions and a
 * sequence of indices into the points that indicate the facets of the
 * polyhedron. An index equal to -1 indicates the end of a facet. Internally,
 * a temporary Polyhedron_3 object is constructed from the input. Then, the
 * Polyhedral_cgm object is constructed from the Polyhedron_3 object. Finally,
 * the temporary Polyhedron_3 object is destructed.
 * 2. Providing a Polyhedron_3 object directly. In this case, the type
 * is constrained. The user is free to define her/his own vertex, halfedge,
 * and face types, but they must derive from the types provided in this file,
 * namely:
 *   Polyhedral_cgm_polyhedron_3_vertex,
 *   Polyhedral_cgm_polyhedron_3_halfedge, and
 *   Polyhedral_cgm_polyhedron_3_face.
 * Notice that the latter is parameterized with a Cgm type, the same type
 * as the type Polyhedral_cgm that, and object of which the user intends to
 * construct.
 */

#include <CGAL/basic.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_face_base.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>

CGAL_BEGIN_NAMESPACE

/*! The extended Polyhedron vertex type */
template <class T_Refs, class T_Point>
class Polyhedral_cgm_polyhedron_3_vertex :
  public HalfedgeDS_vertex_base<T_Refs, Tag_true, T_Point>
{
private:
  typedef HalfedgeDS_vertex_base<T_Refs, Tag_true, T_Point> Base;

  /*! Indicates whether it is a marked vertex */
  bool m_marked;

public:
  typedef typename Base::Point                                  Point;

  /*! Constructor */
  Polyhedral_cgm_polyhedron_3_vertex() : Base(), m_marked(false) {}

  /*! Constructor */
  Polyhedral_cgm_polyhedron_3_vertex(const Point & p) :
    Base(p), m_marked(false) {}

  /*! Obtain the mutable (geometrical) point. Delegate */
  Point & point() { return Base::point(); }

  /*! Obtain the constant (geometrical) point. Delegate */
  const Point & point () const { return Base::point(); }

  /*! Set the "marked" flag */
  void set_marked(bool marked) { m_marked = marked; }

  /*! Obtain the "marked" flag */
  bool get_marked() const { return m_marked; }
};

/*! The extended Polyhedron halfedge type */
template <class T_Refs>
class Polyhedral_cgm_polyhedron_3_halfedge :
  public HalfedgeDS_halfedge_base<T_Refs>
{
private:
  /*! Indicates that the halfedge has been processed already */
  bool m_processed;

  /*! Indicates whether it is a marked vertex */
  bool m_marked;

public:
  /*! Constructor */
  Polyhedral_cgm_polyhedron_3_halfedge() :
    m_processed(false), m_marked(false) {}

  /*! Set the flag */
  void set_processed(bool processed) { m_processed = processed; }

  /*! Obtain the flag */
  bool processed() const { return m_processed; }

  /*! Set the "marked" flag */
  void set_marked(bool marked) { m_marked = marked; }

  /*! Obtain the "marked" flag */
  bool get_marked() const { return m_marked; }
};

/*! The extended Polyhedron face type */
template <class T_Refs, class T_Plane, class Cgm>
class Polyhedral_cgm_polyhedron_3_face :
  public HalfedgeDS_face_base<T_Refs, Tag_true, T_Plane>
{
private:
  typedef HalfedgeDS_face_base<T_Refs, Tag_true, T_Plane> Base;
  typedef typename Cgm::Projected_normal                Projected_normal;
  
  /*! The normal projected onto the unit cube */
  Projected_normal m_dual;

  /*! Indicates whether it is a marked face */
  bool m_marked;

public:
  typedef typename Base::Plane                          Plane;

  /*! Constructor */
  Polyhedral_cgm_polyhedron_3_face() : m_marked(false) {}

  /*! Obtain the mutable plane. Delegate */
  Plane & plane() { return Base::plane(); }

  /*! Obtain the constant plane. Delegate */
  const Plane & plane() const { return Base::plane(); }

  /*! Obtain the dual structure */
  Projected_normal & get_dual() { return m_dual; }

  /*! Set the "marked" flag */
  void set_marked(bool marked) { m_marked = marked; }

  /*! Obtain the "marked" flag */
  bool get_marked() const { return m_marked; }

  /*! Compute the central projection */
  void compute_projection(void)
  {
    m_dual.compute_projection(plane());
    m_dual.set_is_real(true);
  }
};

/*! The "items" type. A model of the PolyhedralCgmPolyhedronItems_3 concept,
 * which is a refinment of the PolyhedronItems_3 concept. Its base class
 * Polyhedron_items_3, a model of the latter concept, provides definitions of
 * vertices with points, halfedges, and faces with normal equations. We extend
 * the definition of each one of the three items with the necessary data
 * required to construct a Polyhedral_cgm object from a Polyhedron_3 object,
 * the type of which, namely Polyhedron_3, is instantiated with this extended
 * items type.
 */
template <class Cgm>
struct Polyhedral_cgm_polyhedron_items : public Polyhedron_items_3 {
  template <class T_Refs, class T_Traits>
  struct Vertex_wrapper {
    typedef typename T_Traits::Point_3                              Point_3;
    typedef Polyhedral_cgm_polyhedron_3_vertex<T_Refs, Point_3>     Vertex;
  };
  template <class T_Refs, class T_Traits>
  struct Halfedge_wrapper {
    typedef Polyhedral_cgm_polyhedron_3_halfedge<T_Refs>            Halfedge;
  };
  template <class T_Refs, class T_Traits>
  struct Face_wrapper {
    typedef typename T_Traits::Plane_3                              Plane_3;
    typedef Polyhedral_cgm_polyhedron_3_face<T_Refs, Plane_3, Cgm>  Face;
  };
};

/*! The default polyhedron type. If the Polyhedral_cgm object is indirectly
 * constructed from the points and the facets provided as indices, then a
 * temporary object of type Polyhedral_cgm_default_polyhedron_3 is constructed
 * internally, and used to represent the polyhedron. Similarly, if the user
 * provides a reference to a polyhedron object as input for the construction
 * of the Polyhedral_cgm object, and she/he has no need to extend the
 * polyhedron features, this type should be used to represent the polyhedron.
 * However, if the user need to extend the vertex, halfedge, or face of the
 * polyhedron, she/he must extend the appropriate type(s), define a new items
 * type that is based on the extended types, and define a new polyhedron type
 * based on the new items type.
 */
template <class Cgm>
struct Polyhedral_cgm_polyhedron_3 :
  public Polyhedron_3<Polyhedron_traits_with_normals_3<typename Cgm::Kernel>,
                      Polyhedral_cgm_polyhedron_items<Cgm> >
{
  /*! Constructor */
  Polyhedral_cgm_polyhedron_3() {}
};

CGAL_END_NAMESPACE

#endif
