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
//
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>

#ifndef CGAL_ARR_TRIANGULATION_POINT_LOCATION_H
#define CGAL_ARR_TRIANGULATION_POINT_LOCATION_H

/*! \file
 * Definition of the Arr_triangulation_point_location<Arrangement> template.
 */

#include <CGAL/Arr_observer.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_point_location_result.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

namespace CGAL {

/*! \class
 * A class that answers point-location queries on an arrangement using the
 * triangulation algorithm.
 */
template <typename Arrangement_>
class Arr_triangulation_point_location : public Arr_observer<Arrangement_>
{
public:
  typedef Arrangement_                                  Arrangement_2;

  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Geometry_traits_2::Kernel            Kernel;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Vertex_handle		Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle	Halfedge_handle;
  typedef typename Arrangement_2::Face_handle		Face_handle;

  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Arrangement_2::Edge_const_iterator   Edge_const_iterator;
  typedef typename Arrangement_2::Face_const_iterator   Face_const_iterator;
  typedef typename Arrangement_2::Halfedge_const_iterator
    Halfedge_const_iterator;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_circulator
    Ccb_halfedge_circulator;
  typedef typename Arrangement_2::Isolated_vertex_const_iterator
    Isolated_vertex_const_iterator;

  typedef typename Geometry_traits_2::Point_2            Point_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2 X_monotone_curve_2;

  typedef std::list<Halfedge_const_handle>               Edge_list;
  typedef typename Edge_list::iterator                   Std_edge_iterator;

  //----------------------------------------------------------
  // Triangulation Types
  //----------------------------------------------------------
  typedef Triangulation_vertex_base_with_info_2<Vertex_const_handle, Kernel>
    Vbb;
  typedef Triangulation_hierarchy_vertex_base_2<Vbb>                  Vb;
  //typedef Triangulation_face_base_with_info_2<CGAL::Color,Kernel>    Fbt;
  typedef Constrained_triangulation_face_base_2<Kernel>               Fb;
  typedef Triangulation_data_structure_2<Vb,Fb>                       TDS;
  typedef Exact_predicates_tag                                        Itag;
  //typedef Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>    CDT;
  typedef Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>     CDT_t;
  typedef Triangulation_hierarchy_2<CDT_t>                            CDTH;
  typedef Constrained_triangulation_plus_2<CDTH>                      CDT;

  typedef typename CDT::Point                    CDT_Point;
  typedef typename CDT::Edge                     CDT_Edge;
  typedef typename CDT::Face_handle              CDT_Face_handle;
  typedef typename CDT::Vertex_handle            CDT_Vertex_handle;
  typedef typename CDT::Finite_faces_iterator    CDT_Finite_faces_iterator;
  typedef typename CDT::Finite_vertices_iterator CDT_Finite_vertices_iterator;
  typedef typename CDT::Finite_edges_iterator    CDT_Finite_edges_iterator;
  typedef typename CDT::Locate_type              CDT_Locate_type;

  typedef Arr_point_location_result<Arrangement_2>       Result;
  typedef typename Result::Type                          Result_type;

  // Support cpp11::result_of
  typedef Result_type                                    result_type;

protected:
  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2>  Traits_adaptor_2;

  // Data members:
  const Traits_adaptor_2* m_traits;     // Its associated traits object.
  bool m_ignore_notifications;
  bool m_ignore_remove_edge;
  CDT m_cdt;

  template<typename T>
  Result_type make_result(T t) const { return Result::make_result(t); }
  inline Result_type default_result() const { return Result::default_result(); }

public:
  /*! Default constructor. */
  Arr_triangulation_point_location() :
    m_traits(NULL),
    m_ignore_notifications(false),
    m_ignore_remove_edge(false)
  {}

  /*! Constructor from an arrangement.
   * \param arr (in) The arrangement.
   */
  Arr_triangulation_point_location(const Arrangement_2& arr) :
    Arr_observer<Arrangement_2>(const_cast<Arrangement_2&>(arr)),
    m_traits(static_cast<const Traits_adaptor_2*>(arr.geometry_traits())),
    m_ignore_notifications(false),
    m_ignore_remove_edge(false)
  { build_triangulation(); }

  /*! Locate the arrangement feature containing the given point.
   * \param p (in) The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type locate(const Point_2& p) const;

  //Observer functions that are relevant to overload
  //-------------------------------------------------

  /*! Attach an arrangement.
   * \param arr (in) The arrangement.
   */
  virtual void before_attach(const Arrangement_2& arr)
  { m_traits = static_cast<const Traits_adaptor_2*>(arr.geometry_traits()); }

  virtual void after_attach() { build_triangulation(); }

  virtual void before_detach() { clear_triangulation(); }

  /// \name Overloaded observer functions on global changes.
  //@{

  /*! Notification after the arrangement has been assigned with another
   * arrangement.
   */
  virtual void after_assign()
  {
    clear_triangulation();
    build_triangulation();
  }

  /*! Notification after the arrangement is cleared.
   */
  virtual void after_clear()
  {
    clear_triangulation();
    build_triangulation();
  }

  /*! Notification before a global operation modifies the arrangement.
   */
  virtual void before_global_change()
  {
    clear_triangulation();
    m_ignore_notifications = true;
  }

  /*! Notification after a global operation is completed.
   */
  virtual void after_global_change()
  {
    build_triangulation();
    m_ignore_notifications = false;
  }
  //@}

  /// \name Overloaded observer functions on local changes.
  //@{

  /*! Notification before the removal of an edge.
   * \param e (in) A handle to one of the twin halfedges to be removed.
   */
  virtual void before_remove_edge(Halfedge_handle /* e */)
  { m_ignore_remove_edge = true; }

  /*! Notification after the creation of a new vertex.
   * \param v (in) A handle to the created vertex.
   */
  virtual void after_create_vertex(Vertex_handle /* v */)
  {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after the creation of a new edge.
   * \param e (in) A handle to one of the twin halfedges that were created.
   */
  virtual void after_create_edge(Halfedge_handle /* e */)
  {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an edge was split.
   * \param e1 (in) A handle to one of the twin halfedges forming the first edge.
   * \param e2 (in) A handle to one of the twin halfedges forming the second edge.
   */
  virtual void after_split_edge(Halfedge_handle /* e1 */,
                                Halfedge_handle /* e2 */)
  {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after a face was split.
   * \param f (in) A handle to the face we have just split.
   * \param new_f (in) A handle to the new face that has been created.
   * \param is_hole (in) Whether the new face forms a hole inside f.
   */
  virtual void after_split_face(Face_handle /* f */,
                                Face_handle /* new_f */,
                                bool /* is_hole */)
  {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an outer CCB was created inside a face.
   * \param h (in) A circulator representing the boundary of the new outer CCB.
   */
  virtual void after_add_outer_ccb(Ccb_halfedge_circulator /* h */)
  {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an edge was merged.
   * \param e (in) A handle to one of the twin halfedges forming the merged edge.
   */
  virtual void after_merge_edge(Halfedge_handle /* e */)
  {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after a face was merged.
   * \param f (in) A handle to the merged face.
   */
  virtual void after_merge_face(Face_handle /* f */)
  {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an outer CCB  is moved from one face to another.
   * \param h (in) A circulator representing the boundary of the component.
   */
  virtual void after_move_outer_ccb(Ccb_halfedge_circulator /* h */)
  {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notificaion before the removal of a vertex.
   * \param v (in) A handle to the vertex to be deleted.
   */
  virtual void after_remove_vertex()
  {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification before the removal of an edge.
   * \param e (in) A handle to one of the twin halfedges to be deleted.
   */
  virtual void after_remove_edge()
  {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
    m_ignore_remove_edge = false;
  }

  /*! Notification before the removal of an outer CCB.
   * \param f (in) The face that used to own the outer CCB.
   */
  virtual void after_remove_outer_ccb(Face_handle /* f */)
  {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an inner CCB was created inside a face.
   * \param h (in) A circulator representing the boundary of the new inner CCB.
   */
  virtual void after_add_inner_ccb(Ccb_halfedge_circulator /* h */)
  {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an inner CCB is moved from one face to another.
   * \param h (in) A circulator representing the boundary of the component.
   */
  virtual void after_move_inner_ccb(Ccb_halfedge_circulator /* h */)
  {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notificaion after the removal of an inner CCB.
   * \param f (in) The face that used to contain the inner CCB.
   */
  virtual void after_remove_inner_ccb(Face_handle /* f */)
  {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_triangulation();
      build_triangulation();
    }
  }

protected:
  /*! Locate the arrangement feature containing the given point in the
   * unbounded face(s).
   * \param p (in) The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle
   *         representing an unbounded face or a Vertex_const_handle
   *         representing an isolated vertex.
   */
  result_type locate_in_unbounded(const Point_2& p) const;

  void clear_triangulation();
  void build_triangulation();
};

} //namespace CGAL

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Arr_triangulation_pl_functions.h>

#endif
