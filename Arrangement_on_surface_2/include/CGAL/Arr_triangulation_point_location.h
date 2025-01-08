// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>

#ifndef CGAL_ARR_TRIANGULATION_POINT_LOCATION_H
#define CGAL_ARR_TRIANGULATION_POINT_LOCATION_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>


/*! \file
 * Definition of the Arr_triangulation_point_location<Arrangement> template.
 */

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
class Arr_triangulation_point_location : public Arrangement_::Observer {
public:
  using Arrangement_2 = Arrangement_;
  using Base_aos = typename Arrangement_2::Base_aos;

  using Geometry_traits_2 = typename Base_aos::Geometry_traits_2;
  using Kernel = typename Geometry_traits_2::Kernel;

  using Vertex_const_handle = typename Base_aos::Vertex_const_handle;
  using Halfedge_const_handle = typename Base_aos::Halfedge_const_handle;
  using Face_const_handle = typename Base_aos::Face_const_handle;
  using Vertex_handle = typename Base_aos::Vertex_handle;
  using Halfedge_handle = typename Base_aos::Halfedge_handle;
  using Face_handle = typename Base_aos::Face_handle;

  using Vertex_const_iterator = typename Base_aos::Vertex_const_iterator;
  using Edge_const_iterator = typename Base_aos::Edge_const_iterator;
  using Face_const_iterator = typename Base_aos::Face_const_iterator;
  using Halfedge_const_iterator = typename Base_aos::Halfedge_const_iterator;
  using Halfedge_around_vertex_const_circulator =
    typename Base_aos::Halfedge_around_vertex_const_circulator;
  using Ccb_halfedge_const_circulator =
    typename Base_aos::Ccb_halfedge_const_circulator;
  using Ccb_halfedge_circulator = typename Base_aos::Ccb_halfedge_circulator;
  using Isolated_vertex_const_iterator =
    typename Base_aos::Isolated_vertex_const_iterator;

  using Point_2 = typename Geometry_traits_2::Point_2;
  using X_monotone_curve_2 = typename Geometry_traits_2::X_monotone_curve_2;

  using Edge_list = std::list<Halfedge_const_handle>;
  using Std_edge_iterator = typename Edge_list::iterator;

  //----------------------------------------------------------
  // Triangulation Types
  //----------------------------------------------------------
  using Vbb = Triangulation_vertex_base_with_info_2<Vertex_const_handle, Kernel>;
  using Vb = Triangulation_hierarchy_vertex_base_2<Vbb>;
  //typedef Triangulation_face_base_with_info_2<CGAL::IO::Color,Kernel>    Fbt;
  using Fb = Constrained_triangulation_face_base_2<Kernel>;
  using TDS = Triangulation_data_structure_2<Vb,Fb>;
  using Itag = Exact_predicates_tag;
  //typedef Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>    CDT;
  using CDT_t = Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>;
  using CDTH = Triangulation_hierarchy_2<CDT_t>;
  using CDT = Constrained_triangulation_plus_2<CDTH>;

  using CDT_Point = typename CDT::Point;
  using CDT_Edge = typename CDT::Edge;
  using CDT_Face_handle = typename CDT::Face_handle;
  using CDT_Vertex_handle = typename CDT::Vertex_handle;
  using CDT_Finite_faces_iterator = typename CDT::Finite_faces_iterator;
  using CDT_Finite_vertices_iterator = typename CDT::Finite_vertices_iterator;
  using CDT_Finite_edges_iterator = typename CDT::Finite_edges_iterator;
  using CDT_Locate_type = typename CDT::Locate_type;

  using Result = Arr_point_location_result<Base_aos>;
  using Result_type = typename Result::Type;

  // Support cpp11::result_of
  using result_type = Result_type;

protected:
  using Traits_adaptor_2 = Arr_traits_basic_adaptor_2<Geometry_traits_2>;

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
    m_traits(nullptr),
    m_ignore_notifications(false),
    m_ignore_remove_edge(false)
  {}

  /*! Constructor from an arrangement.
   * \param arr (in) The arrangement.
   */
  Arr_triangulation_point_location(const Base_aos& arr) :
    Base_aos::Observer(const_cast<Base_aos&>(arr)),
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
  virtual void before_attach(const Base_aos& arr) override
  { m_traits = static_cast<const Traits_adaptor_2*>(arr.geometry_traits()); }

  virtual void after_attach() override { build_triangulation(); }

  virtual void before_detach() override { clear_triangulation(); }

  /// \name Overloaded observer functions on global changes.
  //@{

  /*! Notification after the arrangement has been assigned with another
   * arrangement.
   */
  virtual void after_assign() override {
    clear_triangulation();
    build_triangulation();
  }

  /*! Notification after the arrangement is cleared.
   */
  virtual void after_clear() override {
    clear_triangulation();
    build_triangulation();
  }

  /*! Notification before a global operation modifies the arrangement.
   */
  virtual void before_global_change() override {
    clear_triangulation();
    m_ignore_notifications = true;
  }

  /*! Notification after a global operation is completed.
   */
  virtual void after_global_change() override {
    build_triangulation();
    m_ignore_notifications = false;
  }
  //@}

  /// \name Overloaded observer functions on local changes.
  //@{

  /*! Notification before the removal of an edge.
   * \param e (in) A handle to one of the twin halfedges to be removed.
   */
  virtual void before_remove_edge(Halfedge_handle /* e */) override
  { m_ignore_remove_edge = true; }

  /*! Notification after the creation of a new vertex.
   * \param v (in) A handle to the created vertex.
   */
  virtual void after_create_vertex(Vertex_handle /* v */) override {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after the creation of a new edge.
   * \param e (in) A handle to one of the twin halfedges that were created.
   */
  virtual void after_create_edge(Halfedge_handle /* e */) override {
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
                                Halfedge_handle /* e2 */) override {
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
                                bool /* is_hole */) override {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an outer CCB was created inside a face.
   * \param h (in) A circulator representing the boundary of the new outer CCB.
   */
  virtual void after_add_outer_ccb(Ccb_halfedge_circulator /* h */) override {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an edge was merged.
   * \param e (in) A handle to one of the twin halfedges forming the merged edge.
   */
  virtual void after_merge_edge(Halfedge_handle /* e */) override {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after a face was merged.
   * \param f (in) A handle to the merged face.
   */
  virtual void after_merge_face(Face_handle /* f */) override {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an outer CCB  is moved from one face to another.
   * \param h (in) A circulator representing the boundary of the component.
   */
  virtual void after_move_outer_ccb(Ccb_halfedge_circulator /* h */) override {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification before the removal of a vertex.
   * \param v (in) A handle to the vertex to be deleted.
   */
  virtual void after_remove_vertex() override {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification before the removal of an edge.
   * \param e (in) A handle to one of the twin halfedges to be deleted.
   */
  virtual void after_remove_edge() override {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
    m_ignore_remove_edge = false;
  }

  /*! Notification before the removal of an outer CCB.
   * \param f (in) The face that used to own the outer CCB.
   */
  virtual void after_remove_outer_ccb(Face_handle /* f */) override {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an inner CCB was created inside a face.
   * \param h (in) A circulator representing the boundary of the new inner CCB.
   */
  virtual void after_add_inner_ccb(Ccb_halfedge_circulator /* h */) override {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after an inner CCB is moved from one face to another.
   * \param h (in) A circulator representing the boundary of the component.
   */
  virtual void after_move_inner_ccb(Ccb_halfedge_circulator /* h */) override {
    if (! m_ignore_notifications) {
      clear_triangulation();
      build_triangulation();
    }
  }

  /*! Notification after the removal of an inner CCB.
   * \param f (in) The face that used to contain the inner CCB.
   */
  virtual void after_remove_inner_ccb(Face_handle /* f */) override {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_triangulation();
      build_triangulation();
    }
  }

  // @}

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

#include <CGAL/enable_warnings.h>

#endif
