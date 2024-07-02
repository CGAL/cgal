// Copyright (c) 2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
//                 Ron Wein     <wein@post.tau.ac.il>

#ifndef CGAL_ARR_LANDMARKS_GENERATOR_BASE_H
#define CGAL_ARR_LANDMARKS_GENERATOR_BASE_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
* Definition of the Arr_landmarks_generator_base<Arrangement> template.
*/
#include <CGAL/Arr_point_location_result.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_point_location/Arr_lm_nearest_neighbor.h>
#include <CGAL/Arr_batched_point_location.h>
#include <CGAL/algorithm.h>

#include <list>
#include <algorithm>
#include <vector>

namespace CGAL {

/*! \class Arr_landmarks_generator_base
* A pure virtual class that handles the changes in a the arrangement in a generic
* manner(i.e., it rebuilds the search structure for whenever a small change takes
* place in the arrangement). It serves as a base class for some of the generator
* classes.
* The classes that inherite from it must implement at least one virtual
* function called "void _create_point_list(Point_list &)" which
* actually creates the list of landmarks.
*/
template <typename Arrangement_,
          typename Nearest_neighbor_ =
            Arr_landmarks_nearest_neighbor<Arrangement_> >
class Arr_landmarks_generator_base : public Arrangement_::Observer {
public:
  using Arrangement_2 = Arrangement_;
  using Base_aos = typename Arrangement_2::Base_aos;

  using Nearest_neighbor = Nearest_neighbor_;

  using Geometry_traits_2 = typename Base_aos::Geometry_traits_2;
  using Vertex_const_handle = typename Base_aos::Vertex_const_handle;
  using Halfedge_const_handle = typename Base_aos::Halfedge_const_handle;
  using Face_const_handle = typename Base_aos::Face_const_handle;
  using Vertex_handle = typename Base_aos::Vertex_handle;
  using Halfedge_handle = typename Base_aos::Halfedge_handle;
  using Face_handle = typename Base_aos::Face_handle;
  using Vertex_const_iterator = typename Base_aos::Vertex_const_iterator;
  using Ccb_halfedge_circulator = typename Base_aos::Ccb_halfedge_circulator;

  using Point_2 = typename Base_aos::Point_2;
  using X_monotone_curve_2 = typename Base_aos::X_monotone_curve_2;

  using NN_Point_2 = typename Nearest_neighbor::NN_Point_2;
  using NN_Points_set = std::list<NN_Point_2>;

  using Points_set = std::vector<Point_2>;

  using PL_result = Arr_point_location_result<Base_aos>;
  using PL_result_type = typename PL_result::Type;

  using PL_pair = std::pair<Point_2, PL_result_type>;
  using Pairs_set = std::vector<PL_pair>;
  using Pairs_iterator = typename std::vector<PL_pair>::iterator;

private:
  using Self = Arr_landmarks_generator_base<Arrangement_2, Nearest_neighbor>;

protected:
  using Traits_adaptor_2 = Arr_traits_basic_adaptor_2<Geometry_traits_2>;

  // Data members:
  const Traits_adaptor_2* m_traits;  // The associated traits object.
  Nearest_neighbor nn;      // The associated nearest neighbor object.
  bool m_ignore_notifications;
  bool m_ignore_remove_edge;
  bool updated;
  int num_small_not_updated_changes;

  template<typename T>
  PL_result_type pl_make_result(T t) { return PL_result::make_result(t); }
  inline PL_result_type pl_default_result()
  { return PL_result::default_result(); }

public:
  bool is_empty() const { return nn.is_empty(); }

private:
  /*! Copy constructor - not supported. */
  Arr_landmarks_generator_base(const Self&);

  /*! Assignment operator - not supported. */
  Self& operator=(const Self& );

public:
  /*! Constructor from an arrangement.
   * \param arr (in) The arrangement.
   */
  Arr_landmarks_generator_base(const Base_aos& arr) :
    Base_aos::Observer(const_cast<Base_aos&>(arr)),
    m_traits(static_cast<const Traits_adaptor_2*>(arr.geometry_traits())),
    m_ignore_notifications(false),
    m_ignore_remove_edge(false),
    updated(false),
    num_small_not_updated_changes(0)
  {
    // One needs to call build_landmark_set() in the constructor of any
    // inherited class.
  }

  /*! Create the landmarks set (choosing the landmarks) ,
   * and saving them in the nearest-neighbor search structure.
   */
  virtual void build_landmark_set() {
    // Create the landmark points.
    NN_Points_set nn_points;
    _create_nn_points_set(nn_points);

    // Update the search structure.
    nn.clear();
    nn.init(nn_points.begin(), nn_points.end());

    num_small_not_updated_changes = 0;
    updated = true;
  }

  /*! clear the set of landmarks.
   */
  virtual void clear_landmark_set() {
    nn.clear();
    num_small_not_updated_changes = 0;
    updated = false;
  }

  /*! Obtain the nearest neighbor (landmark) to the given point.
   * \param p The query point.
   * \param obj Output: The location of the nearest landmark point in the
   *                    arrangement (a vertex, halfedge, or face handle).
   * \return The nearest landmark point.
   */
  virtual Point_2 closest_landmark(const Point_2& p, PL_result_type& obj) {
    CGAL_assertion(updated);
    return (nn.find_nearest_neighbor(p, obj));
  }

  /// \name Overloaded observer functions on global changes.
  //@{

  /*! Notification before the arrangement is assigned with the content of
   * another arrangement.
   * \param arr The other arrangement. Notice that the arrangement type is the
   *            type used to instantiate the observer, which is conveniently
   *            defined as `Arrangement_2::Base_aos`.
   */
  virtual void before_assign(const Base_aos& arr) override {
    this->clear_landmark_set();
    m_traits = static_cast<const Traits_adaptor_2*>(arr.geometry_traits());
    m_ignore_notifications = true;
  }

  /*! Notification after the arrangement has been assigned with another
   * arrangement.
   */
  virtual void after_assign() override {
    this->build_landmark_set();
    m_ignore_notifications = false;
  }

  /*! Notification before the observer is attached to an arrangement.
   * \param arr The arrangement we are about to attach the observer to.
   */
  virtual void before_attach(const Base_aos& arr) override {
    this->clear_landmark_set();
    m_traits = static_cast<const Traits_adaptor_2*>(arr.geometry_traits());
    m_ignore_notifications = true;
  }

  /*! Notification after the observer has been attached to an arrangement.
   */
  virtual void after_attach() override {
    this->build_landmark_set();
    m_ignore_notifications = false;
  }

  /*! Notification before the observer is detached from the arrangement.
   */
  virtual void before_detach() override { this->clear_landmark_set(); }

  /*! Notification after the arrangement is cleared.
   * \param u A handle to the unbounded face.
   */
  virtual void after_clear() override {
    this->clear_landmark_set();
    this->build_landmark_set();
  }

  /*! Notification before a global operation modifies the arrangement.
   */
  virtual void before_global_change() override {
    this->clear_landmark_set();
    m_ignore_notifications = true;
  }

  /*! Notification after a global operation is completed.
   */
  virtual void after_global_change() override {
    this->build_landmark_set();
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

  /*! Notification after the creation of a new vertex. */
  virtual void after_create_vertex(Vertex_handle) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after the creation of a new edge. */
  virtual void after_create_edge(Halfedge_handle) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an edge was split. */
  virtual void after_split_edge(Halfedge_handle, Halfedge_handle) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after a face was split. */
  virtual void after_split_face(Face_handle, Face_handle, bool) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an outer CCB was split.*/
  virtual void after_split_outer_ccb(Face_handle,
                                     Ccb_halfedge_circulator,
                                     Ccb_halfedge_circulator) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an inner CCB was split. */
  virtual void after_split_inner_ccb(Face_handle,
                                     Ccb_halfedge_circulator,
                                     Ccb_halfedge_circulator) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an outer CCB was added to a face. */
  virtual void after_add_outer_ccb(Ccb_halfedge_circulator) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an inner CCB was created inside a face. */
  virtual void after_add_inner_ccb(Ccb_halfedge_circulator) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an isolated vertex was created inside a face. */
  virtual void after_add_isolated_vertex(Vertex_handle) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an edge was merged. */
  virtual void after_merge_edge(Halfedge_handle) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after a face was merged. */
  virtual void after_merge_face(Face_handle) override {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an outer CCB was merged. */
  virtual void after_merge_outer_ccb(Face_handle, Ccb_halfedge_circulator)
  override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an inner CCB was merged. */
  virtual void after_merge_inner_ccb(Face_handle, Ccb_halfedge_circulator)
  override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an outer CCB is moved from one face to another. */
  virtual void after_move_outer_ccb(Ccb_halfedge_circulator) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an inner CCB is moved from one face to another. */
  virtual void after_move_inner_ccb(Ccb_halfedge_circulator) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after an isolated vertex is moved. */
  virtual void after_move_isolated_vertex(Vertex_handle) override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after the removal of a vertex. */
  virtual void after_remove_vertex() override {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after the removal of an edge. */
  virtual void after_remove_edge() override {
    if (! m_ignore_notifications) {
      clear_landmark_set();
      build_landmark_set();
    }
    m_ignore_remove_edge = false;
  }

  /*! Notification after the removal of an outer CCB. */
  virtual void after_remove_outer_ccb(Face_handle) override {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_landmark_set();
      build_landmark_set();
    }
  }

  /*! Notification after the removal of an inner CCB. */
  virtual void after_remove_inner_ccb(Face_handle) override {
    if (! m_ignore_notifications && ! m_ignore_remove_edge) {
      clear_landmark_set();
      build_landmark_set();
    }
  }
  //@}

protected:
  /*! Create the list of landmarks with their location.
   * This is a pure virtual function, and the class that inherites from
   * this generator must implement it.
   */
  virtual void _create_points_set(Points_set&) = 0;

  virtual void _create_nn_points_set(NN_Points_set& nn_points) {
    Points_set points;
    Pairs_set pairs;

    // Create the set of landmark points.
    _create_points_set(points);

    // Locate the landmarks in the arrangement using batched point-location
    // global function.
    locate(*(this->arrangement()), points.begin(), points.end(),
           std::back_inserter(pairs));

    // Apply a random shuffle on the points, since the batched point-location
    // returns them sorted.
    CGAL::cpp98::random_shuffle(pairs.begin(), pairs.end());

    // Insert all landmarks (paired with their current location in the
    // arrangement) into the nearest-neighbor search structure.
    for (auto itr = pairs.begin(); itr != pairs.end(); ++itr) {
      NN_Point_2 np(itr->first, itr->second);
      nn_points.push_back(np);
    }
  }
};

} //namespace CGAL

#endif
