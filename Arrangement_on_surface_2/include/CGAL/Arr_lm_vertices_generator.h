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
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
//                 Ron Wein     <wein@post.tau.ac.il>
#ifndef CGAL_ARR_LANDMARKS_VERTICES_GENERATOR_H
#define CGAL_ARR_LANDMARKS_VERTICES_GENERATOR_H

/*! \file
 * Definition of the Arr_landmarks_vertices_generator<Arrangement> template.
 */

#include <list>
#include <CGAL/Arr_observer.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_point_location/Arr_lm_nearest_neighbor.h>

CGAL_BEGIN_NAMESPACE

/*! \class Arr_landmarks_vertices_generator
 * A generator for the landmarks point-locatoion class, which uses the
 * arrangement vertices as its set of landmarks.
*/
template <class Arrangement_,
          class Nearest_neighbor_  =
            Arr_landmarks_nearest_neighbor<typename
                                           Arrangement_::Geometry_traits_2> >
class Arr_landmarks_vertices_generator :
    public Arr_observer<Arrangement_>
{
public:

  typedef Arrangement_                                  Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef Nearest_neighbor_                             Nearest_neighbor;

  typedef Arr_landmarks_vertices_generator<Arrangement_2,
                                           Nearest_neighbor>  Self;
  
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;
  
  typedef typename Arrangement_2::Point_2               Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2    X_monotone_curve_2;

  typedef typename Nearest_neighbor::NN_Point_2         NN_Point_2;
  typedef std::list<NN_Point_2>                         NN_Point_list;

protected:

  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2> Traits_adaptor_2;

  // Data members:
  const Traits_adaptor_2  *traits;  // Its associated traits object.
  Nearest_neighbor         nn;      // The associated nearest neighbor object.
  bool                     ignore_notifications;
  bool                     updated;
  int                      num_small_not_updated_changes;
  int                      num_landmarks;

private:

  /*! Copy constructor - not supported. */
  Arr_landmarks_vertices_generator (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

public: 

  /*! Constructor. */
  Arr_landmarks_vertices_generator (const Arrangement_2& arr) :
    Arr_observer<Arrangement_2> (const_cast<Arrangement_2 &>(arr)),
    ignore_notifications (false),
    updated (false),
    num_small_not_updated_changes(0),
    num_landmarks(0)
  {
    traits = static_cast<const Traits_adaptor_2*> (arr.geometry_traits());
    build_landmark_set();
  }

  /*!
   * Creates the landmark set, using all arrangement vertices.
   */
  void build_landmark_set ()
  {
    // Go over the arrangement, and insert all its vertices as landmarks.
    NN_Point_list           nnp_list; 
    const Arrangement_2    *arr = this->arrangement();
    Vertex_const_iterator   vit;
    Vertex_const_handle     vh;

    num_landmarks = 0;
    for (vit = arr->vertices_begin(); vit != arr->vertices_end(); ++vit)
    {
      vh = vit;
      nnp_list.push_back (NN_Point_2 (vh->point(),
                                      CGAL::make_object (vh)));
      num_landmarks++;
    }

    // Update the search structure.
    nn.clear();
    nn.init (nnp_list.begin(), nnp_list.end());

    num_small_not_updated_changes = 0;
    updated = true;
  }

  /*!
   * Clear the landmark set.
   */
  void clear_landmark_set ()
  {
    nn.clear();

    num_landmarks = 0;
    num_small_not_updated_changes = 0;
    updated = false;
  }

  /*!
   * Get the nearest neighbor (landmark) to the given point.
   * \param q The query point.
   * \param obj Output: The location of the nearest landmark point in the
   *                    arrangement (a vertex, halfedge, or face handle).
   * \return The nearest landmark point.
   */
  const Point_2& get_closest_landmark (const Point_2& q, Object &obj)
  {
    CGAL_assertion(updated);
    return (nn.find_nearest_neighbor (q, obj));
  }
  
  /// \name Overloaded observer functions.
  //@{
  
  /*!
   * Notification after the arrangement has been assigned with another
   * arrangement.
   */
  virtual void after_assign ()
  { 
    clear_landmark_set();
    build_landmark_set();
    ignore_notifications = false;
  }
  
  /*! 
   * Notification before the observer is attached to an arrangement.
   * \param arr The arrangement we are about to attach the observer to.
   */
  virtual void before_attach (const Arrangement_2& arr)
  {
    clear_landmark_set();
    traits = static_cast<const Traits_adaptor_2*> (arr.geometry_traits());
    ignore_notifications = false;
  }
  
  /*!
   * Notification after the observer has been attached to an arrangement.
   */
  virtual void after_attach ()
  {
    build_landmark_set();
  }

  /*! 
   * Notification before the observer is detached from the arrangement.
   */
  virtual void before_detach ()
  {
    clear_landmark_set();
  }

  /*!
   * Notification after the arrangement is cleared.
   */
  virtual void after_clear ()
  { 
    clear_landmark_set();
    build_landmark_set();
  }

  /*! Notification before a global operation modifies the arrangement. */
  virtual void before_global_change ()
  {
    clear_landmark_set();
    ignore_notifications = true;
  }

  /*! Notification after a global operation is completed. */
  virtual void after_global_change ()
  {
    build_landmark_set();
    ignore_notifications = false;
  }

  /*! Notification after the creation of a new vertex. */
  virtual void after_create_vertex (Vertex_handle )
  {
    if (! ignore_notifications)
      _handle_local_change();
  }
  
  /*! Notificaion after the removal of a vertex. */
  virtual void after_remove_vertex ()
  {
    if (! ignore_notifications)
      _handle_local_change();
  }
  //@}

protected:

  /*! Handle a local change. */
  void _handle_local_change ()
  {
    // Rebuild the landmark set only if the number of small
    // changes is greater than sqrt(num_landmarks).
    double    nl = static_cast<double> (num_landmarks);
    const int sqrt_num_landmarks = 
      static_cast<int> (std::sqrt (nl) + 0.5);

    num_small_not_updated_changes++;
    if ((num_landmarks < 10) ||
        (num_small_not_updated_changes >=  sqrt_num_landmarks))
    {
      clear_landmark_set();
      build_landmark_set();
    }

    return;
  }

};

CGAL_END_NAMESPACE

#endif
