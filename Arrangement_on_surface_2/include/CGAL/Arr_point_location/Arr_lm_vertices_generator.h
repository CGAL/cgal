// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
//                 Ron Wein     <wein@post.tau.ac.il>
#ifndef CGAL_ARR_LANDMARKS_VERTICES_GENERATOR_H
#define CGAL_ARR_LANDMARKS_VERTICES_GENERATOR_H

/*! \file
 * Definition of the Arr_landmarks_vertices_generator<Arrangement> template.
 */

#include <CGAL/Arr_point_location/Arr_lm_generator_base.h>

namespace CGAL {

/*! \class Arr_landmarks_vertices_generator
 * A generator for the landmarks point-locatoion class, which uses the
 * arrangement vertices as its set of landmarks.
*/
template <class Arrangement_,
          class Nearest_neighbor_  =
            Arr_landmarks_nearest_neighbor<typename
                                           Arrangement_::Geometry_traits_2> >
class Arr_landmarks_vertices_generator :
    public Arr_landmarks_generator_base<Arrangement_, Nearest_neighbor_>
{
public:

  typedef Arrangement_                                  Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef Nearest_neighbor_                             Nearest_neighbor;

  typedef Arr_landmarks_generator_base<Arrangement_2,
                                       Nearest_neighbor>      Base;
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

  typedef typename Base::Points_set                     Points_set;

  typedef typename Nearest_neighbor::NN_Point_2         NN_Point_2;
  typedef std::list<NN_Point_2>                         NN_Point_list;

protected:

  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2> Traits_adaptor_2;

  // Data members:
  const Traits_adaptor_2  *m_traits;  // Its associated traits object.
  int                      num_landmarks;

private:

  /*! Copy constructor - not supported. */
  Arr_landmarks_vertices_generator (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

public: 

  /*! Constructor. */
  Arr_landmarks_vertices_generator (const Arrangement_2& arr) :
    Base (arr),
    num_landmarks(0)
  {
    m_traits = static_cast<const Traits_adaptor_2*> (arr.geometry_traits());
    build_landmark_set();//this->
  }

  virtual void _create_points_set (Points_set & /* points */)
  {
    std::cerr << "should not reach here!"<< std::endl;
    CGAL_error();
  }

  /*!
   * Creates the landmark set, using all arrangement vertices.
   */
  virtual void build_landmark_set ()
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
    this->nn.clear();
    this->nn.init (nnp_list.begin(), nnp_list.end());

    this->num_small_not_updated_changes = 0;
    this->updated = true;
  }

  /*!
   * Clear the landmark set.
   */
  virtual void clear_landmark_set ()
  {
    this->nn.clear();
    num_landmarks = 0;
    this->num_small_not_updated_changes = 0;
    this->updated = false;
  }

protected:

  /*! Handle a local change. */
  void _handle_local_change_notification ()
  {
    // Rebuild the landmark set only if the number of small
    // changes is greater than sqrt(num_landmarks).
    double    nl = static_cast<double> (num_landmarks);
    const int sqrt_num_landmarks = 
      static_cast<int> (std::sqrt (nl) + 0.5);

    this->num_small_not_updated_changes++;
    if ((num_landmarks < 10) ||
        (this->num_small_not_updated_changes >=  sqrt_num_landmarks))
    {
      clear_landmark_set();
      build_landmark_set();//this->
    }

    return;
  }

};

} //namespace CGAL

#endif
