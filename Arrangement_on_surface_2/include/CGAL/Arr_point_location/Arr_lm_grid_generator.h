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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
//                 Ron Wein     <wein@post.tau.ac.il>

#ifndef CGAL_ARR_LANDMARKS_GRID_GENERATOR_H
#define CGAL_ARR_LANDMARKS_GRID_GENERATOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
* Definition of the Arr_grid_landmarks_generator<Arrangement> template.
*/

#include <CGAL/Arr_point_location/Arr_lm_generator_base.h>

namespace CGAL {

/*! \class Arr_grid_landmarks_generator
 * A generator for the landmarks point-locatoion class, which uses a
 * set of points on a grid as its set of landmarks.
*/
template <typename Arrangement_,
          typename Nearest_neighbor_ =
            Arr_landmarks_nearest_neighbor<Arrangement_> >
class Arr_grid_landmarks_generator :
    public Arr_landmarks_generator_base<Arrangement_, Nearest_neighbor_>
{
public:
  typedef Arrangement_                                  Arrangement_2;
  typedef Nearest_neighbor_                             Nearest_neighbor;

  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Ccb_halfedge_circulator
                                                        Ccb_halfedge_circulator;

  typedef typename Geometry_traits_2::Approximate_number_type
                                                        ANT;

  typedef typename Arrangement_2::Point_2               Point_2;

private:
  typedef Arr_landmarks_generator_base<Arrangement_2, Nearest_neighbor>
                                                        Base;
  typedef Arr_grid_landmarks_generator<Arrangement_2, Nearest_neighbor>
                                                        Self;

protected:
  typedef typename Base::Points_set                     Points_set;
  typedef typename Base::PL_result_type                 PL_result_type;
  typedef std::pair<Point_2, PL_result_type>            PL_pair;
  typedef std::vector<PL_pair>                          Pairs_set;

  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2> Traits_adaptor_2;

  // Data members:
  const Traits_adaptor_2* m_traits;
  unsigned int num_landmarks;
  Pairs_set lm_pairs;

  ANT x_min, y_min;    // Bounding box for the
  ANT x_max, y_max;    // arrangement vertices.
  ANT step_x, step_y;  // Grid step sizes.
  unsigned int sqrt_n;

  bool fixed_number_of_lm; // indicates if the constructor got
                           // number of landmarks as parameter

private:
  /*! Copy constructor - not supported. */
  Arr_grid_landmarks_generator(const Self&);

  /*! Assignment operator - not supported. */
  Self& operator=(const Self&);

public:
  /*! Constructor from an arrangement.
   * \param arr (in) The arrangement.
   */
  Arr_grid_landmarks_generator(const Arrangement_2& arr) :
    Base(arr),
    m_traits(static_cast<const Traits_adaptor_2*>(arr.geometry_traits())),
    num_landmarks(0),
    fixed_number_of_lm(false)
  {
    build_landmark_set();//this->
  }

  Arr_grid_landmarks_generator(const Arrangement_2& arr,
                               unsigned int n_landmarks) :
    Base(arr),
    m_traits(static_cast<const Traits_adaptor_2*>(arr.geometry_traits())),
    num_landmarks(n_landmarks),
    fixed_number_of_lm(true)
  {
    build_landmark_set();//this->
  }

  /*! Create the landmarks set (choosing the landmarks),
   * and store them in the nearest neighbor search structure.
   */
  virtual void build_landmark_set()
  {
    // Create a set of points on a grid.
    Points_set points;
    _create_points_set(points);
    // Locate the landmarks in the arrangement using batched point-location
    // global function. Note that the resulting pairs are returned sorted by
    // their lexicographic xy-order.
    lm_pairs.clear();
    locate(*(this->arrangement()), points.begin(), points.end(),
           std::back_inserter(lm_pairs));
    this->updated = true;
  }

  /*! Clear the set of landmarks.
   */
  virtual void clear_landmark_set()
  {
    lm_pairs.clear();
    this->updated = false;
  }

  /*! Obtain the nearest neighbor (landmark) to the given point.
   * \param q The query point.
   * \param obj (out) The location of the nearest landmark point in the
   *                  arrangement (a vertex, halfedge, or face handle).
   * \return The nearest landmark point.
   */
  virtual Point_2 closest_landmark(const Point_2& q, PL_result_type& obj)
  {
    CGAL_assertion(this->updated);

    // Calculate the index of the nearest grid point point to q.
    typename Geometry_traits_2::Approximate_2 approximate =
      m_traits->approximate_2_object();
    const ANT qx = approximate(q, 0);
    const ANT qy = approximate(q, 1);
    unsigned int i = (CGAL::compare(qx, x_min) == SMALLER) ? 0 :
      (CGAL::compare(qx, x_max) == LARGER) ? (sqrt_n - 1) :
      static_cast<int>(((qx - x_min) / step_x) + 0.5);
    unsigned int j = (CGAL::compare(qy, y_min) == SMALLER) ? 0 :
      (CGAL::compare(qy, y_max) == LARGER) ? (sqrt_n - 1) :
      static_cast<int>(((qy - y_min) / step_y) + 0.5);
    unsigned int index = sqrt_n * i + j;

    // Return the result.
    obj = lm_pairs[index].second;
    return (lm_pairs[index].first);
  }

protected:
  /*! Create a set of landmark points on a grid.
   */
  virtual void _create_points_set(Points_set& points)
  {
    Arrangement_2* arr = this->arrangement();

    if (arr->is_empty())
      return;

    // Locate the arrangement vertices with minimal and maximal x and
    // y-coordinates.
    Vertex_const_iterator    vit = arr->vertices_begin();
    x_min = x_max = m_traits->approximate_2_object()(vit->point(), 0);
    y_min = y_max = m_traits->approximate_2_object()(vit->point(), 1);

    if (arr->number_of_vertices() == 1) {
      // There is only one isolated vertex at the arrangement:
      step_x = step_y = 1;
      sqrt_n = 1;
      points.push_back(Point_2(x_min, y_min));
      return;
    }

    ANT x, y;
    Vertex_const_iterator left, right, top, bottom;

    left = right = top = bottom = vit;

    for (++vit; vit != arr->vertices_end(); ++vit) {
      x = m_traits->approximate_2_object()(vit->point(), 0);
      y = m_traits->approximate_2_object()(vit->point(), 1);

      if (CGAL::compare(x, x_min) == SMALLER) {
        x_min = x;
        left = vit;
      }
      else if (CGAL::compare(x, x_max) == LARGER) {
        x_max = x;
        right = vit;
      }

      if (CGAL::compare(y, y_min) == SMALLER) {
        y_min = y;
        bottom = vit;
      }
      else if (CGAL::compare(y, y_max) == LARGER) {
        y_max = y;
        top = vit;
      }
    }

    // Create N Halton points. If N was not given to the constructor,
    // set it to be the number of vertices V in the arrangement (actually
    // we generate ceiling(sqrt(V))^2 landmarks to obtain a square grid).
    if (!fixed_number_of_lm)
      num_landmarks = static_cast<unsigned int>(arr->number_of_vertices());

    sqrt_n = static_cast<unsigned int>
      (std::sqrt(static_cast<double>(num_landmarks)) + 0.99999);
    num_landmarks = sqrt_n * sqrt_n;

    CGAL_assertion(sqrt_n > 1);

    // Calculate the step sizes for the grid.
    ANT delta_x = m_traits->approximate_2_object()(right->point(), 0) -
      m_traits->approximate_2_object()(left->point(), 0);
    ANT delta_y = m_traits->approximate_2_object()(top->point(), 1) -
      m_traits->approximate_2_object()(bottom->point(), 1);

    if (CGAL::sign(delta_x) == CGAL::ZERO) delta_x = delta_y;
    if (CGAL::sign(delta_y) == CGAL::ZERO) delta_y = delta_x;
    CGAL_assertion((CGAL::sign(delta_x) == CGAL::POSITIVE) &&
                   (CGAL::sign(delta_y) == CGAL::POSITIVE));
    step_x = delta_x / (sqrt_n - 1);
    step_y = delta_y / (sqrt_n - 1);

    // Create the points on the grid.
    const double  x_min =
      CGAL::to_double(m_traits->approximate_2_object()(left->point(), 0));
    const double  y_min =
      CGAL::to_double(m_traits->approximate_2_object()(bottom->point(), 1));
    const double  sx = CGAL::to_double(step_x);
    const double  sy = CGAL::to_double(step_y);
    double        px, py;
    unsigned int  i, j;

    for (i = 0; i< sqrt_n; i++) {
      px = x_min + i*sx;
      for (j = 0; j< sqrt_n; j++) {
        py = y_min + j*sy;
        points.push_back(Point_2(px, py));
      }
    }
  }
};

} //namespace CGAL

#endif
