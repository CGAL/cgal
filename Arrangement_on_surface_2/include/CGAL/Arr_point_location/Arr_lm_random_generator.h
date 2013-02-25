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
#ifndef CGAL_ARR_LM_RANDOM_GENERATOR_H
#define CGAL_ARR_LM_RANDOM_GENERATOR_H

/*! \file
* Definition of the Arr_random_landmarks_generator<Arrangement> template.
*/

#include <CGAL/Arr_point_location/Arr_lm_generator_base.h>
#include <CGAL/Random.h>

namespace CGAL {

/*! \class Arr_random_landmarks_generator
 * A generator for the landmarks point-locatoion class, which uses a
 * random set of points as its set of landmarks.
*/

template <class Arrangement_,
          class Nearest_neighbor_  =
            Arr_landmarks_nearest_neighbor<typename
                                           Arrangement_::Geometry_traits_2> >
class Arr_random_landmarks_generator :
    public Arr_landmarks_generator_base<Arrangement_, Nearest_neighbor_>
{
public:
  typedef Arrangement_                                  Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef Nearest_neighbor_                             Nearest_neighbor;

  typedef Arr_landmarks_generator_base<Arrangement_2,
                                       Nearest_neighbor>    Base;
  typedef Arr_random_landmarks_generator<Arrangement_2,
                                         Nearest_neighbor>  Self;

  typedef typename Arrangement_2::Point_2               Point_2;
  typedef typename Base::Points_set                     Points_set;

  typedef typename Arrangement_2::Vertex_const_iterator
                                                Vertex_const_iterator;

protected:

  // Data members:
  unsigned int     num_landmarks; 

private:

  /*! Copy constructor - not supported. */
  Arr_random_landmarks_generator (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

  
public: 

  /*! Constructor. */
  Arr_random_landmarks_generator (const Arrangement_2& arr,
                                  unsigned int n_landmarks = 0) :
    Base (arr),
    num_landmarks (n_landmarks)
  {
    this->build_landmark_set();
  }

protected:

  /*!
   * Create a set of random points (the number of points is given as a
   * parameter to the constructor, or is taken from the arrangement size).
   * The coordinates of the landmarks are selected randomly in the
   * bounding rectangle of the Arrangement's vertices 
   */
  virtual void _create_points_set (Points_set& points)
  {
    points.clear();

    // Go over the arrangement vertices and construct their boundig box.
    const Arrangement_2    *arr = this->arrangement();
    Vertex_const_iterator   vit; 
    double                  x_min = 0, x_max = 1, y_min = 0, y_max = 1;
    double                  x, y;
    bool                    first = true;

    for (vit=arr->vertices_begin(); vit != arr->vertices_end(); ++vit)
    {
      x = CGAL::to_double(vit->point().x());
      y = CGAL::to_double(vit->point().y());

      if (first)
      {
        x_min = x_max = x;
        y_min = y_max = y;
        first = false;
      }
      else
      {
        if (x < x_min)
          x_min = x;
        else if (x > x_max)
          x_max = x;

        if (y < y_min)
          y_min = y;
        else if (y > y_max)
          y_max = y;
      }
    }

    // Create N random landmarks. If N was not given to the constructor,
    // set it to be the number of vertices in the arrangement.
    if (num_landmarks == 0)
      num_landmarks = static_cast<unsigned int>(arr->number_of_vertices());

    CGAL::Random random;
    for (unsigned int i = 0; i < num_landmarks; ++i) {
      double px = (x_min == x_max) ? x_min : random.get_double(x_min, x_max);
      double py = (y_min == y_max) ? y_min : random.get_double(y_min, y_max);
      points.push_back(Point_2 (px, py)); 
    }
  }

};

} //namespace CGAL


#endif
