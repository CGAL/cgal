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
//                 Ron Wein     <haranidi@post.tau.ac.il>
#ifndef CGAL_ARR_LM_HALTON_GENERATOR_H
#define CGAL_ARR_LM_HALTON_GENERATOR_H

/*! \file
* Definition of the Arr_halton_landmarks_generator<Arrangement> template.
*/

#include <CGAL/Arr_point_location/Arr_lm_generator_base.h>

namespace CGAL {

/*! \class
 * A generator for the landmarks point-locatoion class, which uses a
 * sequence of Halton points as its set of landmarks.
 */
template <class Arrangement_,
          class Nearest_neighbor_  =
            Arr_landmarks_nearest_neighbor<typename
                                           Arrangement_::Geometry_traits_2> >
class Arr_halton_landmarks_generator :
    public Arr_landmarks_generator_base<Arrangement_, Nearest_neighbor_>
{
public:
  typedef Arrangement_                                  Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef Nearest_neighbor_                             Nearest_neighbor;

  typedef Arr_landmarks_generator_base<Arrangement_2,
                                       Nearest_neighbor>    Base;
  typedef Arr_halton_landmarks_generator<Arrangement_2,
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
  Arr_halton_landmarks_generator (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );


public: 
  
  /*! Constructor. */
  Arr_halton_landmarks_generator (const Arrangement_2& arr,
                                  unsigned int n_landmarks = 0) : 
    Base (arr),
    num_landmarks (n_landmarks)
  {
    this->build_landmark_set();
  }

protected:
  
  /*!
   * Create a set of Halton points (the number of points is given as a
   * parameter to the constructor, or is taken from the arrangement size).
   * The Halton points are constructed in the bounding rectangle of the
   * arrangement vertices.
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

    // Create N Halton points. If N was not given to the constructor,
    // set it to be the number of vertices in the arrangement.
    if (num_landmarks == 0)
      num_landmarks = static_cast<unsigned int>(arr->number_of_vertices());

    if (num_landmarks == 0)
      return;

    if (num_landmarks == 1)
    {
      points.push_back (Point_2 (x_max, y_max)); 
      return;
    }

    // Create the Halton sequence.
    const double  x_scale = x_max - x_min;
    const double  y_scale = y_max - y_min;
    double        base_inv;
    int           digit, i, seed2;
    int           base[2], leap[2], seed[2];
    double        r[2];
    double        px, py;
    int           ndim = 2;
    unsigned int  step = 1;

    seed[0] = seed[1] = 0;
    leap[0] = leap[1] = 1;
    base[0] = 2;
    base[1] = 3;
    for (step = 1; step <= num_landmarks; step++)
    {
      for (i = 0; i < ndim; i++)
      {
        seed2 = seed[i] + step * leap[i];
        r[i] = 0;
        base_inv = 1 / static_cast<double> (base[i]);
        while (seed2 != 0)
        {
          digit = seed2 % base[i];
          r[i] = r[i] + static_cast<double> (digit) * base_inv;
          base_inv = base_inv / static_cast<double> (base[i]);
          seed2 = seed2 / base[i];
        }
      }

      // Now r[0] and r[1] are the x and y-coordinates of the Halton point
      // in the unit square. We scale and shift them to fit out bounding box.
      px = r[0] * x_scale + x_min;
      py = r[1] * y_scale + y_min;

      points.push_back (Point_2 (px, py));
    }
    return;
  }

};

} //namespace CGAL

#endif
