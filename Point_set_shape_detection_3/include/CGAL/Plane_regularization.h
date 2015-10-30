// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : 
//

/**
* \ingroup PkgPointSetShapeDetection3
* \file CGAL/Plane_regularization.h
*
*/


#ifndef CGAL_PLANE_REGULARIZATION_H
#define CGAL_PLANE_REGULARIZATION_H

#include <CGAL/Shape_detection_3.h>

#include <boost/foreach.hpp>


namespace CGAL {

template <typename Traits>
class Plane_regularization
{
public:

  typedef Plane_regularization<Traits> Self;

  typedef typename Traits::FT FT;
  typedef typename Traits::Point_3 Point;
  typedef typename Traits::Vector_3 Vector;
  typedef typename Traits::Plane_3 Plane;

  typedef typename Traits::Point_map Point_map;
  typedef typename Traits::Normal_map Normal_map;
  typedef typename Traits::Input_range Input_range;
  typedef typename Input_range::iterator Input_iterator;

  typedef Shape_detection_3::Shape_base<Traits> Shape;
  typedef Shape_detection_3::Plane<Traits> Plane_shape;
  

private:

  Traits m_traits;

  Input_iterator m_input_begin;
  Input_iterator m_input_end;
  Point_map m_point_pmap;
  Normal_map m_normal_pmap;

  std::vector<boost::shared_ptr<Plane_shape> > m_planes;

public:

  Plane_regularization (Traits t = Traits ())
    : m_traits (t)
  {

  }

  Plane_regularization (Input_range& input_range,
                         Shape_detection_3::Efficient_RANSAC<Traits>& shape_detection)
    : m_traits (shape_detection.traits())
  {
    m_input_begin = input_range.begin ();
    m_input_end = input_range.end ();

    BOOST_FOREACH (boost::shared_ptr<Shape> shape, shape_detection.shapes())
      {
        boost::shared_ptr<Plane_shape> pshape
          = boost::dynamic_pointer_cast<Plane_shape>(shape);
        
        // Ignore all shapes other than plane
        if (!pshape)
          continue;

        m_planes.push_back (pshape);

      }

  }

  virtual ~Plane_regularization ()
  {
    clear ();
  }

  void clear ()
  {
    std::vector<boost::shared_ptr<Plane_shape> > ().swap (m_planes);
  }


  
  void run ()
  {

  }

  Vector regularize_normal (const Vector& n, FT cos_vertical)
  {
    FT A = 1 - cos_vertical * cos_vertical;
    FT B = 1 + (n.y() * n.y()) / (n.x() * n.x());
    
    FT vx = std::sqrt (A/B);
    
    if (n.x() < 0)
      vx = -vx;
    
    FT vy = vx * (n.y() / n.x()); 

    Vector res (vx, vy, cos_vertical);

    return res / std::sqrt (res * res);
  }

  Vector regularize_normals_from_prior (const Vector& np,
                                        const Vector& n,
                                        FT cos_vertical)
  {
    FT vx, vy;

    if (np.x() != 0)
      { 
        FT a = (np.y() * np.y()) / (np.x() * np.x()) + 1;
        FT b = 2 * np.y() * np.z() * cos_vertical / np.x();
        FT c= cos_vertical * cos_vertical-1;

        if (4 * a * c > b * b)
          return regularize_normal (n, cos_vertical); 
        else
          {
            FT delta = std::sqrt (b * b-4 * a * c);
            FT vy1= (-b-delta) / (2 * a);
            FT vy2= (-b+delta) / (2 * a);

            vy = (std::fabs(n.y()-vy1) < std::fabs(n.y()-vy2))
              ? vy1 : vy2;

            vx = (-np.y() * vy-np.z() * cos_vertical) / np.x();
          }
      }
    else if (np.y() != 0)
      {
        vy = -np.z() * cos_vertical / np.y();
        vx = std::sqrt (1 - cos_vertical * cos_vertical - vy * vy);
        
        if (n.x() < 0)
          vx = -vx;
      }
    else
      return regularize_normal (n, cos_vertical); 

    Vector_3 res (vx, vy, cos_vertical);
    FT norm = std::max(1e-5, 1. / sqrt(res.squared_length ()));

    return norm * res;
  }

  

};


}; // namespace CGAL

#endif // CGAL_PLANE_REGULARIZATION_H
