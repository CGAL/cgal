// Copyright (c) 2007  INRIA (France).
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
// Author(s)     : Laurent Saboret, Pierre Alliez

#ifndef GYROVIZ_POINT_3_H
#define GYROVIZ_POINT_3_H

#include <CGAL/Point_with_normal_3.h>

#include <vector>
#include <algorithm>


/// The Gyroviz_point_3 class represents a 3D point with:
/// - a position,
/// - a normal (oriented or not),
///  - a list of cameras used to reconstruct the point from an image sequence.
///
/// @heading Is Model for the Concepts: Model of the PointWithNormal_3 concept.
///
/// @heading Parameters:
/// @param Gt   Kernel's geometric traits.

template<class Gt>
class Gyroviz_point_3 : public CGAL::Point_with_normal_3<Gt>
{
// Private types
private:

  typedef typename CGAL::Point_with_normal_3<Gt> Base;

// Public types
public:

    // Repeat Point_with_normal_3 public types
    typedef Gt Geom_traits; ///< Kernel's geometric traits
	  typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3  Point;  ///< Kernel's Point_3 class.
    typedef typename CGAL::Oriented_normal_3<Geom_traits> Normal; ///< Model of OrientedNormal_3 concept.

    /// Iterator over cameras
    typedef typename std::vector<Point>::const_iterator Camera_const_iterator;

// Public methods
public:

		/// Point is (0,0,0) by default.
		/// Normal is (0,0,0) by default.
		/// Normal is oriented by default.
		/// Camera list is empty by default.
    Gyroviz_point_3(const CGAL::Origin& o = CGAL::ORIGIN)
    : Base(o)
		{
		}
    Gyroviz_point_3(FT x, FT y, FT z)
    : Base(x,y,z)
		{
		}
    Gyroviz_point_3(const Point& point,
                    const Normal& normal = CGAL::NULL_VECTOR)
    : Base(point, normal)
		{
		}
    template < class InputIterator >
    Gyroviz_point_3(const Point& point,
	                  const Normal& normal,
										InputIterator first_camera, InputIterator beyond_camera)
    : Base(point, normal)
    {
      std::copy(first_camera, beyond_camera, std::back_inserter(cameras));
    }

    // Default copy constructor and operator =() are fine

    /// Compare positions
    bool operator==(const Gyroviz_point_3& that)
    { 
      return Base::operator==(that); 
    }
    bool operator!=(const Gyroviz_point_3& that)
    { 
      return ! (*this == that); 
    }

    // Get/set cameras
    template < class InputIterator >
    void set_cameras(InputIterator first_camera, InputIterator beyond_camera)
    {
      cameras.clear();
      std::copy(first_camera, beyond_camera, std::back_inserter(cameras));
    }
    Camera_const_iterator cameras_begin() const { return	cameras.begin(); }
    Camera_const_iterator cameras_end  () const { return	cameras.end(); }
    
// Data
private:

  // List of cameras
  std::vector<Point> cameras;
};


#endif //GYROVIZ_POINT_3_H

