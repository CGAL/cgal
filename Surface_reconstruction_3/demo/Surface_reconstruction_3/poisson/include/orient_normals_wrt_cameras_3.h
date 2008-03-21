// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
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
// Author(s) : Pierre Alliez and Laurent Saboret


#ifndef ORIENT_NORMALS_WRT_CAMERAS_3_H
#define ORIENT_NORMALS_WRT_CAMERAS_3_H

#include <CGAL/basic.h>
#include <Gyroviz_point_3.h>

#include <iterator>


/// Orient a 3D point's normal w.r.t. the position of cameras
/// that reconstructed the point by photogrammetry.
/// @return true on success.
template < class Gt, class OrientedNormal_3, class InputIterator >
void
orient_normal_wrt_cameras_3(const typename Gt::Point_3& p, // in
			                      OrientedNormal_3& normal, // in and out
										        InputIterator first_camera, InputIterator beyond_camera) // in
{
  typedef typename Gt::FT       FT;
  typedef typename Gt::Point_3  Point;
  typedef OrientedNormal_3      Normal;
  typedef typename Gt::Vector_3 Vector;

  Vector n = normal.get_vector();

  //                                        ->
	//                                        cp     ->
	// Search for the camera position c / | ------ * n | is maximum
	//                                      ||cp||
	FT max_dot_product = 0;
	Point max_camera;
	for(InputIterator c = first_camera; c != beyond_camera; c++)
	{
		Vector cp = p - *c;
		FT norm_cp = std::sqrt(cp * cp);
		if (norm_cp != 0.0)
		{
      FT dot = (cp * n) / norm_cp;
		  if (dot < 0)
		    dot = -dot;

		  if (max_dot_product < dot)
		  {
		     max_dot_product = dot;
		     max_camera = *c;
		  }
		}
	}

	//        ->             ->
	// Orient n backwards to cp
	if (max_dot_product > 0)
	{
		Vector cp = p - max_camera;
		FT dot = (cp * n);
		if (dot > 0)
      n = -n;
	  normal = Normal(n, true /* oriented */);
	}
	else // if failure
	{
	  normal = Normal(n, false /* non oriented */);
	}
    //// TEST: flag only inverted normals as oriented to see the result in 3D rendering
    //if (max_dot_product > 0)
    //{
	   // Vector cp = p - max_camera;
	   // FT dot = (cp * n);
	   // if (dot > 0)
	   //   normal = Normal(-n, true /* oriented */);
    //  else
  	 //   normal = Normal(n, false /* non oriented */);
    //}
    //else // if failure
    //{
	   // normal = Normal(n, false /* non oriented */);
    //}
}


/// Specific to Gyroviz: orient the normals w.r.t. the position of cameras
/// that reconstructed the points by photogrammetry.
///
/// Preconditions:
/// - ForwardIterator value_type is Gyroviz_point_3.

template < class ForwardIterator >
void
orient_normals_wrt_cameras_3(ForwardIterator first, ForwardIterator beyond)
{
	typedef typename std::iterator_traits<ForwardIterator>::value_type Value_type;
  typedef typename Value_type::Geom_traits Gt;

	// iterate over input points and orient normals
	for (ForwardIterator it = first; it != beyond; it++)
		orient_normal_wrt_cameras_3<Gt>(*it, it->normal(), it->cameras_begin(), it->cameras_end());
}


#endif // ORIENT_NORMALS_WRT_CAMERAS_3_H

