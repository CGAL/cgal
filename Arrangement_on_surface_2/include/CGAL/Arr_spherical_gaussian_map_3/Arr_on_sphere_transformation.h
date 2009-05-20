// Copyright (c) 2009 Tel-Aviv University (Israel).
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
// $URL: svn+ssh://naamamay@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_on_surface_2/include/CGAL/Arr_spherical_gaussian_map_3/Arr_on_sphere_transformation $
// $Id: Arr_on_sphere_transformation.h 43873 2008-06-30 00:07:53Z naamamay $
// 
// Author(s)     : Naama mayer         <naamamay@post.tau.ac.il>

#ifndef CGAL_ARR_ON_SPHERE_TRANSFORMATION_H
#define CGAL_ARR_ON_SPHERE_TRANSFORMATION_H

CGAL_BEGIN_NAMESPACE

/* This function rotates the face when the arrangement is sgm */

template <class Arrangement, class Transformation_3>
class Arr_on_sphere_transformation
{
	void rotate_face(Arrangement & arr, Face_Handle f1, Face_handle f2)
	{			

	}
};

CGAL_END_NAMESPACE

#endif