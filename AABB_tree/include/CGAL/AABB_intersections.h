// Copyright (c) 2008-2009  INRIA Sophia-Antipolis (France), ETHZ (Suisse).
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
// $URL: svn+ssh://lrineau@scm.gforge.inria.fr/svn/cgal/trunk/AABB_tree/include/CGAL/AABB_intersections/Bbox_3_ray_3_do_intersect.h $
// $Id: Bbox_3_ray_3_do_intersect.h 50124 2009-06-26 17:55:17Z palliez $
// 
//
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez

#include <CGAL/intersections.h>

#include <CGAL/AABB_intersections/Bbox_3_ray_3_do_intersect.h>
#include <CGAL/AABB_intersections/Bbox_3_line_3_do_intersect.h>
#include <CGAL/AABB_intersections/Bbox_3_bbox_3_do_intersect.h>
#include <CGAL/AABB_intersections/Bbox_3_plane_3_do_intersect.h>
#include <CGAL/AABB_intersections/Bbox_3_sphere_3_do_intersect.h>
#include <CGAL/AABB_intersections/Bbox_3_segment_3_do_intersect.h>
#include <CGAL/AABB_intersections/Bbox_3_triangle_3_do_intersect.h>

#include <CGAL/AABB_intersections/Triangle_3_ray_3_intersection.h>
#include <CGAL/AABB_intersections/Triangle_3_line_3_intersection.h>
#include <CGAL/AABB_intersections/Triangle_3_plane_3_intersection.h>
#include <CGAL/AABB_intersections/Triangle_3_segment_3_intersection.h>

#include <CGAL/AABB_intersections/nearest_point_segment_3.h>
#include <CGAL/AABB_intersections/nearest_point_triangle_3.h>
