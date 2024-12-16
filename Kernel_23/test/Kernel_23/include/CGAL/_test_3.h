// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Stefan Schirra


#ifndef CGAL__TEST_3_C
#define CGAL__TEST_3_C

#include "_test_cls_vector_3.h"
#include "_test_fct_vector_3.h"
#include "_test_cls_point_3.h"
#include "_test_cls_weighted_point_3.h"
#include "_test_fct_point_vector_3.h"
#include "_test_fct_point_3.h"
#include "_test_fct_weighted_point_3.h"
#include "_test_fct_plane_3.h"
#include "_test_further_fct_point_plane_3.h"
#include "_test_cls_direction_3.h"
#include "_test_cls_plane_3.h"
#include "_test_cls_line_3.h"
#include "_test_cls_segment_3.h"
#include "_test_fct_segment_3.h"
#include "_test_cls_sphere_3.h"
#include "_test_cls_ray_3.h"
#include "_test_cls_triangle_3.h"
#include "_test_cls_tetrahedron_3.h"
#include "_test_cls_aff_transformation_3.h"
#include "_test_cls_iso_cuboid_3.h"

template <class R>
bool
_test_3(const R& r)
{
 return
    _test_cls_vector_3(r)
 && _test_fct_vector_3(r)
 && _test_cls_point_3(r)
 && _test_cls_weighted_point_3(r)
 && _test_fct_point_vector_3(r)
 && _test_fct_point_3(r)
 && _test_fct_weighted_point_3(r)
 && _test_fct_plane_3(r)
 && _test_fct_segment_3(r)
 && _test_further_fct_point_plane_3(r)
 && _test_cls_direction_3(r)
 && _test_cls_plane_3( r )
 && _test_cls_line_3( r )
 && _test_cls_segment_3( r )
 && _test_cls_sphere_3(r)
 && _test_cls_ray_3( r )
 && _test_cls_triangle_3( r )
 && _test_cls_tetrahedron_3( r )
 && _test_cls_aff_transformation_3( r )
 && _test_cls_iso_cuboid_3( r )
 ;
}

#endif // CGAL__TEST_3_C
