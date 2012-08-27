// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL__TEST_2_C
#define CGAL__TEST_2_C

#include "_test_cls_vector_2.h"
#include "_test_fct_vector_2.h"
#include "_test_cls_point_2.h"
#include "_test_fct_point_vector_2.h"
#include "_test_fct_point_2.h"
#include "_test_fct_line_2.h"
#include "_test_fct_segment_2.h"
#include "_test_further_fct_point_2.h"
#include "_test_cls_direction_2.h"
#include "_test_fct_direction_2.h"
#include "_test_fct_point_line_2.h"
#include "_test_fct_point_segment_2.h"
#include "_test_further_fct_point_line_2.h"
#include "_test_cls_line_2.h"
#include "_test_cls_segment_2.h"
#include "_test_cls_ray_2.h"
#include "_test_cls_triangle_2.h"
#include "_test_cls_circle_2.h"
#include "_test_cls_iso_rectangle_2.h"
#include "_test_cls_aff_transformation_2.h"

template <class R>
bool
_test_2(const R& r)
{
 return
    _test_cls_vector_2(r)
 && _test_fct_vector_2(r)
 && _test_cls_point_2(r)
 && _test_fct_point_vector_2(r)
 && _test_fct_point_2(r)
 && _test_fct_line_2(r)
 && _test_fct_segment_2(r)
 && _test_further_fct_point_2(r)
 && _test_cls_direction_2(r)
 && _test_fct_direction_2(r)
 && _test_fct_point_line_2( r )
 && _test_fct_point_segment_2( r )
 && _test_further_fct_point_line_2( r )
 && _test_cls_line_2( r )
 && _test_cls_segment_2( r )
 && _test_cls_ray_2( r )
 && _test_cls_triangle_2( r )
 && _test_cls_circle_2( r )
 && _test_cls_iso_rectangle_2( r )
 && _test_cls_aff_transformation_2( r )
 ;
}
#endif // CGAL__TEST_2_C
