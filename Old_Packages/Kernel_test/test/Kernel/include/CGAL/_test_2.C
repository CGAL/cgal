// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : _test_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_2_C
#define CGAL__TEST_2_C

#include <CGAL/_test_cls_vector_2.C>
#include <CGAL/_test_fct_vector_2.C>
#include <CGAL/_test_cls_point_2.C>
#include <CGAL/_test_fct_point_vector_2.C>
#include <CGAL/_test_fct_point_2.C>
#include <CGAL/_test_further_fct_point_2.C>
#include <CGAL/_test_cls_direction_2.C>
#include <CGAL/_test_fct_direction_2.C>
#include <CGAL/_test_fct_point_line_2.C>
#include <CGAL/_test_further_fct_point_line_2.C>
#include <CGAL/_test_cls_line_2.C>
#include <CGAL/_test_cls_segment_2.C>
#include <CGAL/_test_cls_ray_2.C>
#include <CGAL/_test_cls_triangle_2.C>
#include <CGAL/_test_cls_circle_2.C>
#include <CGAL/_test_cls_iso_rectangle_2.C>
#include <CGAL/_test_cls_aff_transformation_2.C>


template <class R> bool _test_2(const R& r);


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
 && _test_further_fct_point_2(r)
 && _test_cls_direction_2(r)
 && _test_fct_direction_2(r)
 && _test_fct_point_line_2( r )
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
