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
// source        : test_kernel_3.fw
// file          : _test_3.C
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_3_C
#define CGAL__TEST_3_C

#include <CGAL/_test_cls_vector_3.C>
#include <CGAL/_test_fct_vector_3.C>
#include <CGAL/_test_cls_point_3.C>
#include <CGAL/_test_fct_point_vector_3.C>
#include <CGAL/_test_fct_point_3.C>
#include <CGAL/_test_cls_direction_3.C>
#include <CGAL/_test_cls_plane_3.C>
#include <CGAL/_test_cls_line_3.C>
#include <CGAL/_test_cls_segment_3.C>
#include <CGAL/_test_cls_sphere_3.C>
#include <CGAL/_test_cls_ray_3.C>
#include <CGAL/_test_cls_triangle_3.C>
#include <CGAL/_test_cls_tetrahedron_3.C>
#include <CGAL/_test_cls_aff_transformation_3.C>


template <class R> bool _test_3(const R& r);


template <class R>
bool
_test_3(const R& r)
{
 return
    _test_cls_vector_3(r)
 && _test_fct_vector_3(r)
 && _test_cls_point_3(r)
 && _test_fct_point_vector_3(r)
 && _test_fct_point_3(r)
 && _test_cls_direction_3(r)
 && _test_cls_plane_3( r )
 && _test_cls_line_3( r )
 && _test_cls_segment_3( r )
 && _test_cls_sphere_3(r)
 && _test_cls_ray_3( r )
 && _test_cls_triangle_3( r )
 && _test_cls_tetrahedron_3( r )
 && _test_cls_aff_transformation_3( r )
 ;
}

#endif // CGAL__TEST_3_C
