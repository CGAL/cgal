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
// source        :
// file          : _test_fct_points_implicit_sphere.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_FCT_POINTS_IMPLICIT_SPHERE_H
#define CGAL__TEST_FCT_POINTS_IMPLICIT_SPHERE_H

#include <CGAL/Point_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/rational_rotation.h>
#include <CGAL/predicates_on_points_3.h>
#include <CGAL/squared_distance_3.h>


template <class R> bool _test_fct_points_implicit_sphere(const R&);


#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/_test_fct_points_implicit_sphere.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL__TEST_FCT_POINTS_IMPLICIT_SPHERE_H
