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
// file          : _test_orientation_and_bounded_side.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_ORIENTATION_AND_BOUNDED_SIDE_H
#define CGAL__TEST_ORIENTATION_AND_BOUNDED_SIDE_H

#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>


template <class R>
bool
_test_orientation_and_bounded_side(const R&);


#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/_test_orientation_and_bounded_side.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL__TEST_ORIENTATION_AND_BOUNDED_SIDE_H
