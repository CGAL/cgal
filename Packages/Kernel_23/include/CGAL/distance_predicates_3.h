// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : distance_predicates_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_DISTANCE_PREDICATES_3_H
#define CGAL_DISTANCE_PREDICATES_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#ifndef CGAL_DISTANCE_PREDICATESH3_H
#include <CGAL/distance_predicatesH3.h>
#endif // CGAL_DISTANCE_PREDICATESH3_H
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#ifndef CGAL_DISTANCE_PREDICATESC3_H
#include <CGAL/Cartesian/distance_predicates_3.h>
#endif // CGAL_DISTANCE_PREDICATESC3_H
#endif // CGAL_CARTESIAN_H

#ifdef CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/SimpleCartesian/distance_predicatesS3.h>
#endif // CGAL_SIMPLE_CARTESIAN_H


#ifndef CGAL_POINT_3_H
#include <CGAL/Point_3.h>
#endif // CGAL_POINT_3_H
#ifndef CGAL_PLANE_3_H
#include <CGAL/Plane_3.h>
#endif // CGAL_PLANE_3_H

#endif //CGAL_DISTANCE_PREDICATES_3_H
