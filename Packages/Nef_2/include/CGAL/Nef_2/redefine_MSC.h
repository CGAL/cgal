// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/redefine_MSC.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
// ============================================================================

#ifndef CGAL_REDEFINE_MSC_H
#define CGAL_REDEFINE_MSC_H
#if defined( _MSC_VER ) && ! defined(__INTEL_COMPILER)

#define RPolynomial_MSC RPMS
#define HalfEdgeDS_default_MSC HDSdMS
#define Nef_polyhedron_2_rep NP2R
#define Nef_polyhedron_2 NP2H
#define In_place_list_iterator Ipli
#define HalfedgeDS_in_place_list_vertex HDSv
#define Filtered_extended_homogeneous FEH
#define Extended_segment Eseg
#define Extended_point Epnt

#endif
#endif //CGAL_REDEFINE_MSC_H
