// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/IO/cgal_logo.h
// package       : window 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken  
// ======================================================================

#ifndef CGAL_LOGO_H
#define CGAL_LOGO_H

#if !defined(__LEDA__) || (__LEDA__ > 400)
extern const char * cgal_logo[];
#else
extern char * cgal_logo[];
#endif // __LEDA__

#endif // CGAL_LOGO_XPM_H
