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
// release       : $CGAL_Revision: CGAL-2.3-I-75 $
// release_date  : $CGAL_Date: 2001/06/21 $
//
// file          : include/CGAL/LEDA/system.h
// package       : cgal_window (1.0.3)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.3
// revision_date : 25 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#if defined(unix) || defined(__unix__) || defined(__unix) || defined(_AIX)
#include <CGAL/LEDA/sys/unix.h>

#elif defined(__CYGWIN32__)
#include <CGAL/LEDA/sys/cygwin32.h>

#elif defined(__WIN32__) || defined(_WIN32) || defined(__NT__)
#include <CGAL/LEDA/sys/win32.h>



#else
// unknown system

#endif

