// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 21
//
// file          : 
// package       : Convex_hull_3 (2.6)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : LEP dd_geo_kernel
// revision      : 2.1
// revision_date : 
// author(s)     : Kurt Mehlhorn
//                 Michael Seel
//
// coordinator   : MPI, Saarbruecken <Stefan.Schirra@mpi-sb.mpg.de>
// ======================================================================
/*******************************************************************************
+
+  LEP dd_geokernel 2.1
+
+  This file is part of the research version of a LEDA extension package,
+  that can be used free of charge in academic research and teaching. 
+  Any commercial use of this software requires a commercial license,
+  which is distributed by the Algorithmic Solutions GmbH, 
+  Postfach 151101, 66041 Saarbruecken, FRG (fax +49 681 31104).
+
+  Copyright (c) 1997-1998  by  Max-Planck-Institut fuer Informatik
+  Im Stadtwald, 66123 Saarbruecken, Germany     
+  All rights reserved.
+ 
*******************************************************************************/

//#ifndef LEDA_DEBUG_H
//#define LEDA_DEBUG_H

#include <LEDA/stream.h>

#undef TRACE
#undef TRACEN
#undef TRACEV
#undef CTRACE
#undef CTRACEN
#undef ASSERT

#ifdef _DEBUG
#define TRACE(t)   cerr << " " << t  
#else
#define TRACE(t) 
#endif

#ifdef _DEBUG
#define TRACEV(t)   cout << " " << #t << " = " << (t)  << endl
#else
#define TRACEV(t) 
#endif

#ifdef _DEBUG
#define TRACEN(t)   cerr << " " << t << "\n"
#else
#define TRACEN(t) 
#endif

#ifdef _DEBUG
#define CTRACE(b,t)  if(b) cerr << " " << t; else cerr << " 0"
#else
#define CTRACE(b,t) 
#endif

#ifdef _DEBUG
#define CTRACEN(b,t)  if(b) cerr << " " << t << "\n"; else cerr << " 0\n"
#else
#define CTRACEN(b,t) 
#endif

#ifndef _ASSERT
#define  ASSERT(cond,fstr) 
#else
#define ASSERT(cond,fstr)   \
  if (!(cond)) {       \
    cerr << "   ASSERT:   " << #fstr << endl; \
    cerr << "   COND:     " << #cond << endl; \
    cerr << "   POSITION: " << __FILE__ << " at line " << __LINE__ << endl; \
    abort();           \
  }
#endif


//#endif //LEDA_DEBUG_H


