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
// revision      : 2.1.2
// revision_date : 09 Jul 1998
// author(s)     : Kurt Mehlhorn
//                 Michael Seel
//                 Stefan Schirra
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

#ifndef CGAL_DD_GEO_CONFIG_H
#define CGAL_DD_GEO_CONFIG_H

// #ifndef __SUNPRO_CC              /* !!! should not be compiler-dependent !!! */
// #define DDGEO_STL_ITERATORS
// #endif // __SUNPRO_CC

#if ( __LEDA__ < 380 )
#define LEDA_PREFIXLI
#else
#define LEDA_PREFIXLI  LEDA::
#endif 

#if !defined(LEP_DDGEO_INCL_ID)
#define LEDA_ROOT_INCL_ID NOT_ANY_KERNEL_NUMBER
#define LEP_DDGEO_INCL_ID 21027
#include <LEDA/REDEFINE_NAMES.h>
#endif


#ifndef CGAL_CFG_NO_TEMPLATE_FRIEND_DISTINCTION
#define EGCSNTA <>
#else
#define EGCSNTA
#endif


#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
#define CHSIMPLEX  ch_simplex
#define DTVERTEX   dt_vertex
#else
#define CHSIMPLEX  ch_Simplex<CHTRAITS,POINT,PLANE>*
#define DTVERTEX   rc_Vertex<CHTRAITS,LPNT>*
#endif // CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS


#ifndef CGAL_CFG_NO_DEFAULT_TEMPLATE_ARGUMENTS
#define DDGEO_TEMPLATE_DEFAULTS
#endif // CGAL_CFG_NO_DEFAULT_TEMPLATE_ARGUMENTS



#if LEP_DDGEO_INCL_ID == 21027
#undef LEDA_ROOT_INCL_ID
#undef LEP_DDGEO_INCL_ID
#include <LEDA/UNDEFINE_NAMES.h>
#endif

#endif // CGAL_DD_GEO_CONFIG_H
