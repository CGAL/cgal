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
// file          : include/CGAL/LEDA/sys/cygwin32.h
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


#ifndef CGAL_WINDOW_SYS_CYGWIN32_H
#define CGAL_WINDOW_SYS_CYGWIN32_H

#define __win32__

#define __HAS_EXPLICIT_KEYWORD__
#define __HAS_TYPENAME_KEYWORD__
#define __HAS_MEMBER_TEMPLATES__


#define __explicit          explicit
#define __typename          typename
#define __temp_friend_decl  <>
#define __temp_func_inline


#define LITTLE_ENDIAN_MACHINE


//------------------------------------------------------------------------------
//  DLL definitions
//------------------------------------------------------------------------------

#define __exportC
#define __exportF
#define __exportD

#define _exportC
#define _exportF
#define _exportD



#endif
