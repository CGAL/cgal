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
// file          : include/CGAL/LEDA/sys/win32.h
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0
// revision_date : 20 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================

#ifndef CGAL_WINDOW_SYS_WIN32_H
#define CGAL_WINDOW_SYS_WIN32_H

#if !defined(__win32__)
#define __win32__
#endif


#define LITTLE_ENDIAN_MACHINE


#if defined(__BORLANDC__)

// turn of some warnings
#pragma option -w-inl -w-ccc -w-aus -w-eff -w-par -w-rch

#define __HAS_EXPLICIT_KEYWORD__
#define __explicit explicit

#if __BORLANDC__ < 0x540
#define NO_STATIC_DATA_MEMBER
#define __ALL_MEMBERS_INSTANT__
#define __typename
#define __temp_friend_decl
#else
#define __HAS_TYPENAME_KEYWORD__
#define __typename typename
#define __temp_friend_decl <>
#define __HAS_MEMBER_TEMPLATES__
#endif

#define __temp_func_inline


#include <float.h>
/*
static int leda_init_fpu()
{ _control87(PC_53,MCW_PC);
  _control87(63U,MCW_EM);
  return 0;
}

static int setdouble_ieee_for_bcc = leda_init_fpu();
*/
#endif



#if defined(_MSC_VER)


// turn off some warnings

// exception handling turned off when using std headers
#pragma warning(disable:4530)

// no matching delete operator (in presence of -GX)
#pragma warning(disable:4291)

// missing dll-interface
#pragma warning(disable:4251)
#pragma warning(disable:4275)


#define __HAS_MEMBER_TEMPLATES__
#define __HAS_EXPLICIT_KEYWORD__
#define __explicit explicit
#define __temp_friend_decl
#define __temp_func_inline


#define __HAS_TYPENAME_KEYWORD__
#define __typename typename


#include <float.h>
/*
inline int leda_init_fpu() 
{ _control87(_PC_53,_MCW_PC);
  return 0; 
}

static int setdouble_ieee_for_msc = leda_init_fpu();
*/
#endif



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
