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
// file          : include/CGAL/LEDA/sys/unix.h
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0
// revision_date : 20 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


#ifndef CGAL_WINDOW_SYS_UNIX_H
#define CGAL_WINDOW_SYS_UNIX_H


//----------------------------------------------------------------------------
// system/compiler dependent definitions 
//----------------------------------------------------------------------------

// In the first section of this file some flags are defined indicating
// that certain c++ features are supported or not by particular compilers. 
// If your compiler does not appear you may have to add it.
// 
// __HAS_MEMBER_TEMPLATES__   : member function templates
// __HAS_TYPENAME_KEYWORD__   : has keyword "typename"
// __HAS_EXPLICIT_KEYWORD__   : has keyword "explicit"
// __HAS_BUILTIN_BOOL__       : has built-in boolean
// __HAS_FNMATCH_H__          : has <fnmatch.h>
// __USE_VALUES_H__           : use <values.h> header file 
// __NO_TEMPLATE_ARG_BASE__   : template arguments may not be used as base class
// __NO_FUNCTEMP_DEF_ARGS__   : function templates may not have default args
// __NO_EXPLICIT_DESTRUCTION__: no explicit destructor call
// __NO_CAST_TO_LOCAL_TYPE__  : no casting to a local type of a class
// __ELSE_SCOPE_BUG__         : else-part does not start new scope
// __ALL_MEMBERS_INSTANT__    : instantiates always ALL members of a template


#if !defined(__unix__) 
#define __unix__
#endif


#if defined(__SVR4) || defined(_SYSTYPE_SVR4) || defined(__SYSTYPE_SVR4)
#if !defined(__svr4__)
#define __svr4__
#endif
#endif



#if defined(mips) && !defined(__GNUC__)
#if defined(_COMPILER_VERSION)
#define __mipspro__
#else
#define __ELSE_SCOPE_BUG__
#define __NO_EXPLICIT_DESTRUCTION__
#endif
#endif


#if defined(__GNUC__) || defined(__mipspro__) || defined(__KCC)\
 || __SUNPRO_CC >= 0x500
#define __HAS_TYPENAME_KEYWORD__
#define __HAS_EXPLICIT_KEYWORD__
#define __typename         typename
#define __explicit         explicit
#else
#define __typename
#define __explicit
#endif


#if defined(__GNUC__) || defined(__mipspro__) || defined(__KCC)
// || __SUNPRO_CC >= 0x510
#define __HAS_MEMBER_TEMPLATES__
#endif


 
#if defined(__GNUC__) || defined(__KCC)\
                      || (defined(__mipspro__) && _COMPILER_VERSION >= 730)
#define __temp_friend_decl <>
#else
#define __temp_friend_decl
#endif


#if defined(__DECCXX) || defined(__KCC) || defined(__hpuxcc__)\
                                        || defined(__SUNPRO_CC) 
#define __temp_func_inline inline
#else
#define __temp_func_inline
#endif



#if defined(__svr4__) || defined(__linux__) || defined(__hpux)
#define __HAS_FNMATCH_H__
#endif



#if defined(__alpha) || defined(__sparcv9) || _MIPS_SZLONG == 64
#define WORD_LENGTH_64
#endif


#if !defined(__GNUC__) && !defined(__hpux)
#define __USE_VALUES_H__
#endif


#if defined(__ALL_MEMBERS_INSTANT__) || defined(__NO_EXPLICIT_DESTRUCTION__)
#define __NO_TEMPLATE_ALGS__
#endif


#if defined(i386) || defined(__i386)
#define LITTLE_ENDIAN_MACHINE
#endif

#if defined(__linux__) && (defined(i386) || defined(__i386))

/*
inline int leda_init_fpu()
{ // set internal rounding of FPU to double (not extended) format
  // and leave rounding to nearest as well as exceptions unchanged
  int cw = 4735;
  __asm__("fldcw %0" : : "m" (cw));
  return 0; 
}

static int set_double_ieee_for_linux = leda_init_fpu();
*/

#endif


#define __exportC
#define __exportF
#define __exportD

#define _exportC
#define _exportF
#define _exportD

#endif
